# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 07:09:03 2017

Author: Paul TO, Ka Po
Contact: kpto@ust.hk
Institude: Hong Kong University of Science and Technology
Department: Division of Biomedical Engineering

Supervisor: Prof. Henry LAM, H. N.
Contact: kehlam@ust.hk
Institude: Hong Kong University of Science and Technology
Department: Department of Chemical and Biomolecular Engineering

Desciption of this module:
"""


# ====BEGIN OF MODULE IMPORT====
import logging
from threading import Lock
import re
import zlib
import base64

from ClusterSheep.envr.session import get_session
from ClusterSheep.share.internal_index import Index
from ClusterSheep.property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    ignore_errors = None
    min_num_peaks = None
    true_precursor_mass = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def build_index(file, log_lock=Lock()):
    total_ms2_count = 0
    valid_ms2_count = 0
    fail_count = 0

    # filtering
    with file.open('rb') as fp:
        offset_list = _get_offset_list(fp, log_lock=log_lock)
        scan_count = _get_scan_count(fp, log_lock=log_lock)

        precursor_mass_arr = np.empty(shape=scan_count, dtype=II_PRECURSOR_MASS_DATA_TYPE)
        precursor_charge_arr = np.empty(shape=scan_count, dtype=II_PRECURSOR_CHARGE_DATA_TYPE)
        file_id_arr = np.zeros(shape=scan_count, dtype=II_FILE_ID_DATA_TYPE)
        offset_arr = np.empty(shape=(scan_count, 2), dtype=II_OFFSET_DATA_TYPE)
        internal_id_arr = np.empty(shape=scan_count, dtype=II_INTERNAL_ID_DATA_TYPE)
        native_id_arr = np.empty(shape=scan_count, dtype=II_NATIVE_ID_DATA_TYPE)

        for i in range(len(offset_list)-1):
            start_pos = offset_list[i]
            end_pos = offset_list[i+1]
            fp.seek(start_pos)
            block = fp.read(end_pos-start_pos)

            scan_num = None
            ms_level = None
            peaks_count = None
            precursor_mz = None
            precursor_charge = None

            try:
                ms_level = int(re.search(b'msLevel="([^"]+)', block).groups()[0])
                if ms_level != 2:
                    continue
                scan_num = int(re.search(b'num="([^"]+)', block).groups()[0])
                peaks_count = int(re.search(b'peaksCount="([^"]+)', block).groups()[0])
                precursor_mz = float(re.search(b'>([^>]+)</precursorMz>', block).groups()[0])
                precursor_charge = re.search(b'precursorCharge="([^"]+)', block)
                precursor_charge = int(precursor_charge.groups()[0]) if precursor_charge else 1
            except Exception:
                if ignore_errors:
                    fail_count += 1
                    continue
                else:
                    err_msg = '\nUnable to get the necessary information of the scan located at {}.'\
                              '\nScan num: {}    MS level: {}    Peaks count: {}    Pre m/z: {}    Pre charge: {}'\
                              .format(start_pos, scan_num, ms_level, peaks_count, precursor_mz, precursor_charge)
                    with log_lock:
                        logging.error(err_msg)
                    raise

            if peaks_count >= min_num_peaks:
                if true_precursor_mass:
                    precursor_mass = precursor_mz * precursor_charge
                else:
                    precursor_mass = precursor_mz

                precursor_mass_arr[valid_ms2_count] = precursor_mass
                precursor_charge_arr[valid_ms2_count] = precursor_charge
                offset_arr[valid_ms2_count] = (start_pos, end_pos)
                internal_id_arr[valid_ms2_count] = valid_ms2_count
                native_id_arr[valid_ms2_count] = scan_num

                valid_ms2_count += 1
            total_ms2_count += 1

    precursor_mass_arr = precursor_mass_arr[:valid_ms2_count]
    precursor_charge_arr = precursor_charge_arr[:valid_ms2_count]
    file_id_arr = file_id_arr[:valid_ms2_count]
    offset_arr = offset_arr[:valid_ms2_count]
    internal_id_arr = internal_id_arr[:valid_ms2_count]
    native_id_arr = native_id_arr[:valid_ms2_count]

    if fail_count > 0:
        wrn_msg = '{} MS scan were skipped while reading file \"{}\" due to reading errors.'\
                  .format(fail_count, file)
        with log_lock:
            logging.warning(wrn_msg)

    with log_lock:
        logging.debug('{} valid MS2 scans (out of {}) were indexed from file: {}'\
                      .format(valid_ms2_count, total_ms2_count, file))

    new_index = Index()
    new_index.precursor_mass = precursor_mass_arr
    new_index.precursor_charge = precursor_charge_arr
    new_index.file_id = file_id_arr
    new_index.offset = offset_arr
    new_index.internal_id = internal_id_arr
    new_index.native_id = native_id_arr
    new_index.length = valid_ms2_count
    return new_index, valid_ms2_count, total_ms2_count


def _get_offset_list(fp, log_lock=Lock()):
    index_offset = _get_index_offset(fp)
    fail_count = 0
    offset_list = []
    fp.seek(index_offset)
    for l in fp.readlines():
        try:
            temp = re.search(b'>([^>]+)</offset>', l)
            if temp:
                offset_list.append(int(temp.groups()[0]))
        except Exception:
            if ignore_errors:
                fail_count += 1
                continue
            else:
                err_msg = '\nUnable to extract the ms scan offset value from line: {}'.format(l)
                with log_lock:
                    logging.error(err_msg)
            raise
    if len(offset_list) == 0:
        wrn_msg = 'Offset list is empty. File \"{}\" contains no MS scan.'.format(fp.name)
        with log_lock:
            logging.warning(wrn_msg)
    with log_lock:
        logging.debug('Extracted offset list with length of {} from file: {}'.format(len(offset_list), fp.name))
    if fail_count > 0:
        wrn_msg = '{} MS scan offsets were skipped while reading file \"{}\" due to reading errors.'\
                  .format(fail_count, fp.name)
        with log_lock:
            logging.warning(wrn_msg)
    # append the ending offset (nearly) of the last spectrum block
    offset_list.append(index_offset)
    return offset_list


def _get_scan_count(fp, log_lock=Lock()):
    end = fp.seek(0, 2)
    fp.seek(0)
    while fp.tell() < end:
        line = fp.readline()
        try:
            temp = re.search(b'scanCount="([^"]+)', line)
            if temp:
                scan_count = int(temp.groups()[0])
                break
        except Exception:
            if not ignore_errors:
                err_msg = '\nUnable to extract the scan count value from line: {}'.format(line)
                with log_lock:
                    logging.error(err_msg)
            raise
    else:
        if not ignore_errors:
            err_msg = '\nCould not find keyword \"scanCount\".'
            with log_lock:
                logging.error(err_msg)
        raise SyntaxError(err_msg)
    return scan_count


# find the offset of the index, search the last 500 bytes only
def _get_index_offset(fp, log_lock=Lock()):
    step = II_INDEX_OFFSET_SEARCH_AREA
    cursor = fp.seek(0, 2) - step
    fp.seek(cursor)
    block = fp.read()
    try:
        index_offset = int(re.search(b'<indexOffset>(.+)</indexOffset>', block).groups()[0])
    except Exception:
        if not ignore_errors:
            err_msg = '\nUnable to extract index offset.'\
                      '\nKeyword searching returned {}'.format(re.search(b'<indexOffset>(.+)</indexOffset>', block))
            with log_lock:
                logging.error(err_msg)
        raise
    return index_offset


def get_raw_text(file, offset):
    with file.open('rb') as fp:
        start_pos, end_pos = offset
        fp.seek(start_pos)
        block = fp.read(end_pos - start_pos)
        try:
            block = block[:re.search(b'</scan>', block).span()[1]]
        except Exception:
            wrr_msg = 'Unable to find the ending tag "</scan>". It may be a broken file.'
            logging.warning(wrr_msg)
    return block


def _get_peaks_raw(file, offset, log_lock=Lock()):
    with file.open('rb') as fp:
        start_pos, end_pos = offset
        fp.seek(start_pos)
        block = fp.read(end_pos - start_pos)
        try:
            start_pos = re.search(b'<peaks', block).span()[0]
            end_pos = re.search(b'</peaks>', block).span()[1]
            block = block[start_pos:end_pos]
        except Exception:
            if not ignore_errors:
                err_msg = '\nUnable to find the peak element. It may be a broken file.'
                with log_lock:
                    logging.error(err_msg)
            raise

        precision = None
        compression = None

        try:
            precision = re.search(b'precision="([^"]+)', block).groups()[0].decode()
            compression = re.search(b'compressionType="([^"]+)', block)
            if compression:
                compression = compression.groups()[0].decode()
        except Exception:
            if not ignore_errors:
                err_msg = '\nUnable to get the necessary information.'\
                          '\nPrecision: {}    Compression: {}'.format(precision, compression)
                with log_lock:
                    logging.error(err_msg)
            raise

        peaks_raw = re.search(b'>([^>]+)</peaks>', block).groups()[0]

    return peaks_raw, precision, compression


def get_peaks(file, offset, log_lock=Lock()):
    peaks_raw, precision, compression = _get_peaks_raw(file, offset, log_lock)
    try:
        binary = base64.decodebytes(peaks_raw)
        if compression == 'zlib':
            binary = zlib.decompress(binary)
        elif compression == 'none' or compression is None:
            pass
        else:
            if not ignore_errors:
                err_msg = '\nInvalid compression type: {}'.format(compression)
                with log_lock:
                    logging.error(err_msg)
            raise TypeError(err_msg)
        if precision == '32': dtype = np.dtype('f4')
        elif precision == '64': dtype = np.dtype('f8')
        else:
            if not ignore_errors:
                err_msg = '\nInvalid precision: {}'.format(precision)
                with log_lock:
                    logging.error(err_msg)
            raise TypeError(err_msg)
        if len(binary) % dtype.itemsize != 0:
            if not ignore_errors:
                err_msg = '\nThe length of bytes {} does not match the precision {}'\
                          .format(len(binary), precision)
                with log_lock:
                    logging.error(err_msg)
        temp = np.frombuffer(binary, dtype=dtype)
        temp = temp.byteswap()
        intensity, mz = np.rot90(temp.reshape((int(temp.size/2), 2)))
    except Exception:
        if not ignore_errors:
            err_msg = '\nUnable to decode the text.'\
                      '\nRaw text:'\
                      '\n{}'.format(peaks_raw)
            with log_lock:
                logging.error(err_msg)
        raise
    return mz, intensity


def _refresh_session():
    global ignore_errors, min_num_peaks, true_precursor_mass
    ignore_errors = session.flags.ignore_errors
    min_num_peaks = session.config.ii_min_num_peaks.value
    true_precursor_mass = session.config.ii_true_precursor_mass.value
    return


_refresh_session()
# ====END OF CODE====

