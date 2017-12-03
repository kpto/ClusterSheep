# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 06:28:56 2017

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
import multiprocessing as mp
from threading import Thread
import time
from pathlib import Path
from uuid import uuid4
import subprocess
import pickle
import sys
import ctypes

from envr.session import get_session
from prcs.parallel.logging_setup import logging_setup
import reader.mzxml as mzxml
import reader.mzml as mzml
from property import *
from share.misc import _clean_temporary
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    ignore_errors = None
    num_of_threads = None
    true_precursor_mass = None
    input_limit = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====

# ====BEGIN OF CODE====
def build_index():
    _refresh_session()
    if len(session.ms_exp_files) == 0:
        err_msg = '\nNo ms exp files specified.'
        logging.error(err_msg)
        raise Exception(err_msg)

    mzxml._refresh_session()
    mzml._refresh_session()
    processes = []
    q = mp.Queue()
    log_lock = mp.Lock()
    merge_lock = mp.Lock()
    exit_signal = mp.Value(ctypes.c_bool, False)
    valid_ms2_count = mp.Value(ctypes.c_uint64, 0)
    total_ms2_count = mp.Value(ctypes.c_uint64, 0)
    finish_count = mp.Value(ctypes.c_uint64, 0)
    num_files = len(session.ms_exp_files)

    temp_storage = Path.cwd().joinpath(str(uuid4()))
    temp_storage.mkdir()
    temp_storage = temp_storage.resolve()
    Path(temp_storage, 'prm').touch()
    Path(temp_storage, 'chg').touch()
    Path(temp_storage, 'fli').touch()
    Path(temp_storage, 'off').touch()
    Path(temp_storage, 'iid').touch()
    Path(temp_storage, 'nvi').touch()

    for n, path in enumerate(session.ms_exp_files): q.put((n, path))
    for i in range(num_of_threads): q.put(None)
    for pid in range(num_of_threads):
        processes.append(mp.Process(target=_worker, args=(pid, q, temp_storage, valid_ms2_count, total_ms2_count,
                                                          finish_count, log_lock, merge_lock, exit_signal)))
    reporter = Thread(target=_reporter, args=(finish_count, num_files, log_lock, exit_signal))

    logging.info('......Start building index with {} subprocesses......'.format(num_of_threads))

    try:
        for i in range(num_of_threads):
            processes[i].start()
            time.sleep(0.1)
        reporter.start()
        for i in range(num_of_threads): processes[i].join()
        exit_signal.value = True
        reporter.join()
    except (Exception, KeyboardInterrupt) as e:
        if type(e) is KeyboardInterrupt:
            with log_lock:
                logging.info('Received KeyboardInterrupt, index building canceled.')
        exit_signal.value = True
        with log_lock:
            logging.debug('Waiting processes to exit.')
        reporter.join()
        for i in range(num_of_threads):
            if processes[i].is_alive():
                processes[i].join()
        logging.debug('All processes exited.')
        _error_clean(temp_storage, q)
        raise

    # in case that subprocesses were exited with errors
    if 1 in [p.exitcode for p in processes]:
        err_msg = '\nSubprocesses exited with abnormal exitcode.'
        logging.error(err_msg)
        logging.debug('Emptying queue')
        _error_clean(temp_storage, q)
        raise mp.ProcessError(err_msg)

    if valid_ms2_count.value == 0:
        err_msg = '\nNo valid MS2 scan has been found.'
        logging.error(err_msg)
        _error_clean(temp_storage, q)
        raise Exception(err_msg)

    logging.info('......Finish building index, found {} valid MS2 scan (out of {} in total)......'
                 .format(valid_ms2_count.value, total_ms2_count.value))
    try:
        _integrate(temp_storage, valid_ms2_count.value)
    except Exception:
        err_msg = '\nFailed to integrate index data.'
        logging.error(err_msg)
        raise
    _clean_temporary(temp_storage)

    session.config.ii_finished.value = True
    logging.debug('Mounting internal index file')
    session.mount_internal_index()
    return


def _reporter(finish_count, num_files, log_lock, exit_signal):
    progress = 0
    while True:
        try:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Reporter thread exits now.')
                break
            time.sleep(0.2)
            old_progress = progress
            progress = int(finish_count.value / num_files * 100)
            if progress > old_progress:
                with log_lock:
                    logging.info('Building progress: {}% ({} out of {})' \
                                 .format(progress, finish_count.value, num_files))
        except KeyboardInterrupt:
            with log_lock:
                logging.debug('Reporter thread received KeyboardInterrupt, exits now.')
            break
    return


def _error_clean(temp_storage, q):
    while q.qsize() > 0:
        q.get()
    _clean_temporary(temp_storage)
    return


def _integrate(temp_storage, length):
    prm_path = Path(temp_storage, 'prm')
    chg_path = Path(temp_storage, 'chg')
    fli_path = Path(temp_storage, 'fli')
    off_path = Path(temp_storage, 'off')
    iid_path = Path(temp_storage, 'iid')
    nvi_path = Path(temp_storage, 'nvi')

    prm = np.memmap(str(prm_path), dtype=II_PRECURSOR_MASS_DATA_TYPE, mode='r+', shape=length)
    chg = np.memmap(str(chg_path), dtype=II_PRECURSOR_CHARGE_DATA_TYPE, mode='r+', shape=length)
    fli = np.memmap(str(fli_path), dtype=II_FILE_ID_DATA_TYPE, mode='r+', shape=length)
    off = np.memmap(str(off_path), dtype=II_OFFSET_DATA_TYPE, mode='r+', shape=(length, 2))
    iid = np.memmap(str(iid_path), dtype=II_INTERNAL_ID_DATA_TYPE, mode='r+', shape=length)
    nvi = np.memmap(str(nvi_path), dtype=II_NATIVE_ID_DATA_TYPE, mode='r+', shape=length)

    logging.debug('Sorting index by precursor mass')
    sorted_idx = np.argsort(prm)
    prm = prm[sorted_idx]
    chg = chg[sorted_idx]
    fli = fli[sorted_idx]
    off = off[sorted_idx]
    iid = iid[sorted_idx]
    nvi = nvi[sorted_idx]

    logging.debug('Integrating index data')
    data_order = (DATA_TYPE_SIZE[II_PRECURSOR_MASS_DATA_TYPE],
                  DATA_TYPE_SIZE[II_PRECURSOR_CHARGE_DATA_TYPE],
                  DATA_TYPE_SIZE[II_FILE_ID_DATA_TYPE],
                  DATA_TYPE_SIZE[II_OFFSET_DATA_TYPE] * 2,
                  DATA_TYPE_SIZE[II_INTERNAL_ID_DATA_TYPE],
                  DATA_TYPE_SIZE[II_NATIVE_ID_DATA_TYPE])
    offsets = [FILE_HEADER_SIZE]
    for n in range(len(data_order) - 1):
        offsets.append(sum(data_order[:n+1])*length + FILE_HEADER_SIZE)

    total_size = sum(data_order) * length + FILE_HEADER_SIZE
    file_path = Path.cwd().joinpath(session.name + FILE_EXTENSION_INTERNAL_INDEX)
    if file_path.exists():
        file_path.unlink()
    file_path.touch()
    try:
        subprocess.check_output(('fallocate', '-l', str(total_size), str(file_path)), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        err_msg = '\nFailed to allocate space for file integration.'\
                  '\nstderr: {}'.format(e.output)
        logging.error(err_msg)
        raise

    metadata = {'version': VERSION,
                'file_type': HEADER_LABEL_INTERNAL_INDEX,
                'magic_label': uuid4(),
                'session_label': session.magic_label,
                'length': length,
                'true_precursor_mass': true_precursor_mass,
                'precursor_mass_data_type': II_PRECURSOR_MASS_DATA_TYPE,
                'precursor_charge_data_type': II_PRECURSOR_CHARGE_DATA_TYPE,
                'file_id_data_type': II_FILE_ID_DATA_TYPE,
                'offset_data_type': II_OFFSET_DATA_TYPE,
                'internal_id_data_type': II_INTERNAL_ID_DATA_TYPE,
                'native_id_data_type': II_NATIVE_ID_DATA_TYPE,
                'precursor_mass_offset': offsets[0],
                'precursor_charge_offset': offsets[1],
                'file_id_offset': offsets[2],
                'offset_offset': offsets[3],
                'internal_id_offset': offsets[4],
                'native_id_offset': offsets[5]}
    header = pickle.dumps(metadata) + b'%end%'
    if len(header) > FILE_HEADER_SIZE:
        err_msg = '\nHeader having size larger than the reserved area.'
        logging.error(err_msg)
        raise Exception(err_msg)
    with file_path.open('r+b') as fp:
        fp.write(header)

    temp = np.memmap(str(file_path), offset=offsets[0], dtype=II_PRECURSOR_MASS_DATA_TYPE, mode='r+', shape=length)
    temp[:] = prm
    del temp
    temp = np.memmap(str(file_path), offset=offsets[1], dtype=II_PRECURSOR_CHARGE_DATA_TYPE, mode='r+', shape=length)
    temp[:] = chg
    del temp
    temp = np.memmap(str(file_path), offset=offsets[2], dtype=II_FILE_ID_DATA_TYPE, mode='r+', shape=length)
    temp[:] = fli
    del temp
    temp = np.memmap(str(file_path), offset=offsets[3], dtype=II_OFFSET_DATA_TYPE, mode='r+', shape=(length, 2))
    temp[:] = off
    del temp
    temp = np.memmap(str(file_path), offset=offsets[4], dtype=II_INTERNAL_ID_DATA_TYPE, mode='r+', shape=length)
    temp[:] = iid
    del temp
    temp = np.memmap(str(file_path), offset=offsets[5], dtype=II_NATIVE_ID_DATA_TYPE, mode='r+', shape=length)
    temp[:] = nvi
    del temp

    return


def _worker(pid, q, temp_storage, valid_ms2_count, total_ms2_count, finish_count, log_lock, merge_lock, exit_signal):
    try:
        logging_setup()
        with log_lock:
            logging.debug('Index building subprocess {} started.'.format(pid))

        prm_path = Path(temp_storage, 'prm')
        chg_path = Path(temp_storage, 'chg')
        fli_path = Path(temp_storage, 'fli')
        off_path = Path(temp_storage, 'off')
        iid_path = Path(temp_storage, 'iid')
        nvi_path = Path(temp_storage, 'nvi')

        while True:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Subprocess {}: Received exit signal, exits now.'.format(pid))
                break
            if input_limit > 0 and valid_ms2_count.value >= input_limit:
                with log_lock:
                    logging.debug('Subprocess {}: Exceeded input limit, purges queue and exits now.'.format(pid))
                item = 0
                while item is not None:
                    item = q.get()
                break
            item = q.get()
            if item:
                file_id, path = item
                try:
                    index, valid, total = mzxml.build_index(path, log_lock=log_lock)
                except Exception:
                    if ignore_errors:
                        wrn_msg = 'Subprocess {}: Failed to parse MS exp file: {}'.format(pid, path)
                        with log_lock:
                            logging.warning(wrn_msg)
                        continue
                    else:
                        raise
                if valid != 0:
                    index.file_id = np.full(len(index.file_id), file_id, dtype=II_FILE_ID_DATA_TYPE)
                    with merge_lock:
                        index.internal_id += valid_ms2_count.value
                        valid_ms2_count.value += valid
                        total_ms2_count.value += total
                        prm = np.memmap(str(prm_path), dtype=II_PRECURSOR_MASS_DATA_TYPE, mode='r+', shape=valid_ms2_count.value)
                        chg = np.memmap(str(chg_path), dtype=II_PRECURSOR_CHARGE_DATA_TYPE, mode='r+', shape=valid_ms2_count.value)
                        fli = np.memmap(str(fli_path), dtype=II_FILE_ID_DATA_TYPE, mode='r+', shape=valid_ms2_count.value)
                        off = np.memmap(str(off_path), dtype=II_OFFSET_DATA_TYPE, mode='r+', shape=(valid_ms2_count.value, 2))
                        iid = np.memmap(str(iid_path), dtype=II_INTERNAL_ID_DATA_TYPE, mode='r+', shape=valid_ms2_count.value)
                        nvi = np.memmap(str(nvi_path), dtype=II_NATIVE_ID_DATA_TYPE, mode='r+', shape=valid_ms2_count.value)

                        prm[-valid:] = index.precursor_mass
                        chg[-valid:] = index.precursor_charge
                        fli[-valid:] = index.file_id
                        off[-valid:] = index.offset
                        iid[-valid:] = index.internal_id
                        nvi[-valid:] = index.native_id
                else:
                    wrn_msg = 'Subprocess {}: Found 0 valid MS2 scan in file: {}'.format(pid, path)
                    with log_lock:
                        logging.warning(wrn_msg)

                finish_count.value += 1
            else:
                logging.debug('Subprocess {}: Got None object, work done.'.format(pid))
                break

    except (Exception, KeyboardInterrupt) as e:
        if type(e) is KeyboardInterrupt:
            with log_lock:
                logging.debug('Subprocess {}: Received KeyboardInterrupt, exits now.'.format(pid))
        else:
            with log_lock:
                logging.exception('\nSubprocess {}: Ended unexpectedly. Logging traceback:\n'
                                  '==========TRACEBACK==========\n'.format(pid))
            exit_signal.value = True
            sys.exit(1)
    return


def _refresh_session():
    global ignore_errors, num_of_threads, true_precursor_mass, input_limit
    ignore_errors = session.flags.ignore_errors
    num_of_threads = session.config.ii_num_of_threads.value
    true_precursor_mass = session.config.ii_true_precursor_mass.value
    input_limit = session.config.ii_input_limit.value
    return
# ====END OF CODE====

