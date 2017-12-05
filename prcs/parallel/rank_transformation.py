# -*- coding: utf-8 -*-
"""
Created on 13:36:24 10/09/2017

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
from threading import Thread, Lock
import time
from pathlib import Path
from uuid import uuid4
import subprocess
import pickle
import sys
import ctypes

from envr.session import get_session
from prcs.parallel.logging_setup import logging_setup
from prcs.parallel.binning import _binning
from property import *
import reader.mzxml as mzxml
import reader.mzml as mzml
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    ignore_errors = None
    true_precursor_mass = None
    num_of_peaks = None
    precursor_removal_range = None
    mz_range = None
    bins_per_th = None
    num_of_threads = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
class _Dispatcher:

    def __init__(self, start, end, step):
        self.lock = mp.Lock()
        self.current = mp.Value(ctypes.c_uint64, 0)
        self.start = start
        self.end = end
        self.step = step
        return

    def __iter__(self):
        return self

    def __next__(self):
        with self.lock:
            if self.current.value == self.end:
                raise StopIteration
            start = self.current.value
            end = start + self.step
            if end > self.end: end = self.end
            self.current.value = end
        return start, end

    def now(self):
        return self.current.value
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def rank_transform():
    _refresh_session()
    processes = []
    index_length = len(session.internal_index)
    dispatcher = _Dispatcher(0, index_length, 10000)
    log_lock = mp.Lock()
    merge_lock = mp.Lock()
    exit_signal = mp.Value(ctypes.c_bool, False)

    file_path = Path.cwd().joinpath(session.name + FILE_EXTENSION_RANKED_SPECTRA)
    if file_path.exists(): file_path.unlink()

    data_order = (DATA_TYPE_SIZE[RT_MZ_DATA_TYPE] * num_of_peaks,
                  DATA_TYPE_SIZE[RT_INTENSITY_DATA_TYPE] * num_of_peaks)
    offsets = [FILE_HEADER_SIZE]
    for n in range(len(data_order) - 1):
        offsets.append(sum(data_order[:n+1])*index_length + FILE_HEADER_SIZE)

    total_size = sum(data_order) * index_length + FILE_HEADER_SIZE

    try:
        subprocess.check_output(('fallocate', '-l', str(total_size), str(file_path)), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        err_msg = '\nFailed to allocate space for ranked spectra storage.'\
                  '\nstderr: {}'.format(e.output)
        logging.error(err_msg)
        raise

    for pid in range(num_of_threads):
        processes.append(mp.Process(target=_worker, args=(pid, dispatcher, offsets, log_lock, merge_lock, exit_signal)))
    reporter = Thread(target=_reporter, args=(dispatcher, index_length, log_lock, exit_signal))

    logging.info('......Start rank transformation with {} subprocesses......'.format(num_of_threads))

    try:
        for i in range(num_of_threads):
            processes[i].start()
            time.sleep(0.1)
        reporter.start()
        for i in range(num_of_threads): processes[i].join()
        exit_signal.value = True
        if reporter.is_alive():
            reporter.join()
    except (Exception, KeyboardInterrupt) as e:
        if type(e) is KeyboardInterrupt:
            with log_lock:
                logging.info('Received KeyboardInterrupt, index building canceled.')
        exit_signal.value = True
        with log_lock:
            logging.debug('Waiting processes to exit.')
        if reporter.is_alive():
            reporter.join()
        for i in range(num_of_threads):
            if processes[i].is_alive():
                processes[i].join()
        logging.debug('All processes exited.')
        raise

    # in case that subprocesses were exited with errors
    if 1 in [p.exitcode for p in processes]:
        err_msg = '\nSubprocesses exited with abnormal exitcode.'
        logging.error(err_msg)
        raise mp.ProcessError(err_msg)

    logging.info('......Finish rank transformation......')

    _write_metadata(file_path, offsets)

    session.config.rt_finished.value = True
    logging.debug('Mounting ranked spectra file')
    session.mount_ranked_spectra()
    return


def _reporter(dispatcher, index_length, log_lock, exit_signal):
    progress = 0
    while True:
        try:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Reporter thread exits now.')
                break
            time.sleep(0.2)
            old_progress = progress
            now = dispatcher.now()
            progress = int(now / index_length * 100)
            if progress > old_progress:
                with log_lock:
                    logging.info('Transforming progress: {}% ({} out of {})' \
                                 .format(progress, now, index_length))
        except KeyboardInterrupt:
            with log_lock:
                logging.debug('Reporter thread received KeyboardInterrupt, exits now.')
            break
    return


def _write_metadata(file_path, offsets):
    metadata = {'version': VERSION,
                'file_type': HEADER_LABEL_RANKED_SPECTRA,
                'magic_label': uuid4(),
                'session_label': session.magic_label,
                'index_label': session.internal_index.magic_label,
                'length': len(session.internal_index),
                'num_of_peaks': num_of_peaks,
                'mz_data_type': RT_MZ_DATA_TYPE,
                'intensity_data_type': RT_INTENSITY_DATA_TYPE,
                'mz_offset': offsets[0],
                'intensity_offset': offsets[1]}
    header = pickle.dumps(metadata) + b'%end%'
    if len(header) > FILE_HEADER_SIZE:
        err_msg = '\nHeader having size larger than the reserved area.'
        logging.error(err_msg)
        raise Exception(err_msg)
    with file_path.open('r+b') as fp:
        fp.write(header)
    return


def _worker(pid, dispatcher, offsets, log_lock, merge_lock, exit_signal):
    try:
        logging_setup()
        with log_lock:
            logging.debug('Rank transformation subprocess {} started.'.format(pid))

        file_path = Path.cwd().joinpath(session.name + FILE_EXTENSION_RANKED_SPECTRA).resolve()
        index_length = len(session.internal_index)
        mz_arr = np.memmap(str(file_path), offset=offsets[0], dtype=RT_MZ_DATA_TYPE, mode='r+',\
                           shape=(index_length, num_of_peaks))
        intensity_arr = np.memmap(str(file_path), offset=offsets[1], dtype=RT_INTENSITY_DATA_TYPE, mode='r+',\
                                  shape=(index_length, num_of_peaks))
        for chunk in dispatcher:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Subprocess {}: Received exit signal, exits now.'.format(pid))
                break

            try:
                start, end = chunk
                size = end - start
                with log_lock:
                    logging.debug('Subprocess {}: Transforming spectra {} - {} with size of {}.'\
                                  .format(pid, start, end, size))
                temp_index = session.internal_index[start:end]
                temp_mz = np.empty((size, num_of_peaks), dtype=np.int32)
                temp_intensity = np.empty((size, num_of_peaks), dtype=np.float32)

                for i in range(size):
                    file, offset = session.ms_exp_files[temp_index.file_id[i]], temp_index.offset[i]
                    precursor_mass, precursor_charge = temp_index.precursor_mass[i], temp_index.precursor_charge[i]
                    format_ = ''.join(file.suffixes).lower()
                    if format_ == '.mzxml': peaks = mzxml.get_peaks(file, offset, log_lock)
                    elif format_ == '.mzml': peaks = mzml.get_peaks(file, offset, log_lock)
                    temp_mz[i], temp_intensity[i] = _rank_transform(precursor_mass, precursor_charge,
                                                                    peaks[0], peaks[1], log_lock)

                with merge_lock:
                    mz_arr[start:end], intensity_arr[start:end] = temp_mz, temp_intensity
                    mz_arr.flush()
                    intensity_arr.flush()
            except Exception:
                err_msg = '\nSubprocess {}: Failed to transform spectrum located at {}'.format(pid, start+i)
                with log_lock:
                    logging.error(err_msg)
                raise
        with log_lock:
            logging.debug('Subprocess {}: Reached the end of iteration, work done.'.format(pid))
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


def _rank_transform(precursor_mass, precursor_charge, mz, intensity, log_lock=Lock()):
    if RT_INPUT_HIGH_PRECISION:
        mz = mz if mz.dtype == np.float64 else mz.astype(np.float64)
        intensity = intensity if intensity.dtype == np.float64 else intensity.astype(np.float64)
    else:
        mz = mz if mz.dtype == np.float32 else mz.astype(np.float32)
        intensity = intensity if intensity.dtype == np.float32 else intensity.astype(np.float32)
    # check is sorted
    if len(mz) > 1 and not np.all(mz[1:] > mz[:-1]):
        err_msg = '\nm/z values are not sorted. This may be a broken spectrum.'
        if not ignore_errors:
            with log_lock:
                logging.error(err_msg)
        raise ValueError(err_msg)

    if mz_range != (0, 0):
        idx = np.logical_and(mz >= mz_range[0], mz <= mz_range[1])
        mz, intensity = mz[idx], intensity[idx]

    if precursor_removal_range != (0.0, 0.0):
        precursor_mass = precursor_mass / precursor_charge if true_precursor_mass else precursor_mass
        idx = np.logical_or(mz < (precursor_removal_range[0] + precursor_mass),
                            mz > (precursor_removal_range[1] + precursor_mass))
        mz, intensity = mz[idx], intensity[idx]

    # mimic the effect of smaller bin size by multiplying the mz values.
    if bins_per_th != 1: mz *= bins_per_th
    new_mz = np.empty(len(mz), dtype=np.int32)
    new_intensity = np.empty(len(intensity), dtype=np.float32)
    _binning(mz, intensity, new_mz, new_intensity)
    mz, intensity = new_mz, new_intensity

    temp = max(0, num_of_peaks - len(mz))
    new_intensity = np.arange(num_of_peaks, temp, -1)
    idx = np.argsort(intensity)[:-(num_of_peaks+1):-1]
    intensity[idx] = new_intensity
    idx.sort()

    mz = mz[idx]
    intensity = intensity[idx]
    intensity /= np.linalg.norm(intensity)

    if temp != 0:
        mz = np.append(np.full(temp, -1, dtype=np.int32), mz)
        intensity = np.append(np.full(temp, 0, dtype=np.float32), intensity)
    return mz, intensity


# @nb.jit(nopython=True)
# def _binning(mz, intensity, new_mz, new_intensity):
#     previous = -1
#     cursor = -1
#     length = len(mz)
#     for i in range(length):
#         mz_int = int(mz[i])
#         if mz_int == previous:
#             new_intensity[cursor] += intensity[i]
#         else:
#             cursor += 1
#             new_mz[cursor] = previous = mz_int
#             new_intensity[cursor] = intensity[i]
#     cursor += 1
#     return


def _refresh_session():
    global ignore_errors, true_precursor_mass, num_of_peaks, remove_precursor
    global precursor_removal_range, mz_range, bins_per_th, num_of_threads
    ignore_errors = session.flags.ignore_errors
    true_precursor_mass = session.internal_index.true_precursor_mass
    num_of_peaks = session.config.rt_num_of_peaks.value
    precursor_removal_range = session.config.rt_precursor_removal_range.value
    mz_range = session.config.rt_mz_range.value
    bins_per_th = session.config.rt_bins_per_th.value
    num_of_threads = session.config.rt_num_of_threads.value
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
