# -*- coding: utf-8 -*-
"""
Created on 13:37:09 10/09/2017

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
from pathlib import Path
import sys
import math
import time
import ctypes

from envr.session import get_session
from property import *
from prcs.parallel.logging_setup import logging_setup
from prcs.parallel.cpu_kernel import _cpu_kernel
from share.misc import _clean_temporary
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    ignore_errors = None
    num_of_peaks = None
    precursor_tolerance = None
    dot_product_threshold = None
    block_dimensions = None
    num_of_threads = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def clustering_cpu(dispatcher, temp_storage):
    _refresh_session()
    processes = []
    log_lock = mp.Lock()
    merge_lock = mp.Lock()
    exit_signal = mp.Value(ctypes.c_bool, False)
    total_edge_count = mp.Value(ctypes.c_uint64, 0)

    for pid in range(num_of_threads):
        processes.append(mp.Process(target=_worker, args=(pid, dispatcher, temp_storage, total_edge_count,
                                                          log_lock, merge_lock, exit_signal)))
    reporter = Thread(target=_reporter, args=(dispatcher, math.ceil(len(session.internal_index)/block_dimensions[0]),
                                              log_lock, exit_signal))

    logging.info('......Start clustering using CPU with {} subprocesses......'.format(num_of_threads))

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
                logging.info('Received KeyboardInterrupt, clustering canceled.')
        exit_signal.value = True
        with log_lock:
            logging.debug('Waiting processes to exit.')
        reporter.join()
        for i in range(num_of_threads):
            if processes[i].is_alive():
                processes[i].join()
        logging.debug('All processes exited.')
        _clean_temporary(temp_storage)
        raise

    # in case that subprocesses were exited with errors
    if 1 in [p.exitcode for p in processes]:
        err_msg = '\nSubprocesses exited with abnormal exitcode.'
        logging.error(err_msg)
        _clean_temporary(temp_storage)
        raise mp.ProcessError(err_msg)
    return total_edge_count


def _reporter(dispatcher, num_rows, log_lock, exit_signal):
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
            progress = min(100, int(now / num_rows * 100))
            if progress > old_progress:
                with log_lock:
                    logging.info('Clustering progress: {}% ({} out of {})' \
                                 .format(progress, min(num_rows, now), num_rows))
        except KeyboardInterrupt:
            with log_lock:
                logging.debug('Reporter thread received KeyboardInterrupt, exits now.')
            break
    return


def _worker(pid, dispatcher, temp_storage, total_edge_count, log_lock, merge_lock, exit_signal):
    try:
        logging_setup()
        with log_lock:
            logging.debug('Clustering subprocess {} started.'.format(pid))

        edg_path = Path(temp_storage, 'edg')
        dps_path = Path(temp_storage, 'dps')
        ranked_spectra = session.ranked_spectra

        with log_lock:
            logging.debug('Clustering subprocess {}: Start iterating dispatcher.'.format(pid))
        previous_row_id = -1
        dispatcher.connect(pid, 0)

        # iterate dispatcher to get blocks
        for row_id, column_id, block in dispatcher.iterate(pid, 0):
            if exit_signal.value:
                with log_lock:
                    logging.debug('Subprocess {}: Received exit signal, exits now.'.format(pid))
                break

            try:
                y_range, x_range = block
                block_height = y_range[1] - y_range[0]
                block_width = x_range[1] - x_range[0]
                if row_id != previous_row_id:
                    with log_lock:
                        logging.debug('Subprocess {}: Processing row {} (y:{}->{}).'.
                                      format(pid, row_id, *y_range))
                    previous_row_id = row_id

                # get necessary data
                precursor_mass = np.concatenate((ranked_spectra.precursor_mass[y_range[0]:y_range[1]],
                                                 ranked_spectra.precursor_mass[x_range[0]:x_range[1]]))
                mz = np.concatenate((ranked_spectra.mz[y_range[0]:y_range[1]],
                                     ranked_spectra.mz[x_range[0]:x_range[1]]))
                intensity = np.concatenate((ranked_spectra.intensity[y_range[0]:y_range[1]],
                                            ranked_spectra.intensity[x_range[0]:x_range[1]]))
                block_dimensions = np.array((block_height, block_width), dtype=CG_BLOCK_DIMENSIONS_DATA_TYPE)
                offset = np.array((y_range[0], x_range[0]), dtype=CG_OFFSET_DATA_TYPE)
                edge = np.empty((block_height*block_width, 2), dtype=CG_EDGE_DATA_TYPE)
                dot_product = np.empty(block_height*block_width, dtype=CG_DOT_PRODUCT_DATA_TYPE)

                edge_list_size = _cpu_kernel(precursor_mass,
                                             mz,
                                             intensity,
                                             block_dimensions,
                                             offset,
                                             precursor_tolerance,
                                             dot_product_threshold,
                                             num_of_peaks,
                                             edge,
                                             dot_product)
                edge = edge[:edge_list_size]
                dot_product = dot_product[:edge_list_size]

                if abs(precursor_mass[block_height - 1] - precursor_mass[
                                    block_height + block_width - 1]) > precursor_tolerance:
                    dispatcher.next_row(pid, 0)
                with merge_lock:
                    if edge_list_size != 0:
                        total_edge_count.value += edge_list_size
                        edg = np.memmap(str(edg_path), dtype=CG_EDGE_DATA_TYPE, mode='r+', shape=(total_edge_count.value, 2))
                        dps = np.memmap(str(dps_path), dtype=CG_DOT_PRODUCT_DATA_TYPE, mode='r+', shape=total_edge_count.value)
                        edg[-edge_list_size:] = edge[:edge_list_size]
                        dps[-edge_list_size:] = dot_product[:edge_list_size]

            except Exception:
                err_msg = '\nSubprocess {}: Failed to clustering block (y:{}->{}, x:{}->{}).' \
                    .format(pid, y_range[0], y_range[1], x_range[0], x_range[1])
                with log_lock:
                    logging.error(err_msg)
                raise

        with log_lock:
            if not exit_signal.value:
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


# @nb.jit(nopython=True)
# def _cpu_kernel(precursor_mass, mz, intensity, block_dimensions, offset, dp, edge, dot_product):
#     count = 0
#     for location_y in range(block_dimensions[0]):
#         for location_x in range(block_dimensions[1]):
#             if abs(precursor_mass[location_y] - precursor_mass[block_dimensions[0] + location_x] <= precursor_tolerance):
#                 temp_mz_y = mz[location_y]
#                 temp_mz_x = mz[block_dimensions[0] + location_x]
#                 temp_intensity_y = intensity[location_y]
#                 temp_intensity_x = intensity[block_dimensions[0] + location_x]
#
#                 dp[0] = 0.0
#                 y_ptr = 0
#                 x_ptr = 0
#
#                 while y_ptr < num_of_peaks and x_ptr < num_of_peaks:
#                     if temp_mz_y[y_ptr] == temp_mz_x[x_ptr]:
#                         dp += temp_intensity_y[y_ptr] * temp_intensity_x[x_ptr]
#                         y_ptr += 1
#                         x_ptr += 1
#                     elif temp_mz_y[y_ptr] < temp_mz_x[x_ptr]:
#                         y_ptr += 1
#                     else:
#                         x_ptr += 1
#
#                 if dp[0] > dot_product_threshold:
#                     global_location_y = location_y + offset[0]
#                     global_location_x = location_x + offset[1]
#                     if global_location_x > global_location_y:
#                         dot_product[count] = dp[0]
#                         edge[count][0] = global_location_y
#                         edge[count][1] = global_location_x
#                         count += 1
#     return edge[:count], dot_product[:count], count


def _refresh_session():
    global ignore_errors, num_of_peaks, precursor_tolerance, dot_product_threshold
    global block_dimensions, num_of_threads
    ignore_errors = session.flags.ignore_errors
    num_of_peaks = session.ranked_spectra.num_of_peaks
    precursor_tolerance = session.config.cg_precursor_tolerance.value
    dot_product_threshold = session.config.cg_dot_product_threshold.value
    block_dimensions = session.config.cg_block_dimensions.value
    num_of_threads = session.config.cpu_num_of_threads.value
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
