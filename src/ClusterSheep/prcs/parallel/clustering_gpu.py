# -*- coding: utf-8 -*-
"""
Created on 13:37:16 10/09/2017

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
import pycuda.driver as drv
from pycuda.compiler import SourceModule
import math
import time
import ctypes

from ClusterSheep.envr.session import get_session
from ClusterSheep.property import *
from ClusterSheep.prcs.parallel.logging_setup import logging_setup
from ClusterSheep.prcs.parallel.cuda_source import get_source_code
from ClusterSheep.share.misc import _clean_temporary
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
    allocation_size_initial_divisor = None
    gpus_used = None
    cuda_block_dimensions = None
    use_fmad = None
    processes_per_device = None
    threads_per_device = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def clustering_gpu(dispatcher, temp_storage):
    _refresh_session()
    processes = []
    log_lock = mp.Lock()
    merge_lock = mp.Lock()
    exit_signal = mp.Value(ctypes.c_bool, False)
    total_edge_count = mp.Value(ctypes.c_uint64, 0)

    gpus = _get_gpus()
    if len(gpus) == 0: raise Exception('\nNo suitable CUDA device available.')

    num_of_threads = processes_per_device * len(gpus)
    for pid in range(num_of_threads):
        processes.append(mp.Process(target=_worker, args=(pid, gpus[pid%len(gpus)], dispatcher, temp_storage, total_edge_count,
                                                          log_lock, merge_lock, exit_signal)))
    reporter = Thread(target=_reporter, args=(dispatcher, math.ceil(len(session.internal_index)/block_dimensions[0]),
                                              log_lock, exit_signal))

    logging.info('......Start clustering using GPU with {} subprocesses......'.format(num_of_threads))

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
                logging.info('Received KeyboardInterrupt, clustering canceled.')
        exit_signal.value = True
        with log_lock:
            logging.debug('Waiting processes to exit.')
        if reporter.is_alive():
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


def _worker(pid, did, dispatcher, temp_storage, total_edge_count, log_lock, merge_lock, exit_signal):
    try:
        logging_setup()
        with log_lock:
            logging.debug('Clustering subprocess {} acquiring device {} started.'.format(pid, did))

        drv.init()
        cuda_device = drv.Device(did)
        cuda_context = cuda_device.make_context()
        compiler_option = ['--fmad=true',] if use_fmad else ['--fmad=false',]
        cuda_module = SourceModule(get_source_code(), options=compiler_option)
        cuda_kernel = cuda_module.get_function('compute_dot_product')
        cuda_kernel.prepare('PPPPPPPPPP')

        threads = []
        exit_state = mp.Value(ctypes.c_uint32, 0)
        with log_lock:
            logging.debug('Subprocess {}: Spawning {} threads for CUDA stream concurrency.'.format(pid, threads_per_device))
        for tid in range(threads_per_device):
            threads.append(Thread(target=_thread,
                                  args=(pid, tid, cuda_context, cuda_kernel, dispatcher, temp_storage,
                                        total_edge_count, log_lock, merge_lock, exit_signal, exit_state)))
        for t in threads: t.start()
        for t in threads: t.join()

        if exit_state.value == 1:
            err_msg = '\nSubprocess {}: Threads exited with abnormal exitcode.'.format(pid)
            logging.error(err_msg)
            raise Exception(err_msg)

        drv.stop_profiler()
        cuda_context.pop()
    except (Exception, KeyboardInterrupt) as e:
        if type(e) is KeyboardInterrupt:
            with log_lock:
                logging.debug('Subprocess {}: Received KeyboardInterrupt, exits now.'.format(pid))
                logging.debug('Subprocess {}: Waiting threads to exit.'.format(pid))
            for t in threads: t.join()
            drv.stop_profiler()
            cuda_context.pop()
        else:
            with log_lock:
                logging.exception('\nSubprocess {}: Ended unexpectedly. Logging traceback:\n'
                                  '==========TRACEBACK==========\n'.format(pid))
            drv.stop_profiler()
            cuda_context.pop()
            exit_signal.value = True
            sys.exit(1)
    return


def _thread(pid, tid, cuda_context, cuda_kernel, dispatcher, temp_storage, total_edge_count, log_lock, merge_lock, exit_signal, exit_state):
    try:
        with log_lock:
            logging.debug('Clustering subprocess {} thread {} started.'.format(pid, tid))

        cuda_context.push()

        ref_block_height, ref_block_width = block_dimensions
        edg_path = Path(temp_storage, 'edg')
        dps_path = Path(temp_storage, 'dps')
        ranked_spectra = session.ranked_spectra
        cuda_stream = drv.Stream()

        allocation_size_divisor = allocation_size_initial_divisor
        allocation_size = int(ref_block_height * ref_block_width / allocation_size_divisor)
        reallocated = False

        with log_lock:
            logging.debug('Clustering subprocess {} thread {}: Allocating host and device memory.'.format(pid, tid))
        # allocate host pagelocked memory
        # input
        plm_precursor_mass = drv.pagelocked_empty(ref_block_height + ref_block_width, dtype=CG_PRECURSOR_MASS_DATA_TYPE)
        plm_mz = drv.pagelocked_empty((ref_block_height + ref_block_width, num_of_peaks), dtype=CG_MZ_DATA_TYPE)
        plm_intensity = drv.pagelocked_empty((ref_block_height + ref_block_width, num_of_peaks), dtype=CG_INTENSITY_DATA_TYPE)
        plm_block_dimensions = drv.pagelocked_empty(2, dtype=CG_BLOCK_DIMENSIONS_DATA_TYPE)
        plm_offset = drv.pagelocked_empty(2, dtype=CG_OFFSET_DATA_TYPE)
        plm_allocation_size = drv.pagelocked_empty(1, dtype=CG_ALLOCATION_SIZE_DATA_TYPE)
        # output
        plm_counter = drv.pagelocked_empty(1, dtype=CG_COUNTER_DATA_TYPE)
        plm_edge = drv.pagelocked_empty((allocation_size, 2), dtype=CG_EDGE_DATA_TYPE)
        plm_dot_product = drv.pagelocked_empty(allocation_size, dtype=CG_DOT_PRODUCT_DATA_TYPE)
        plm_overflowed = drv.pagelocked_empty(1, dtype=CG_OVERFLOWED_DATA_TYPE)

        # allocate device memory
        # input
        dvp_precursor_mass = drv.mem_alloc_like(plm_precursor_mass)
        dvp_mz = drv.mem_alloc_like(plm_mz)
        dvp_intensity = drv.mem_alloc_like(plm_intensity)
        dvp_block_dimensions = drv.mem_alloc_like(plm_block_dimensions)
        dvp_offset = drv.mem_alloc_like(plm_offset)
        dvp_allocation_size = drv.mem_alloc_like(plm_allocation_size)
        # output
        dvp_counter = drv.mem_alloc_like(plm_counter)
        dvp_edge = drv.mem_alloc_like(plm_edge)
        dvp_dot_product = drv.mem_alloc_like(plm_dot_product)
        dvp_overflowed = drv.mem_alloc_like(plm_overflowed)

        with log_lock:
            logging.debug('Clustering subprocess {} thread {}: Start iterating dispatcher.'.format(pid, tid))
        previous_row_id = -1
        dispatcher.connect(pid, tid)

        # iterate dispatcher to get blocks
        for row_id, column_id, block in dispatcher.iterate(pid, tid):
            if exit_signal.value:
                with log_lock:
                    logging.debug('Subprocess {} thread {}: Received exit signal, exits now.'.format(pid, tid))
                break

            try:
                y_range, x_range = block
                block_height = y_range[1] - y_range[0]
                block_width = x_range[1] - x_range[0]
                if row_id != previous_row_id:
                    with log_lock:
                        logging.debug('\033[92mSubprocess {} thread {}: Processing row {} (y:{}->{}).\033[0m'.
                                      format(pid, tid, row_id, *y_range))
                    previous_row_id = row_id

                # get necessary data
                plm_precursor_mass[:block_height] = ranked_spectra.precursor_mass[y_range[0]:y_range[1]]
                plm_precursor_mass[block_height:block_height+block_width] = ranked_spectra.precursor_mass[x_range[0]:x_range[1]]
                plm_mz[:block_height] = ranked_spectra.mz[y_range[0]:y_range[1]]
                plm_mz[block_height:block_height + block_width] = ranked_spectra.mz[x_range[0]:x_range[1]]
                plm_intensity[:block_height] = ranked_spectra.intensity[y_range[0]:y_range[1]]
                plm_intensity[block_height:block_height + block_width] = ranked_spectra.intensity[x_range[0]:x_range[1]]
                plm_block_dimensions[:] = (block_height, block_width)
                plm_offset[:] = (y_range[0], x_range[0])
                # upload data
                drv.memcpy_htod_async(dvp_precursor_mass, plm_precursor_mass, cuda_stream)
                drv.memcpy_htod_async(dvp_mz, plm_mz, cuda_stream)
                drv.memcpy_htod_async(dvp_intensity, plm_intensity, cuda_stream)
                drv.memcpy_htod_async(dvp_block_dimensions, plm_block_dimensions, cuda_stream)
                drv.memcpy_htod_async(dvp_offset, plm_offset, cuda_stream)

                if reallocated:
                    allocation_size_divisor = allocation_size_initial_divisor
                    allocation_size = int(ref_block_height * ref_block_width / allocation_size_divisor)
                    # reallocate host pagelocked memory
                    del plm_edge
                    del plm_dot_product
                    plm_edge = drv.pagelocked_empty((allocation_size, 2), dtype=CG_EDGE_DATA_TYPE)
                    plm_dot_product = drv.pagelocked_empty(allocation_size, dtype=CG_DOT_PRODUCT_DATA_TYPE)
                    # reallocate device memory
                    del dvp_edge
                    del dvp_dot_product
                    dvp_edge = drv.mem_alloc_like(plm_edge)
                    dvp_dot_product = drv.mem_alloc_like(plm_dot_product)
                    with log_lock:
                        logging.debug('\033[92mSubprocess {} thread {}: Reset memory allocation size divisor to {}.\033[0m'.
                                      format(pid, tid, allocation_size_divisor))
                    reallocated = False

                cublockdim = (cuda_block_dimensions[1], cuda_block_dimensions[0], 1)
                cugriddim = (math.ceil(block_width / cuda_block_dimensions[1]),
                             math.ceil(block_height / cuda_block_dimensions[0]))
                while True:
                    plm_allocation_size[0] = allocation_size
                    plm_counter[0] = 0
                    plm_overflowed[0] = False
                    drv.memcpy_htod_async(dvp_allocation_size, plm_allocation_size, cuda_stream)
                    drv.memcpy_htod_async(dvp_counter, plm_counter, cuda_stream)
                    drv.memcpy_htod_async(dvp_overflowed, plm_overflowed, cuda_stream)

                    cuda_kernel.prepared_async_call(cugriddim, cublockdim, cuda_stream,
                                                    dvp_precursor_mass,
                                                    dvp_mz,
                                                    dvp_intensity,
                                                    dvp_block_dimensions,
                                                    dvp_offset,
                                                    dvp_allocation_size,
                                                    dvp_counter,
                                                    dvp_edge,
                                                    dvp_dot_product,
                                                    dvp_overflowed)

                    # transfer computation result from device to host
                    drv.memcpy_dtoh_async(plm_edge, dvp_edge, cuda_stream)
                    drv.memcpy_dtoh_async(plm_counter, dvp_counter, cuda_stream)
                    drv.memcpy_dtoh_async(plm_overflowed, dvp_overflowed, cuda_stream)
                    drv.memcpy_dtoh_async(plm_dot_product, dvp_dot_product, cuda_stream)
                    cuda_stream.synchronize()

                    if plm_overflowed[0]:
                        allocation_size_divisor = int(allocation_size_divisor / 2)
                        if allocation_size_divisor < 1:
                            err_msg = ('\nSubprocess {} thread {}: Allocation size divisor reached to the impossible value of {}.'
                                       .format(pid, tid, allocation_size_divisor))
                            with log_lock:
                                logging.error(err_msg)
                            raise Exception(err_msg)
                        with log_lock:
                            logging.debug('\033[92mSubprocess {} thread {}: Edge list overflowed, '
                                          'decreases allocation size divisor to {}.\033[0m'
                                          .format(pid, tid, allocation_size_divisor))
                        allocation_size = int(block_width * block_height / allocation_size_divisor)
                        # reallocate host pagelocked memory
                        del plm_edge
                        del plm_dot_product
                        plm_edge = drv.pagelocked_empty((allocation_size, 2), dtype=CG_EDGE_DATA_TYPE)
                        plm_dot_product = drv.pagelocked_empty(allocation_size, dtype=CG_DOT_PRODUCT_DATA_TYPE)
                        # reallocate device memory
                        del dvp_edge
                        del dvp_dot_product
                        dvp_edge = drv.mem_alloc_like(plm_edge)
                        dvp_dot_product = drv.mem_alloc_like(plm_dot_product)
                        reallocated = True
                        continue
                    else:
                        break

                if abs(plm_precursor_mass[block_height-1] - plm_precursor_mass[block_height+block_width-1]) > precursor_tolerance:
                    dispatcher.next_row(pid, tid)
                with merge_lock:
                    edge_list_size = int(plm_counter[0])
                    if edge_list_size != 0:
                        total_edge_count.value += edge_list_size
                        edg = np.memmap(str(edg_path), dtype=CG_EDGE_DATA_TYPE, mode='r+', shape=(total_edge_count.value, 2))
                        dps = np.memmap(str(dps_path), dtype=CG_DOT_PRODUCT_DATA_TYPE, mode='r+', shape=total_edge_count.value)
                        edg[-edge_list_size:] = plm_edge[:edge_list_size]
                        dps[-edge_list_size:] = plm_dot_product[:edge_list_size]

            except Exception:
                err_msg = '\nSubprocess {} thread {}: Failed to clustering block (y:{}->{}, x:{}->{}).' \
                    .format(pid, tid, y_range[0], y_range[1], x_range[0], x_range[1])
                with log_lock:
                    logging.error(err_msg)
                raise

        with log_lock:
            if not exit_signal.value:
                logging.debug('Subprocess {} thread {}: Reached the end of iteration, work done.'.format(pid, tid))
        cuda_context.pop()

    except (Exception, KeyboardInterrupt) as e:
        if type(e) is KeyboardInterrupt:
            with log_lock:
                logging.debug('Subprocess {} thread {}: Received KeyboardInterrupt, exits now.'.format(pid, tid))
        else:
            with log_lock:
                logging.exception('\nSubprocess {} thread {}: Ended unexpectedly. Logging traceback:\n'
                                  '==========TRACEBACK==========\n'.format(pid, tid))
            exit_signal.value = True
            exit_state.value = 1
        cuda_context.pop()
    return


def _get_gpus():

    def get_information(q):
        drv.init()
        count = drv.Device.count()
        gpus = list(range(count)) if len(gpus_used) == 0 else sorted(list(set([d for d in gpus_used if d < count])))
        info = []
        for id_ in gpus:
            d = drv.Device(id_)
            info.append((id_, d.name(), d.total_memory(), d.compute_capability()))
        q.put(info)
        return

    q = mp.Queue()
    p = mp.Process(target=get_information, args=(q,))
    p.start()
    p.join()
    info = q.get()
    gpus = []
    for d in info:
        id_, name, mem, cc = d
        if cc[0] < 2:
            wrn_msg = 'Device "{}" does not meet the compute capability requirement (>= 2.0) and has been ignored.'\
                      .format(name)
            logging.warning(wrn_msg)
            continue
        logging.info('Selected device: {}    Memory: {} MiB    Compute capability: {}'\
                     .format(name, int(mem/(1024**2)), '.'.join(map(str, cc))))
        gpus.append(id_)
    return gpus


def _refresh_session():
    global ignore_errors, num_of_peaks, precursor_tolerance, dot_product_threshold
    global block_dimensions, allocation_size_initial_divisor, gpus_used, cuda_block_dimensions
    global use_fmad, processes_per_device, threads_per_device
    ignore_errors = session.flags.ignore_errors
    num_of_peaks = session.ranked_spectra.num_of_peaks
    precursor_tolerance = session.config.cg_precursor_tolerance.value
    dot_product_threshold = session.config.cg_dot_product_threshold.value
    block_dimensions = session.config.cg_block_dimensions.value
    allocation_size_initial_divisor = session.config.cg_allocation_size_initial_divisor.value
    gpus_used = session.config.gpu_gpus_used.value
    cuda_block_dimensions = session.config.gpu_cuda_block_dimensions.value
    use_fmad = session.config.gpu_use_fmad.value
    processes_per_device = session.config.gpu_processes_per_device.value
    threads_per_device = session.config.gpu_threads_per_device.value
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
