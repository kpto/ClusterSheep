# -*- coding: utf-8 -*-
"""
Created on 01:56:11 01/12/2017

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
import pickle
import io
import os
import multiprocessing as mp
from threading import Thread
import time
import ctypes
import sys
import sqlite3
import numpy as np

from ClusterSheep.envr.session import get_session
from ClusterSheep.prcs.parallel.logging_setup import logging_setup
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def export_cluster(file, num_of_threads=os.cpu_count(), specific_cluster=None):
    _refresh_session()

    processes = []
    log_lock = mp.Lock()
    count_lock = mp.Lock()
    read_lock = mp.Lock()
    merge_lock = mp.Lock()
    exit_signal = mp.Value(ctypes.c_bool, False)
    finish_count = mp.Value(ctypes.c_uint64, 0)

    if not file.exists():
        with file.open('w', encoding='utf-8') as fp:
            fp.write(_get_header_line('Columns'))
            fp.write(_get_header_line('Cluster'))
            fp.write('ID\tNum of nodes\tNum of edges\tNum of identifications\t' +
                     'Major identification\tIdentified ratio\tAverage precursor mass\n')
            fp.write(_get_header_line('Nodes'))
            fp.write('ID\tFile\tScan num\tIdentification\tProbability\n')
            fp.write(_get_header_line('Edges'))
            fp.write('ID of source\tID of target\tDot product\n')
            fp.write('\n' + _get_header_line('Clusters'))

    if specific_cluster:
        _write_cluster(None, session.clusters.cursor, specific_cluster, file, log_lock, read_lock, merge_lock)
        return

    session.clusters.connect()
    num_clusters = session.clusters.cursor.execute('SELECT "num_of_clusters" FROM "metadata"').fetchone()[0]
    session.clusters.disconnect()

    for pid in range(num_of_threads):
        processes.append(mp.Process(target=_worker, args=(pid, file, num_clusters, finish_count, log_lock,
                                                          count_lock, read_lock, merge_lock, exit_signal)))
    reporter = Thread(target=_reporter, args=(finish_count, num_clusters, log_lock, exit_signal))

    logging.info('......Start cluster export with {} subprocesses......'.format(num_of_threads))

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
                logging.info('Received KeyboardInterrupt, cluster enrichment canceled.')
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

    logging.info('......Finish cluster export......')

    logging.debug('Reconnecting clusters file')
    session.clusters.connect()
    return


def _reporter(finish_count, num_clusters, log_lock, exit_signal):
    progress = 0
    while True:
        try:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Reporter thread exits now.')
                break
            time.sleep(0.2)
            old_progress = progress
            progress = int(finish_count.value / num_clusters * 100)
            if progress > old_progress:
                with log_lock:
                    logging.info('Export progress: {}% ({} out of {})' \
                                 .format(progress, finish_count.value, num_clusters))
        except KeyboardInterrupt:
            with log_lock:
                logging.debug('Reporter thread received KeyboardInterrupt, exits now.')
            break
    return


def _worker(pid, file, num_clusters, finish_count, log_lock, count_lock, read_lock, merge_lock, exit_signal):
    try:
        logging_setup()
        with log_lock:
            logging.debug('Cluster export subprocess {} started.'.format(pid))

        clusters_conn = sqlite3.connect(str(session.clusters.file_path))
        clusters_cur = clusters_conn.cursor()

        while True:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Subprocess {}: Received exit signal, exits now.'.format(pid))
                break

            with count_lock:
                cluster_id = finish_count.value + 1
                if cluster_id > num_clusters: break
                finish_count.value += 1
            _write_cluster(pid, clusters_cur, cluster_id, file, log_lock, read_lock, merge_lock)

        clusters_conn.close()
        logging.debug('Subprocess {}: Reached the end of iteration, work done.'.format(pid))
    except (Exception, KeyboardInterrupt) as e:
        if type(e) is KeyboardInterrupt:
            clusters_conn.close()
            with log_lock:
                logging.debug('Subprocess {}: Received KeyboardInterrupt, exits now.'.format(pid))
        else:
            clusters_conn.close()
            with log_lock:
                logging.exception('\nSubprocess {}: Ended unexpectedly. Logging traceback:\n'
                                  '==========TRACEBACK==========\n'.format(pid))
            exit_signal.value = True
            sys.exit(1)
    return


def _write_cluster(pid, clusters_cur, cluster_id, file, log_lock, read_lock, merge_lock):
    with read_lock:
        cluster = clusters_cur.execute('SELECT * FROM "clusters" WHERE "cluster_id"=?',
                                       (cluster_id,)).fetchone()
    if not cluster:
        header = 'Subprocess {}: '.format(pid) if pid is not None else ''
        wrn_msg = '{}Cluster with id "{}" does not exist.'.format(header, cluster_id)
        with log_lock:
            logging.warning(wrn_msg)
        return

    graph = pickle.loads(cluster[-1])
    internal_ids = graph.vp['iid'].a
    files = session.ms_exp_files[session.internal_index.file_id[internal_ids]]
    native_ids = session.internal_index.native_id[internal_ids]
    if 'ide' in graph.gp:
        idens = graph.gp['ide']
        idens[-1] = 'None'
        iden_ids = graph.vp['ide'].a
        prb = graph.vp['prb'].a
        iden_provider = lambda i: idens[iden_ids[i]]
        prb_provider = lambda i: prb[i]
    else:
        iden_provider = lambda i: None
        prb_provider = lambda i: None
    with merge_lock:
        with file.open('a', encoding='utf-8') as fp:
            fp.write(_get_header_line('Cluster'))
            fp.write(_join_elements(cluster[:-1]))
            fp.write(_get_header_line('Nodes'))
            lines = []
            for i in range(graph.num_vertices()):
                lines.append((internal_ids[i], files[i], native_ids[i], iden_provider(i), prb_provider(i)))
            lines = sorted(lines, key=lambda x: x[2])
            lines = sorted(lines, key=lambda x: x[1])
            fp.writelines(map(_join_elements, lines))
            fp.write(_get_header_line('Edges'))
            edges = graph.get_edges()
            roted = np.rot90(edges)
            sa = np.empty(len(edges), dtype=[('s', np.uint64), ('t', np.uint64), ('dp', np.float32)])
            sa['s'] = internal_ids[roted[1]]
            sa['t'] = internal_ids[roted[0]]
            sa['dp'] = graph.ep['dps'].a
            sio = io.StringIO()
            np.savetxt(sio, sa, fmt=['%u', '%u', '%.8f'], delimiter='\t')
            sio.seek(0)
            fp.write(sio.read())
            fp.write('\n')
    return


def _get_header_line(header):
    return '# {}\n'.format(header)


def _join_elements(line):
    return '\t'.join(map(str, line)) + '\n'


def _refresh_session():
    return


_refresh_session()
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
