# -*- coding: utf-8 -*-
"""
Created on 06:58:39 30/11/2017

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
import os
import multiprocessing as mp
from threading import Thread
import time
import ctypes
import sys
import sqlite3

from envr.session import get_session
from prcs.parallel.logging_setup import logging_setup
from property import *
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
def enrich_clusters(update=False, num_of_threads=os.cpu_count()):
    _refresh_session()
    if not session.iden_lut:
        logging.info('No identification lookup table is mounted, cannot proceed.')
        return
    processes = []
    log_lock = mp.Lock()
    count_lock = mp.Lock()
    read_lock = mp.Lock()
    merge_lock = mp.Lock()
    exit_signal = mp.Value(ctypes.c_bool, False)
    finish_count = mp.Value(ctypes.c_uint64, 1)

    session.clusters.connect()
    num_clusters = session.clusters.cursor.execute('SELECT "num_of_clusters" FROM "metadata"').fetchone()[0]
    session.clusters.disconnect()
    session.iden_lut.disconnect()

    for pid in range(num_of_threads):
        processes.append(mp.Process(target=_worker, args=(pid, update, num_clusters, finish_count, log_lock,
                                                          count_lock, read_lock, merge_lock, exit_signal)))
    reporter = Thread(target=_reporter, args=(finish_count, num_clusters, log_lock, exit_signal))

    logging.info('......Start cluster enrichment with {} subprocesses......'.format(num_of_threads))

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

    logging.info('......Finish cluster enrichment......')

    logging.debug('Reconnecting clusters  and identification lookup table file')
    session.clusters.connect()
    session.iden_lut.connect()
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
                    logging.info('Enrichment progress: {}% ({} out of {})' \
                                 .format(progress, finish_count.value, num_clusters))
        except KeyboardInterrupt:
            with log_lock:
                logging.debug('Reporter thread received KeyboardInterrupt, exits now.')
            break
    return


def _worker(pid, update, num_clusters, finish_count, log_lock, count_lock, read_lock, merge_lock, exit_signal):
    try:
        logging_setup()
        with log_lock:
            logging.debug('Cluster enrichment subprocess {} started.'.format(pid))

        clusters_conn = sqlite3.connect(str(session.clusters.file_path))
        clusters_cur = clusters_conn.cursor()
        session.iden_lut.connect()

        while True:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Subprocess {}: Received exit signal, exits now.'.format(pid))
                break

            with count_lock:
                cluster_id = finish_count.value
                if cluster_id > num_clusters: break
                finish_count.value += 1
            _add_iden_stat(pid, clusters_conn, clusters_cur, cluster_id, update, log_lock, read_lock, merge_lock)

        clusters_conn.close()
        session.iden_lut.disconnect()
        logging.debug('Subprocess {}: Reached the end of iteration, work done.'.format(pid))
    except (Exception, KeyboardInterrupt) as e:
        if type(e) is KeyboardInterrupt:
            clusters_conn.close()
            session.iden_lut.disconnect()
            with log_lock:
                logging.debug('Subprocess {}: Received KeyboardInterrupt, exits now.'.format(pid))
        else:
            clusters_conn.close()
            session.iden_lut.disconnect()
            with log_lock:
                logging.exception('\nSubprocess {}: Ended unexpectedly. Logging traceback:\n'
                                  '==========TRACEBACK==========\n'.format(pid))
            exit_signal.value = True
            sys.exit(1)
    return


def _add_iden_stat(pid, clusters_conn, clusters_cur, cluster_id, update, log_lock, read_lock, merge_lock):
    with read_lock:
        query = clusters_cur.execute('SELECT "num_idens", "pickled" FROM "clusters" WHERE "cluster_id"=?',
                                     (cluster_id,)).fetchone()
    if not query:
        wrn_msg = 'Subprocess {}: Cluster with id "{}" does not exist.'.format(pid, cluster_id)
        with log_lock:
            logging.warning(wrn_msg)
        return
    num_idens, pickled = query
    if update or not num_idens:
        graph = pickle.loads(pickled)
        num_vertices = graph.num_vertices()
        iid_arr = graph.vp['iid'].a
        ide = graph.new_vp('int', val=-1)
        prb = graph.new_vp('float', val=0.0)
        precursor_mass = np.empty(num_vertices, dtype=np.float32)
        ide_arr = ide.a
        prb_arr = prb.a
        idens = {}

        for i in range(num_vertices):
            entry = session.internal_index[int(iid_arr[i])]
            precursor_mass[i] = entry.precursor_mass
            with read_lock:
                identification = entry.get_identification()
            if identification and not identification.is_decoy:
                tpp_string = identification.to_tpp_string()
                if tpp_string not in idens:
                    idens[tpp_string] = len(idens)
                ide_arr[i] = idens[tpp_string]
                prb_arr[i] = identification.probability
        pre_mass_avg = float(np.mean(precursor_mass))

        if len(idens) == 0:
            num_idens = 0
            iden_ratio = 0.0
            with merge_lock:
                clusters_cur.execute('UPDATE "clusters" '
                                     'SET "num_idens"=?, "iden_ratio"=?, "pre_mass_avg"=? '
                                     'WHERE "cluster_id"=?', (num_idens, iden_ratio, pre_mass_avg, cluster_id))
                clusters_conn.commit()
        else:
            num_idens = len(idens)
            graph.vp['ide'] = ide
            graph.vp['prb'] = prb
            idens = {v: k for k, v in idens.items()}
            graph.gp['ide'] = graph.new_gp('object')
            graph.gp['ide'] = idens
            iden_id, counts = np.unique(ide.a, return_counts=True)
            idx = np.where(iden_id == -1)[0]
            if len(idx) == 0:
                num_uniden = 0
            else:
                idx = idx[0]
                num_uniden = counts[idx]
                iden_id = np.append(iden_id[:idx], iden_id[idx+1:])
                counts = np.append(counts[:idx], counts[idx+1:])
            iden_id = iden_id[np.argsort(counts)[::-1]]
            graph.gp['ord'] = graph.new_gp('object')
            graph.gp['ord'] = list(iden_id)
            major_iden = idens[iden_id[0]]
            iden_ratio = (num_vertices - num_uniden) / num_vertices
            pickled = pickle.dumps(graph)
            with merge_lock:
                clusters_cur.execute('UPDATE "clusters" '
                                     'SET "num_idens"=?, "major_iden"=?, "iden_ratio"=?, "pre_mass_avg"=?, "pickled"=? '
                                     'WHERE "cluster_id"=?',
                                     (num_idens, major_iden, iden_ratio, pre_mass_avg, pickled, cluster_id))
                clusters_conn.commit()
    return


def _refresh_session():
    return


_refresh_session()
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
