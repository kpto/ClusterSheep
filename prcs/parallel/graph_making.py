# -*- coding: utf-8 -*-
"""
Created on 12:16:58 23/11/2017

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
import sqlite3
import sys
import ctypes
import pickle
import graph_tool.all as gt
from uuid import uuid4

from envr.session import get_session
from prcs.parallel.logging_setup import logging_setup
from prcs.parallel.find_cluster import _scan
from property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    num_of_threads = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def make_graphs(temp_storage, edge_list_length):
    _refresh_session()
    processes = []
    log_lock = mp.Lock()
    scan_lock = mp.Lock()
    count_lock = mp.Lock()
    merge_lock = mp.Lock()
    exit_signal = mp.Value(ctypes.c_bool, False)
    scan_starting_idx = mp.Value(ctypes.c_uint64, 0)
    finish_count = mp.Value(ctypes.c_uint64, 0)

    clusters = Path.cwd().joinpath(session.name + FILE_EXTENSION_CLUSTERS).absolute()
    if clusters.exists():
        clusters.unlink()
    clusters = sqlite3.connect(str(clusters))
    cur = clusters.cursor()
    cur.execute('''
        CREATE TABLE "metadata" (
        "version" TEXT,
        "file_type" TEXT,
        "num_of_clusters" INTEGER,
        "true_precursor_mass" BOOL,
        "magic_label" BLOB,
        "session_label" BLOB,
        "index_label" BLOB
        )''')
    metadata = (VERSION, HEADER_LABEL_CLUSTERS, 0, session.internal_index.true_precursor_mass, pickle.dumps(uuid4()),
                pickle.dumps(session.magic_label), pickle.dumps(session.internal_index.magic_label))
    cur.execute('INSERT INTO "metadata" VALUES (?, ?, ?, ?, ?, ?, ?)', metadata)
    cur.execute('''
        CREATE TABLE "clusters" (
        "cluster_id" INTEGER PRIMARY KEY,
        "num_nodes" INTEGER,
        "num_edges" INTEGER,
        "num_idens" INTEGER,
        "major_iden" TEXT,
        "iden_ratio" REAL,
        "pre_mass_avg" REAL,
        "pickled" BLOB
        )''')
    clusters.commit()
    clusters.close()

    for pid in range(num_of_threads):
        processes.append(mp.Process(target=_worker, args=(pid, scan_starting_idx, temp_storage, edge_list_length,
                                                          finish_count, log_lock, scan_lock, count_lock, merge_lock, exit_signal)))
    reporter = Thread(target=_reporter, args=(finish_count, log_lock, exit_signal))

    logging.info('......Start graph making with {} subprocesses......'.format(num_of_threads))

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
                logging.info('Received KeyboardInterrupt, identification import canceled.')
        exit_signal.value = True
        with log_lock:
            logging.debug('Waiting processes to exit.')
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

    if finish_count.value == 0:
        err_msg = '\nNo cluster has been made.'
        logging.error(err_msg)
        raise Exception(err_msg)

    logging.info('......Finish graph making......')
    return finish_count.value


def _reporter(finish_count, log_lock, exit_signal):
    last = 1000
    while True:
        try:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Reporter thread exits now.')
                break
            time.sleep(0.2)
            current = finish_count.value
            if current > (last * 2):
                with log_lock:
                    logging.info('Make progress: {} graphs'.format(current))
                last = current
        except KeyboardInterrupt:
            with log_lock:
                logging.debug('Reporter thread received KeyboardInterrupt, exits now.')
            break
    return


def _worker(pid, scan_starting_idx, temp_storage, edge_list_length, finish_count, log_lock, scan_lock, count_lock, merge_lock, exit_signal):
    try:
        logging_setup()
        with log_lock:
            logging.debug('Graph making subprocess {} started.'.format(pid))

        clusters_path = Path.cwd().joinpath(session.name + FILE_EXTENSION_CLUSTERS).absolute()
        conn = sqlite3.connect(str(clusters_path))
        cur = conn.cursor()
        edg_path = temp_storage.joinpath('edg')
        dps_path = temp_storage.joinpath('dps')
        nbg_path = temp_storage.joinpath('nbg')
        edge = np.memmap(str(edg_path), dtype=CG_EDGE_DATA_TYPE, mode='c', shape=(edge_list_length, 2))
        dps = np.memmap(str(dps_path), dtype=CG_DOT_PRODUCT_DATA_TYPE, mode='c', shape=edge_list_length)
        belonging = np.memmap(str(nbg_path), dtype=FG_PARENT_DATA_TYPE, mode='c', shape=edge_list_length)

        # push to database if graphs made storing in RAM exceed GM_BUFFER_SIZE
        bytes_count = 0
        local_finish_count = 0
        data = []

        while True:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Subprocess {}: Received exit signal, exits now.'.format(pid))
                break

            with scan_lock:
                start = scan_starting_idx.value
                end = _scan(belonging, start)
                scan_starting_idx.value = end
            if start == end:
                _push(pid, conn, cur, data, log_lock, merge_lock)
                with count_lock:
                    finish_count.value += local_finish_count
                break

            temp_graph = gt.Graph(directed=False)
            temp_graph.vp['iid'] = temp_graph.add_edge_list(edge[start:end], hashed=True)
            temp_graph.ep['dps'] = dps_prop = temp_graph.new_ep('float')
            dps_prop.a = dps[start:end]
            pickled = pickle.dumps(temp_graph)
            bytes_count += len(pickled)
            data.append((temp_graph.num_vertices(), temp_graph.num_edges(), pickled))
            local_finish_count += 1

            if local_finish_count >= 1000:
                with count_lock:
                    finish_count.value += local_finish_count
                local_finish_count = 0

            if bytes_count >= GM_BUFFER_SIZE:
                _push(pid, conn, cur, data, log_lock, merge_lock)
                data.clear()
                bytes_count = 0

        conn.close()
        logging.debug('Subprocess {}: Reached the end of iteration, work done.'.format(pid))
    except (Exception, KeyboardInterrupt) as e:
        if type(e) is KeyboardInterrupt:
            conn.close()
            with log_lock:
                logging.debug('Subprocess {}: Received KeyboardInterrupt, exits now.'.format(pid))
        else:
            conn.close()
            with log_lock:
                logging.exception('\nSubprocess {}: Ended unexpectedly. Logging traceback:\n'
                                  '==========TRACEBACK==========\n'.format(pid))
            exit_signal.value = True
            sys.exit(1)
    return


def _push(pid, conn, cur, data, log_lock, merge_lock):
    if len(data) == 0: return
    with log_lock:
        logging.debug('Subprocess {}: Pushing {} clusters'.format(pid, len(data)))
    with merge_lock:
        cur.execute('BEGIN TRANSACTION')
        cur.executemany('INSERT INTO "clusters" (num_nodes, num_edges, pickled) VALUES (?, ?, ?)', data)
        num_of_clusters = cur.execute('SELECT "num_of_clusters" FROM "metadata"').fetchone()[0] + len(data)
        cur.execute('UPDATE "metadata" SET "num_of_clusters"=?', (num_of_clusters,))
        conn.commit()
    return


def _refresh_session():
    global num_of_threads
    num_of_threads = session.config.gm_num_of_threads.value
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
