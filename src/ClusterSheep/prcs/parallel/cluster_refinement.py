# -*- coding: utf-8 -*-
"""
Created on 13:37:26 10/09/2017

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
import math

from ClusterSheep.envr.session import get_session
from ClusterSheep.prcs.parallel.logging_setup import logging_setup
from ClusterSheep.property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    outlier_threshold = None
    keep_raw = None
    set_num_of_threads = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def refine_cluster():
    _refresh_session()
    if len(outlier_threshold) == 0:
        logging.info('Outlier threshold list is empty, skip cluster refinement.')
        return

    processes = []
    log_lock = mp.Lock()
    count_lock = mp.Lock()
    read_lock = mp.Lock()
    merge_lock = mp.Lock()
    exit_signal = mp.Value(ctypes.c_bool, False)
    finish_count = mp.Value(ctypes.c_uint64, 0)
    total_count = mp.Value(ctypes.c_uint64, 0)

    raw_clusters = session.clusters.file_path
    raw_clusters = sqlite3.connect(str(raw_clusters.absolute()))
    cur = raw_clusters.cursor()
    num_raw_clusters = cur.execute('SELECT "num_of_clusters" FROM "metadata"').fetchone()[0]
    raw_clusters.close()

    refined_clusters = Path.cwd().joinpath(session.name + FILE_EXTENSION_TEMPORARY).absolute()
    if refined_clusters.exists():
        refined_clusters.unlink()
    refined_clusters = sqlite3.connect(str(refined_clusters))
    cur = refined_clusters.cursor()
    cur.execute('''
        CREATE TABLE "metadata" (
        "version" TEXT,
        "file_type" TEXT,
        "num_of_clusters" INTEGER,
        "magic_label" BLOB,
        "session_label" BLOB,
        "index_label" BLOB
        )''')
    metadata = (VERSION, HEADER_LABEL_CLUSTERS, 0, pickle.dumps(uuid4()),
                pickle.dumps(session.magic_label), pickle.dumps(session.internal_index.magic_label))
    cur.execute('INSERT INTO "metadata" VALUES (?, ?, ?, ?, ?, ?)', metadata)
    cur.execute('''
        CREATE TABLE "clusters" (
        "cluster_id" INTEGER PRIMARY KEY AUTOINCREMENT,
        "num_nodes" INTEGER,
        "num_edges" INTEGER,
        "num_idens" INTEGER,
        "major_iden" TEXT,
        "iden_ratio" REAL,
        "pre_mass_avg" REAL,
        "pickled" BLOB
        )''')
    refined_clusters.commit()
    refined_clusters.close()

    chunk_size = math.ceil(num_raw_clusters / set_num_of_threads)
    ranges = [(start, min(num_raw_clusters, start+chunk_size)) for start in range(0, num_raw_clusters, chunk_size)]
    num_of_threads = set_num_of_threads if len(ranges) == set_num_of_threads else len(ranges)
    for pid in range(num_of_threads):
        processes.append(mp.Process(target=_worker, args=(pid, ranges[pid], finish_count, total_count, log_lock,
                                                          count_lock, read_lock, merge_lock, exit_signal)))
    reporter = Thread(target=_reporter, args=(finish_count, num_raw_clusters, log_lock, exit_signal))

    logging.info('......Start cluster refinement with {} subprocesses......'.format(num_of_threads))

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
                logging.info('Received KeyboardInterrupt, identification import canceled.')
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

    if finish_count.value == 0:
        err_msg = '\nNo cluster has been made.'
        logging.error(err_msg)
        raise Exception(err_msg)

    if keep_raw:
        Path.cwd().joinpath(session.name + FILE_EXTENSION_CLUSTERS).rename(session.name + FILE_EXTENSION_RAW_CLUSTERS)
    Path.cwd().joinpath(session.name + FILE_EXTENSION_TEMPORARY).rename(session.name + FILE_EXTENSION_CLUSTERS)

    session.config.cr_finished.value = True
    logging.info('......Finish cluster refinement, formed {} clusters......'.format(total_count.value))
    return finish_count.value


def _reporter(finish_count, num_raw_clusters, log_lock, exit_signal):
    progress = 0
    while True:
        try:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Reporter thread exits now.')
                break
            time.sleep(0.2)
            old_progress = progress
            progress = int(finish_count.value / num_raw_clusters * 100)
            if progress > old_progress:
                with log_lock:
                    logging.info('Refinement progress: {}% ({} out of {})'
                                 .format(progress, finish_count.value, num_raw_clusters))
        except KeyboardInterrupt:
            with log_lock:
                logging.debug('Reporter thread received KeyboardInterrupt, exits now.')
            break
    return


def _worker(pid, range_, finish_count, total_count, log_lock, count_lock, read_lock, merge_lock, exit_signal):
    try:
        logging_setup()
        with log_lock:
            logging.debug('Cluster refinement subprocess {} started.'.format(pid))
            logging.debug('Received range: {}'.format(range_))

        raw_clusters = session.clusters.file_path
        raw_clusters = sqlite3.connect(str(raw_clusters))
        raw_cur = raw_clusters.cursor()
        with read_lock:
            raw_cur.execute('SELECT "pickled" FROM "clusters" ORDER BY "rowid" LIMIT ? OFFSET ?',
                            (range_[1] - range_[0], range_[0]))
        refined_clusters = Path.cwd().joinpath(session.name + FILE_EXTENSION_TEMPORARY).absolute()
        refined_clusters = sqlite3.connect(str(refined_clusters))
        refined_cur = refined_clusters.cursor()

        # push to database if graphs made storing in RAM exceed GM_BUFFER_SIZE
        bytes_count = 0
        data = []

        while True:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Subprocess {}: Received exit signal, exits now.'.format(pid))
                break

            with read_lock:
                clusters = raw_cur.fetchone()
            if clusters:
                clusters = [pickle.loads(clusters[0])]
            else:
                _push(pid, refined_clusters, refined_cur, data, log_lock, merge_lock)
                break

            for threshold in outlier_threshold:
                for i in range(len(clusters)):
                    clusters[i] = _refinement(clusters[i], threshold)
                clusters = [c for sublist in clusters for c in sublist]
            for c in clusters:
                pickled = pickle.dumps(c)
                bytes_count += len(pickled)
                data.append((c.num_vertices(), c.num_edges(), pickled))

            with count_lock:
                finish_count.value += 1
                total_count.value += len(clusters)

            if bytes_count >= GM_BUFFER_SIZE:
                _push(pid, refined_clusters, refined_cur, data, log_lock, merge_lock)
                data.clear()
                bytes_count = 0

        raw_clusters.close()
        refined_clusters.close()
        logging.debug('Subprocess {}: Reached the end of iteration, work done.'.format(pid))
    except (Exception, KeyboardInterrupt) as e:
        if type(e) is KeyboardInterrupt:
            raw_clusters.close()
            refined_clusters.close()
            with log_lock:
                logging.debug('Subprocess {}: Received KeyboardInterrupt, exits now.'.format(pid))
        else:
            raw_clusters.close()
            refined_clusters.close()
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


def _refinement(graph, threshold):
    vertex_betweenness_value = gt.betweenness(graph)[0].get_array()

    d = np.abs(vertex_betweenness_value - np.median(vertex_betweenness_value))
    mdev = np.median(d)
    s = d / mdev if mdev else np.zeros_like(d)
    vfilt = s < threshold

    graph =  gt.GraphView(graph, vfilt=vfilt)
    comp, hist = gt.label_components(graph)
    temp = []
    for i in range(len(hist)):
        if hist[i] > 1:
            temp.append(gt.Graph(gt.GraphView(graph, vfilt=(comp.a == i)), prune=True, directed=False))
    return temp


def _refresh_session():
    global outlier_threshold, keep_raw, set_num_of_threads
    outlier_threshold = session.config.cr_outlier_threshold.value
    keep_raw = session.config.cr_keep_raw.value
    set_num_of_threads = session.config.cr_num_of_threads.value
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
