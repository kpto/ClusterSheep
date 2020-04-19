# -*- coding: utf-8 -*-
"""
Created on 13:53:23 10/09/2017

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
from uuid import uuid4

from ClusterSheep.envr.session import get_session
from ClusterSheep.prcs.parallel.logging_setup import logging_setup
import ClusterSheep.reader.pepxml as pepxml
from ClusterSheep.property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    ignore_errors = None
    num_of_threads = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def import_identification():
    _refresh_session()
    if len(session.iden_files) == 0:
        logging.info('No identification file specified, skip identification import.')
        return

    pepxml._refresh_session()
    processes = []
    q = mp.Queue()
    log_lock = mp.Lock()
    merge_lock = mp.Lock()
    exit_signal = mp.Value(ctypes.c_bool, False)
    valid_iden_count = mp.Value(ctypes.c_uint64, 0)
    total_iden_count = mp.Value(ctypes.c_uint64, 0)
    finish_count = mp.Value(ctypes.c_uint64, 0)
    num_files = len(session.iden_files)

    iden_lut = Path.cwd().joinpath(session.name + FILE_EXTENSION_IDEN_LUT)
    if iden_lut.exists():
        iden_lut.unlink()
    iden_lut = sqlite3.connect(str(iden_lut.absolute()))
    cur = iden_lut.cursor()
    cur.execute('''
        CREATE TABLE "metadata" (
        "version" TEXT,
        "file_type" TEXT,
        "mods_pos_data_type" TEXT,
        "mods_mass_data_type" TEXT,
        "magic_label" BLOB,
        "session_label" BLOB
        )''')
    metadata = (VERSION, HEADER_LABEL_IDEN_LUT, ID_MODS_POS_DATA_TYPE, ID_MODS_MASS_DATA_TYPE,
                pickle.dumps(uuid4()), pickle.dumps(session.magic_label))
    cur.execute('INSERT INTO "metadata" VALUES (?, ?, ?, ?, ?, ?)', metadata)
    iden_lut.commit()
    iden_lut.close()

    for n, path in enumerate(session.iden_files): q.put((n, path))
    for i in range(num_of_threads): q.put(None)
    for pid in range(num_of_threads):
        processes.append(mp.Process(target=_worker, args=(pid, q, valid_iden_count, total_iden_count,
                                                          finish_count, log_lock, merge_lock, exit_signal)))
    reporter = Thread(target=_reporter, args=(finish_count, num_files, log_lock, exit_signal))

    logging.info('......Start identification import with {} subprocesses......'.format(num_of_threads))

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
        _error_clean(q)
        raise

    # in case that subprocesses were exited with errors
    if 1 in [p.exitcode for p in processes]:
        err_msg = '\nSubprocesses exited with abnormal exitcode.'
        logging.error(err_msg)
        logging.debug('Emptying queue')
        _error_clean(q)
        raise mp.ProcessError(err_msg)

    if valid_iden_count.value == 0:
        err_msg = '\nNo valid identification has been found.'
        logging.error(err_msg)
        raise Exception(err_msg)

    logging.info('......Finish identification import, found {} valid identifications (out of {} in total)......'
                 .format(valid_iden_count.value, total_iden_count.value))

    session.config.id_finished.value = True
    logging.debug('Mounting identification lookup table')
    session.mount_identification_lut()
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
                    logging.info('Import progress: {}% ({} out of {})' \
                                 .format(progress, finish_count.value, num_files))
        except KeyboardInterrupt:
            with log_lock:
                logging.debug('Reporter thread received KeyboardInterrupt, exits now.')
            break
    return


def _error_clean(q):
    while q.qsize() > 0:
        q.get()
    return


def _worker(pid, q, valid_iden_count, total_iden_count, finish_count, log_lock, merge_lock, exit_signal):
    try:
        logging_setup()
        with log_lock:
            logging.debug('Identification import subprocess {} started.'.format(pid))

        table_path = Path.cwd().joinpath(session.name + FILE_EXTENSION_IDEN_LUT).absolute()
        conn = sqlite3.connect(str(table_path))
        cur = conn.cursor()

        while True:
            if exit_signal.value:
                with log_lock:
                    logging.debug('Subprocess {}: Received exit signal, exits now.'.format(pid))
                break
            item = q.get()
            if item:
                file_id, path = item
                try:
                    iden_lut, valid, total = pepxml.import_identification(path, log_lock=log_lock)
                except Exception:
                    if ignore_errors:
                        wrn_msg = 'Subprocess {}: Failed to parse identification file: {}'.format(pid, path)
                        with log_lock:
                            logging.warning(wrn_msg)
                        continue
                    else:
                        raise
                if valid != 0:
                    with merge_lock:
                        for f in iden_lut.keys():
                            cur.execute('''
                                CREATE TABLE IF NOT EXISTS "{}" (
                                "native_id" INTEGER UNIQUE,
                                "peptide" TEXT,
                                "charge" INTEGER,
                                "probability" REAL,
                                "source" TEXT,
                                "is_decoy" INTEGER,
                                "prev_aa" TEXT,
                                "next_aa" TEXT,
                                "mods_pos" BLOB,
                                "mods_mass" BLOB,
                                "nterm_mod" REAL,
                                "cterm_mod" REAL,
                                "iden_file_id" INTEGER,
                                "l_offset" INTEGER,
                                "r_offset" INTEGER
                                )'''.format(f))
                            conn.commit()
                        iden_lut = _dict_to_list(iden_lut, file_id)
                        for f in iden_lut.keys():
                            cur.execute('BEGIN TRANSACTION')
                            for id_ in iden_lut[f]:
                                try:
                                    cur.execute('INSERT INTO "{}" VALUES ({})'.format(f, ', '.join(['?'] * 15)), id_)
                                except sqlite3.IntegrityError:
                                    if ignore_errors:
                                        wrn_msg = 'Duplicated identification: scan {} in MS file {}.'.format(id_[0], f)
                                        with log_lock:
                                            logging.warning(wrn_msg)
                                    else:
                                        err_msg = '\nDuplicated identification: scan {} in MS file {}.'\
                                                  '\nIdentification file: {}'.format(id_[0], f, path)
                                        with log_lock:
                                            logging.error(err_msg)
                            conn.commit()
                        valid_iden_count.value += valid
                        total_iden_count.value += total
                else:
                    wrn_msg = 'Subprocess {}: Found 0 valid identification in file: {}'.format(pid, path)
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


def _dict_to_list(iden_lut, iden_file_id):
    for f in iden_lut.keys():
        temp = []
        for id_ in iden_lut[f].keys():
            iden = iden_lut[f][id_]
            temp.append((id_, iden.peptide, iden.charge, iden.probability, iden.source, int(iden.is_decoy),
                         iden.prev_aa, iden.next_aa, iden.mods_pos.tostring(), iden.mods_mass.tostring(),
                         iden.nterm_mod, iden.cterm_mod, iden_file_id, iden.l_offset, iden.r_offset))
        iden_lut[f] = temp
    return iden_lut


def _refresh_session():
    global ignore_errors, num_of_threads
    ignore_errors = session.flags.ignore_errors
    num_of_threads = session.config.id_num_of_threads.value
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
