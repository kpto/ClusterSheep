# -*- coding: utf-8 -*-
"""
Created on 17:54:29 23/11/2017

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
from uuid import uuid4
import pickle
from pathlib import Path
import sqlite3
from pprint import pformat

from envr.session import get_session
from property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    integrity_check = session.config.gr_integrity_check.value
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
class Clusters:

    def __init__(self):
        self.file_path = None
        self.magic_label = uuid4()
        self.session_label = None
        self.index_label = None
        self.connection = None
        self.cursor = None
        self.is_connected = False
        return

    @staticmethod
    def mount_file():
        file = Path.cwd().joinpath(session.name + FILE_EXTENSION_CLUSTERS)
        try:
            file = file.resolve()
        except Exception:
            err_msg = '\nCould not resolve clusters file: {}'.format(file)
            logging.error(err_msg)
            raise
        try:
            clusters = sqlite3.connect(str(file))
            cur = clusters.cursor()
            file_type = cur.execute('SELECT "file_type" FROM "metadata"').fetchone()[0]
            magic_label = pickle.loads(cur.execute('SELECT "magic_label" FROM "metadata"').fetchone()[0])
            session_label = pickle.loads(cur.execute('SELECT "session_label" FROM "metadata"').fetchone()[0])
            index_label = pickle.loads(cur.execute('SELECT "index_label" FROM "metadata"').fetchone()[0])

            if file_type != HEADER_LABEL_CLUSTERS:
                err_msg = '\nTarget file is not a cluster file. It is \"{}\".'.format(file_type)
                logging.error(err_msg)
                raise Exception(err_msg)

            if integrity_check and session_label != session.magic_label:
                err_msg = '\nSession label does not match the current session.'\
                          '\nThis cluster file was not produced from the current session.'
                logging.error(err_msg)
                raise Exception(err_msg)

            if integrity_check and index_label != session.internal_index.magic_label:
                err_msg = '\nIndex label does not match the current internal index.'\
                          '\nThis cluster file was not produced from the current internal index.'
                logging.error(err_msg)
                raise Exception(err_msg)

            clusters.close()
        except Exception:
            err_msg = '\nIncorrect file type or corrupted file.'
            logging.error(err_msg)
            clusters.close()
            raise

        clusters = Clusters()
        clusters.file_path = file
        clusters.magic_label = magic_label
        clusters.index_label = index_label
        clusters.session_label = session_label
        return clusters

    def connect(self):
        if self.is_connected: return
        self.connection = sqlite3.connect(str(self.file_path))
        self.cursor = self.connection.cursor()
        self.is_connected = True
        return

    def disconnect(self):
        if not self.is_connected: return
        self.connection.close()
        self.cursor = None
        self.connection = None
        self.is_connected = False
        return

    def __repr__(self):
        return pformat(vars(self))
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
