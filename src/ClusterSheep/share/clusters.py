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

from ClusterSheep.envr.session import get_session
from ClusterSheep.property import *
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

    def get_cluster(self, cluster_id):
        self.connect()
        row = self.cursor.execute('SELECT * FROM "clusters" WHERE "cluster_id"=?',
                                  (cluster_id,)).fetchone()
        if row is None:
            err_msg = '\nCluster {} does not exist.'.format(cluster_id)
            logging.error(err_msg)
        return Clusters._row_to_cluster(row[:-1], row[-1])

    def exists(self, cluster_id):
        self.connect()
        return self.cursor.execute('SELECT EXISTS(SELECT * FROM "clusters" WHERE "cluster_id"=?)',
                                   (cluster_id,)).fetchone()[0]

    def _get_cluster_no_binary(self, cluster_id):
        row = self.cursor.execute('SELECT "cluster_id", "num_nodes", "num_edges", "num_idens", ' +
                                  '"idens", "major_iden", "iden_ratio", "pre_mass_avg" FROM "clusters"' +
                                  'WHERE "cluster_id" = ?', (cluster_id,)).fetchone()
        return Clusters._row_to_cluster(row, self._get_graph_binary_provider(row[0]))

    def _get_graph_binary_provider(self, cluster_id):
        return lambda: self._provide_graph_binary(cluster_id)

    def _provide_graph_binary(self, cluster_id):
        cur = self.connection.cursor()
        return self.cursor.execute('SELECT "pickled" FROM "clusters" WHERE "cluster_id" = ?',
                                   (cluster_id,)).fetchone()[0]

    @staticmethod
    def _row_to_cluster(row, graph_binary):
        cluster = Cluster()
        cluster.id,\
        cluster.num_of_nodes,\
        cluster.num_of_edges,\
        cluster.num_of_identifications,\
        cluster._identifications,\
        cluster.major_identification,\
        cluster.identified_ratio,\
        cluster.average_precursor_mass = row
        cluster._graph = graph_binary
        return cluster

    def __iter__(self):
        self.cursor.execute('SELECT "cluster_id", "num_nodes", "num_edges", "num_idens", ' +
                            '"idens", "major_iden", "iden_ratio", "pre_mass_avg" FROM "clusters"')
        for row in self.cursor:
            graph_provider = self._get_graph_binary_provider(row[0])
            yield Clusters._row_to_cluster(row, graph_provider)

    def __repr__(self):
        return pformat(vars(self))


class Cluster:

    def __init__(self):
        self.id = None
        self.num_of_nodes = None
        self.num_of_edges = None
        self.num_of_identifications = None
        self.major_identification = None
        self.identified_ratio = None
        self.average_precursor_mass = None
        self._graph = None
        self._graph_loaded = False
        self._identifications = None
        self._identification_splitted = False
        return

    @property
    def graph(self):
        if not self._graph_loaded:
            # Value is a provider, get the binary first
            if callable(self._graph):
                self._graph = self._graph()
            self._graph = pickle.loads(self._graph)
            self._graph_loaded = True
        return self._graph

    @property
    def identifications(self):
        if not self._identification_splitted:
            if self.num_of_identifications is not None:
                if self._identifications is None:
                    self._identifications = []
                else:
                    self._identifications = self._identifications.split(';')
            self._identification_splitted = True
        return self._identifications

    def __repr__(self):
        return pformat({ 'id': self.id, 'num_of_nodes': self.num_of_nodes, 'num_of_edges': self.num_of_edges })
# ====END OF CLASS DEFINITION====

# ====BEGIN OF CODE====
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
