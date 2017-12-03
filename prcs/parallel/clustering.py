# -*- coding: utf-8 -*-
"""
Created on 12:15:06 17/11/2017

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
from pathlib import Path
import ctypes
from uuid import uuid4
import math

from envr.session import get_session
from prcs.parallel.clustering_gpu import clustering_gpu
from prcs.parallel.clustering_cpu import clustering_cpu
from prcs.parallel.find_cluster import _union_find, _find_parent, _find_belonging, _argsort
from prcs.parallel.graph_making import make_graphs
from share.misc import _clean_temporary
from property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    ignore_errors = None
    block_dimensions = None
    use_gpu = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
class _ClientState:

    def __int__(self):
        self.row_id = 0
        self.column_id = 0
        self.y_offset = 0
        self.x_offset = 0
        return


class _Dispatcher:

    def __init__(self, length, block_dimensions):
        self.lock = mp.Lock()
        self.current_row = mp.Value(ctypes.c_uint64, 0)
        self.length = length
        self.block_height, self.block_width = block_dimensions
        self.processes = {}
        return

    def connect(self, pid, tid):
        self.processes[(pid, tid)] = _ClientState()
        return

    def iterate(self, pid, tid):
        client = self.processes[(pid, tid)]
        self.next_row(pid, tid)
        while True:
            if client.x_offset >= self.length: self.next_row(pid, tid)
            if client.y_offset >= self.length: break
            y_range = (client.y_offset, min([self.length, (client.y_offset + self.block_height)]))
            x_range = (client.x_offset, min([self.length, (client.x_offset + self.block_width)]))
            client.x_offset += self.block_width
            yield client.row_id, client.column_id, (y_range, x_range)
            client.column_id += 1
        return

    def next_row(self, pid, tid):
        client = self.processes[(pid, tid)]
        with self.lock:
            client.row_id = self.current_row.value
            client.column_id = 0
            client.y_offset = client.x_offset = self.current_row.value * self.block_height
            self.current_row.value += 1
        return

    def now(self):
        return self.current_row.value - 1
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def clustering():
    _refresh_session()
    if session.ranked_spectra is None:
        err_msg = '\nNo ranked spectra mounted, cannot proceed.'
        logging.error(err_msg)
        raise Exception(err_msg)

    dispatcher = _Dispatcher(len(session.internal_index), block_dimensions)
    temp_storage = _new_temp_storage()
    gpu_failed = False

    if use_gpu:
        try:
            total_edge_count = clustering_gpu(dispatcher, temp_storage)
        except Exception:
            if not session.flags.force:
                err_msg = '\nUnable to do clustering with GPU.'
                logging.error(err_msg)
                raise
            else:
                gpu_failed = True

        if gpu_failed and session.flags.force:
            wrn_msg = 'Unable to do clustering with GPU. Degrades to CPU.'
            logging.warning(wrn_msg)
            dispatcher = _Dispatcher(len(session.internal_index), block_dimensions)
            temp_storage = _new_temp_storage()
            total_edge_count = _clustering_cpu(dispatcher, temp_storage)
    else:
        total_edge_count = _clustering_cpu(dispatcher, temp_storage)

    if total_edge_count.value <= 0:
        err_msg = '\nClustering produced no edge.'
        logging.error(err_msg)
        raise Exception(err_msg)
    logging.info('......Finish similarity computation, start graph formation......')

    try:
        num_of_clusters = _edge_list_to_graphs(temp_storage)
    except Exception:
        err_msg = '\nUnable to transform edge list to graphs.'
        logging.error(err_msg)
        _clean_temporary(temp_storage)
        raise

    logging.info('......Finish clustering, formed {} clusters......'.format(num_of_clusters))
    _clean_temporary(temp_storage)

    session.config.cg_finished.value = True
    logging.debug('Mounting clusters file')
    session.mount_clusters()
    return


def _clustering_cpu(dispatcher, temp_storage):
    try:
        total_edge_count = clustering_cpu(dispatcher, temp_storage)
    except Exception:
        err_msg = '\nUnable to do clustering with CPU.'
        logging.error(err_msg)
        raise
    return total_edge_count


def _new_temp_storage():
    temp_storage = Path.cwd().joinpath(str(uuid4()))
    temp_storage.mkdir()
    temp_storage = temp_storage.resolve()
    temp_storage.joinpath('edg').touch()
    temp_storage.joinpath('dps').touch()
    return temp_storage


# def _ref_edge_list_to_graphs(temp_storage):
#     edg_path = temp_storage.joinpath('edg')
#     edge_list_length = int(temp_storage.joinpath('dps').stat().st_size / 4)
#     edge = np.memmap(str(edg_path), dtype=CG_EDGE_DATA_TYPE, mode='c', shape=(edge_list_length, 2))
#     g = gt.Graph(directed=False)
#     g.add_edge_list(edge)
#     comp, hist = gt.label_components(g)
#     clusters = [[] for i in range(len(hist))]
#     for v in g.vertices():
#         clusters[comp[v]].append(int(v))
#     return clusters


def _edge_list_to_graphs(temp_storage):
    edg_path = temp_storage.joinpath('edg')
    dps_path = temp_storage.joinpath('dps')
    edge_list_length = int(temp_storage.joinpath('dps').stat().st_size / 4)
    logging.debug('The edge list file contains {} edges.'.format(edge_list_length))
    # somehow cython only accepts writable array
    # mode 'c' ensures that no change of actual data in disk will be made
    num_of_nodes = len(session.internal_index)
    edge = np.memmap(str(edg_path), dtype=CG_EDGE_DATA_TYPE, mode='r+', shape=(edge_list_length, 2))
    dps = np.memmap(str(dps_path), dtype=CG_DOT_PRODUCT_DATA_TYPE, mode='r+', shape=edge_list_length)
    par_path = temp_storage.joinpath('par')
    rnk_path = temp_storage.joinpath('rnk')
    parent = np.memmap(str(par_path), dtype=FG_PARENT_DATA_TYPE, mode='w+', shape=num_of_nodes)
    rank = np.memmap(str(rnk_path), dtype=FG_RANK_DATA_TYPE, mode='w+', shape=num_of_nodes)

    chunk_size = int(FG_DATA_COPY_CHUNK_SIZE / DATA_TYPE_SIZE[FG_PARENT_DATA_TYPE])
    # fill range chunk by chunk
    for start in range(0, math.ceil(num_of_nodes/chunk_size)):
        end = min((start + chunk_size), num_of_nodes)
        parent[start:end] = np.arange(start, end, dtype=FG_PARENT_DATA_TYPE)[:]
    # fill zeros
    for start in range(0, math.ceil(num_of_nodes/chunk_size)):
        end = min((start + chunk_size), num_of_nodes)
        rank[start:end] = np.zeros(end-start, dtype=FG_RANK_DATA_TYPE)[:]
    logging.debug('Finding connected components.')

    _union_find(edge, parent, rank)

    npr_path = temp_storage.joinpath('npr')
    new_parent = np.memmap(str(npr_path), dtype=FG_PARENT_DATA_TYPE, mode='w+', shape=num_of_nodes)
    _find_parent(parent, new_parent)
    del parent, rank
    par_path.unlink()
    rnk_path.unlink()

    blg_path = temp_storage.joinpath('blg')
    belonging = np.memmap(str(blg_path), dtype=FG_PARENT_DATA_TYPE, mode='w+', shape=edge_list_length*2)
    logging.debug('Assigning cluster id to each edge.')
    _find_belonging(edge, belonging, new_parent)
    del new_parent
    npr_path.unlink()

    _argsort(belonging)
    idx_path = temp_storage.joinpath('idx')
    nbg_path = temp_storage.joinpath('nbg')
    idx = np.memmap(str(idx_path), dtype=FG_PARENT_DATA_TYPE, mode='w+', shape=edge_list_length)
    nbg = np.memmap(str(nbg_path), dtype=FG_PARENT_DATA_TYPE, mode='w+', shape=edge_list_length)
    idx[:] = belonging[::2]
    idx.flush()
    nbg[:] = belonging[1::2]
    nbg.flush()
    del belonging
    blg_path.unlink()
    edge[:] = edge[idx]
    edge.flush()
    dps[:] = dps[idx]
    dps.flush()
    del idx
    idx_path.unlink()
    num_of_clusters = make_graphs(temp_storage, edge_list_length)
    return num_of_clusters


def _refresh_session():
    global ignore_errors, block_dimensions, use_gpu
    ignore_errors = session.flags.ignore_errors
    block_dimensions = session.config.cg_block_dimensions.value
    use_gpu = session.config.gpu_use_gpu.value and not session.flags.use_cpu
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
