#cython: boundscheck=False, wraparound=False, nonecheck=False
# -*- coding: utf-8 -*-
"""
Created on 08:14:33 20/11/2017

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
from libc.stdlib cimport qsort
cimport numpy as np
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
ctypedef np.uint64_t CG_EDGE_DATA_TYPE
ctypedef np.float32_t CG_DOT_PRODUCT_DATA_TYPE
ctypedef CG_EDGE_DATA_TYPE FG_PARENT_DATA_TYPE
ctypedef np.uint32_t FG_RANK_DATA_TYPE
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
cdef struct Element:
    FG_PARENT_DATA_TYPE index
    FG_PARENT_DATA_TYPE value
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def _union_find(CG_EDGE_DATA_TYPE[:,:] edge, FG_PARENT_DATA_TYPE[:] parent, FG_RANK_DATA_TYPE[:] rank):

    cdef CG_EDGE_DATA_TYPE source, target, i, length
    length = edge.shape[0]

    for i in range(length):
        source = edge[i][0]
        target = edge[i][1]

        while True:
            if parent[source] == source:
                break
            else:
                source = parent[source]

        while True:
            if parent[target] == target:
                break
            else:
                target = parent[target]

        if source == target:
            continue
        if rank[source] < rank[target]:
            parent[source] = target
        elif rank[source] > rank[target]:
            parent[target] = source
        else:
            parent[target] = source
            rank[source] += 1

    return


def _find_parent(FG_PARENT_DATA_TYPE[:] parent, FG_PARENT_DATA_TYPE[:] new_parent):

    cdef FG_PARENT_DATA_TYPE i, temp, length
    length = parent.shape[0]

    for i in range(length):
        temp = i
        while True:
            if parent[temp] == temp:
                break
            else:
                temp = parent[temp]
        new_parent[i] = temp

    return


def _find_belonging(CG_EDGE_DATA_TYPE[:,:] edge, FG_PARENT_DATA_TYPE[:] belonging, FG_PARENT_DATA_TYPE[:] parent):

    cdef FG_PARENT_DATA_TYPE i, idx, length
    length = edge.shape[0]

    for i in range(0, length):
        idx = i * 2
        belonging[idx] = i
        belonging[idx+1] = parent[edge[i][0]]

    return


cdef int _compare(const void *a, const void *b) nogil:
    if (<Element *> a).value < (<Element *> b).value: return -1
    if (<Element *> a).value == (<Element *> b).value: return 0
    if (<Element *> a).value > (<Element *> b).value: return 1


def _argsort(np.ndarray[FG_PARENT_DATA_TYPE, ndim=1] belonging):
    cdef Element *belonging_view = <Element *> belonging.data
    qsort(<void *> belonging_view, int(belonging.shape[0]/2), sizeof(Element), &_compare)
    return


def _unique(np.ndarray[FG_PARENT_DATA_TYPE, ndim=1] belonging, FG_PARENT_DATA_TYPE[:] unique_index):
    if belonging.shape[0] == 0:
        return

    cdef Element *belonging_view = <Element *> belonging.data
    cdef FG_PARENT_DATA_TYPE limit = int(belonging.shape[0]/2)
    cdef FG_PARENT_DATA_TYPE count = 0
    cdef FG_PARENT_DATA_TYPE cursor = 0
    cdef FG_PARENT_DATA_TYPE previous = belonging_view[0].index

    while True:
        if cursor == limit:
            unique_index[count] = cursor
            count += 1
            break
        if belonging_view[cursor].index != previous:
            unique_index[count] = cursor
            count += 1
            previous = belonging_view[cursor].index
        cursor += 1

    return count


def _scan(FG_PARENT_DATA_TYPE[:] belonging, FG_PARENT_DATA_TYPE start):

    if start >= belonging.shape[0]:
        return start
    cdef FG_PARENT_DATA_TYPE id_ = belonging[start]

    while True:
        if start == belonging.shape[0] or belonging[start] != id_:
            break
        start += 1

    return start
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
