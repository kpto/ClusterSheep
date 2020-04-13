# -*- coding: utf-8 -*-
"""
Created on 16:35:20 05/12/2017

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

from ClusterSheep.envr.session import get_session
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
try:
    session = get_session()
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def get_graph(cluster_id=None):
    import pickle
    query = session.clusters.cursor.execute('SELECT * FROM "clusters" WHERE "cluster_id"=?',
                                            (cluster_id,)).fetchone()
    return query[:-1] + (pickle.loads(query[-1]),)


def clustering():
    precursor_tolerance = session.config.cg_precursor_tolerance.value
    dot_product_threshold = session.config.cg_dot_product_threshold.value
    index_length = len(session.internal_index)
    dps = []
    edg = []
    for i, ei in enumerate(session.internal_index):
        print(i, end='\r')
        for j, ej in enumerate(session.internal_index[i+1:index_length]):
            if abs(ei.precursor_mass - ej.precursor_mass) <= precursor_tolerance:
                dp = ei.get_spectrum().verificative_ranked_dp(ej.get_spectrum())
                if dp > dot_product_threshold:
                    edg.append((i, j+i+1))
                    dps.append(dp)
    print()
    import pprint
    pprint.pprint(edg)
    pprint.pprint(dps)
    import graph_tool.all as gt
    import numpy as np
    g = gt.Graph(directed=False)
    g.add_edge_list(edg)
    g.vp['iid'] = g.new_vp('int')
    for i in range(g.num_vertices()):
        g.vp['iid'][i] = i
    g.ep['dps'] = g.new_ep('float')
    g.ep['dps'].a = np.array(dps, dtype=np.float32)
    comp, hist = gt.label_components(g)
    graphs = []
    for i in range(len(hist)):
        graphs.append(gt.Graph(gt.GraphView(g, vfilt=(comp.a == i)), directed=False, prune=True))
    return graphs


def mount_edge_list(path):
    import numpy as np
    length = int(path.joinpath('dps').stat().st_size / 4)
    edge = np.memmap(str(path.joinpath('edg')), dtype=np.uint64, mode='c', shape=(length, 2))
    dps = np.memmap(str(path.joinpath('dps')), dtype=np.float32, mode='c', shape=length)
    return edge, dps


def dp_spectrum_pair(idx_a, idx_b):
    from ClusterSheep.share.spectrum import Spectrum
    index = session.internal_index
    s_a = Spectrum(index[idx_a]).clip(copy=False).remove_precursor(copy=False).rank_transform(copy=False)
    s_b = Spectrum(index[idx_b]).clip(copy=False).remove_precursor(copy=False).rank_transform(copy=False)
    return s_a.ranked_dot_product(s_b)


def check_dot_product_correctness(edge, dps, sampling_size=1000, rel_tol=1e-06):
    import random, math
    length = len(dps)
    try:
        for i in range(sampling_size):
            choice = random.randrange(length)
            record = dps[choice]
            idx_a, idx_b = int(edge[choice][0]), int(edge[choice][1])
            calc = dp_spectrum_pair(idx_a, idx_b)
            if not math.isclose(record, calc, rel_tol=rel_tol):
                print('Incorrect value at {}, spectrum pair {}-{}.'.format(choice, idx_a, idx_b))
                print('Record: {}    Calc: {}'.format(record, calc))
                break
    except Exception:
        print('Exception occurs, choice = {}'.format(choice))
        raise
    print('Check completed, no error.')
    return


def check_precursor_tolerance_correctness(edge, dps, sampling_size=1000):
    import random
    length = len(dps)
    index = session.internal_index
    precursor_tolerance = session.config.cg_precursor_tolerance.value
    try:
        for i in range(sampling_size):
            choice = random.randrange(length)
            idx_a, idx_b = int(edge[choice][0]), int(edge[choice][1])
            difference = abs(index.precursor_mass[idx_a]-index.precursor_mass[idx_b])
            if not difference <= precursor_tolerance:
                print('Incorrect value at {}, spectrum pair {}-{}.'.format(choice, idx_a, idx_b))
                print('Set tolerance: {}    Difference: {}'.format(precursor_tolerance, difference))
                break
    except Exception:
        print('Exception occurs, choice = {}'.format(choice))
        raise
    print('Check completed, no error.')
    return


def check_dot_product_threshold_correctness(edge, dps, sampling_size=1000, find_close_tol=1e-06):
    import random, math
    length = len(dps)
    dot_product_threshold = session.config.cg_dot_product_threshold.value
    found_close = []
    try:
        for i in range(sampling_size):
            choice = random.randrange(length)
            dp = dps[choice]
            if not dp > dot_product_threshold:
                idx_a, idx_b = int(edge[choice][0]), int(edge[choice][1])
                print('Incorrect value at {}, spectrum pair {}-{}.'.format(choice, idx_a, idx_b))
                print('Set threshold: {}    Dot product: {}'.format(dot_product_threshold, dp))
                break
            if math.isclose(dp, dot_product_threshold, rel_tol=find_close_tol):
                found_close.append(dp)
    except Exception:
        print('Exception occurs, choice = {}'.format(choice))
        raise
    print('Check completed, no error.')
    if found_close:
        print('Found close values:')
        print(found_close)
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
