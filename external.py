# -*- coding: utf-8 -*-
"""
Created on 02:03:20 04/11/2019

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
import dev

from envr.session import get_session
from share.identification_lut import Identification
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
try:
    session = get_session()

    iid = session.internal_index
    ide = session.iden_lut
    rks = session.ranked_spectra

    identifications_dict = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def _get_all_identifications_python_dict():
    ide = session.iden_lut
    main = {}
    for f in session.ms_exp_files:
        sub = {}
        if ide.cursor.execute('SELECT name FROM sqlite_master WHERE type="table" AND name="{}"'.format(f.stem))\
                     .fetchone() is not None:
            for i in ide.cursor.execute('SELECT * FROM "{}"'.format(f.stem)).fetchall():
                mapped = _sqlite_to_iden_mapping(i)
                sub[mapped.native_id] = mapped.to_tpp_string_integer()
        main[f.stem] = sub

    return main


def _sqlite_to_iden_mapping(query):
    identification = Identification()
    identification.native_id = query[0]
    identification.peptide = query[1]
    identification.charge = query[2]
    identification.probability = query[3]
    identification.source = query[4]
    identification.is_decoy = query[5]
    identification.prev_aa = query[6]
    identification.next_aa = query[7]
    identification.mods_pos = np.frombuffer(query[8], dtype=ID_MODS_POS_DATA_TYPE)
    identification.mods_mass = np.frombuffer(query[9], dtype=ID_MODS_MASS_DATA_TYPE)
    identification.nterm_mod = query[10]
    identification.cterm_mod = query[11]
    identification.iden_file_id = query[12]
    identification.l_offset = query[13]
    identification.r_offset = query[14]
    return identification


def peptide_entropy(full_clusters_set):

    if identifications_dict is None:
        logging.info("Converting identification sqlite database to Python dictionary...")
        identifications_dict = _get_all_identifications_python_dict()

    import math

    ##Equation:
    ## p_ = [# of spectra id as peptide i in cluster C]/[Total # of spectra id as peptide i]
    ##peptide entropy = -SUM_clusterC{p_log2(p_)}

    temp_dict = {}  # for temporary storage, clear after every cluster to release space.
    pep_entropy_dict = {}  # true dictionary, should be ready for csv-writing.

    for cluster_id in session.clusters.cursor.execute('SELECT "cluster_id" FROM "clusters"').fetchmany():
        g = dev.get_graph(cluster_id)[-1]

        for internal_id in g.vp['iid'].a:
            entry = iid[internal_id]
            v_pep = identifications_dict.get(entry.get_file_path().stem)
            v_pep = v_pep.get(entry.native_id) if v_pep is not None else None
            if v_pep is not None:
                if v_pep not in temp_dict:
                    temp_dict[v_pep] = 1
                else:
                    temp_dict[v_pep] += 1

        for key in temp_dict.keys():  # for-loop peptide entropy
            if key not in pep_entropy_dict:
                p_ = temp_dict[key] / cluster.num_vertices()
                pep_entropy_dict[key] = -(p_ * math.log2(p_))
            else:
                pep_entropy_dict[key] += -(p_ * math.log2(p_))
    return pep_entropy_dict


def compute_external_indices():
    import math

    global identifications_dict
    if identifications_dict is None:
        logging.info("Converting identification sqlite database to Python dictionary...")
        identifications_dict = _get_all_identifications_python_dict()

    def _cluster_entropy_calculator_withoutUN(pep_countdict): #internal function to compute cluster entropy
        entropy = 0
        if sum(pep_countdict.values()) == 0:
            pass
        else:
            denominator = sum(pep_countdict.values())
            for key in pep_countdict.keys():
                numerator = pep_countdict[key]
                probability = numerator / denominator
                entropy = entropy - probability * math.log2(probability)
        return entropy

    overall_peptide_entropy = 0
    sum_clusters_entropy = 0
    num_pep_dict = {}

    clusters_set = []
    for cluster_id in session.clusters.cursor.execute('SELECT "cluster_id" FROM "clusters"').fetchmany():
        g = dev.get_graph(cluster_id)[-1]
        if g.num_vertices() == 0: continue
        clusters_set.append(g)

    for g in clusters_set:
        num_id_vertices = 0
        clus_vert_dict = {}

        for internal_id in g.vp['iid'].a:
            entry = iid[internal_id]
            v_pep = identifications_dict.get(entry.get_file_path().stem)
            v_pep = v_pep.get(entry.native_id) if v_pep is not None else None
            if v_pep is not None:
                num_id_vertices += 1
                if v_pep not in temp_dict:
                    clus_vert_dict[v_pep] = 1
                else:
                    clus_vert_dict[v_pep] += 1
                if v_pep not in num_pep_dict:
                    num_pep_dict[v_pep] = 1
                else:
                    num_pep_dict[v_pep] += 1
        cluster_entropy = _cluster_entropy_calculator_withoutUN(clus_vert_dict)
        sum_clusters_entropy += cluster_entropy * num_id_vertices

    if len(num_pep_dict) == 0: #exit function if nothing is marked by prior for-loop.
        print("No valid peptide identification found, program exited!")
        return nan, nan

    overall_cluster_entropy = sum_clusters_entropy / sum(num_pep_dict.values())  # ID mixedness (smaller the better)

    for peptide in num_pep_dict.keys(): #for-loop peptide entropy: "splitness"
        peptides_entropy = 0
        for g in clusters_set:
            count_ = 0
            for internal_id in g.vp['iid'].a:
                entry = iid[internal_id]
                v_pep = identifications_dict.get(entry.get_file_path().stem)
                v_pep = v_pep.get(entry.native_id) if v_pep is not None else None
                if v_pep == peptide:
                    count_ += 1
            if count_ > 0:
                probability_ = count_ / num_pep_dict[peptide]
                peptides_entropy += - probability_ * math.log2(probability_)
        # print(peptides_entropy)
        overall_peptide_entropy += peptides_entropy * num_pep_dict[peptide] / sum(num_pep_dict.values()) #ID splitness (smaller the better)
    print("Overall cluster entropy = {}\nOverall peptide entropy = {}".format(overall_cluster_entropy,overall_peptide_entropy))
    return overall_cluster_entropy, overall_peptide_entropy
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
