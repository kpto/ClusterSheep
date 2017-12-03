# -*- coding: utf-8 -*-
"""
Created on 13:38:30 10/09/2017

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
from threading import Lock
import mmap
import re

from envr.session import get_session
from share.identification_lut import Identification
from property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    ignore_errors = None
    prob_threshold = None
    include_decoy = None
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def import_identification(file, log_lock=Lock()):
    total_iden_count = 0
    valid_iden_count = 0
    fail_count = 0
    iden_lut = {}
    with file.open('rb') as fp:
        offset_list = _get_offset_list(fp, log_lock)
        for offset in offset_list:
            start_pos, end_pos = offset
            fp.seek(start_pos)
            block = fp.read(end_pos-start_pos)

            new_iden = Identification()
            temp_mods_pos = []
            temp_mods_mass = []

            try:
                if re.search(b'<interprophet_result', block):
                    new_iden.probability = float(re.search(b'<interprophet.+probability="([^"]+)', block).groups()[0])
                    new_iden.source = 'ip'
                else:
                    new_iden.probability = float(re.search(b'<peptideprophet.+probability="([^"]+)', block).groups()[0])
                    new_iden.source = 'pp'
                if new_iden.probability < prob_threshold:
                    total_iden_count += 1
                    continue
                if re.search(b'protein="(.+)"', block).groups()[0].startswith(b'DECOY'):
                    if include_decoy:
                        new_iden.is_decoy = True
                    else:
                        total_iden_count += 1
                        continue
                else:
                    new_iden.is_decoy = False
                ms_file_name = re.search(b'spectrum="([^.]+)', block).groups()[0].decode()
                native_id = int(re.search(b'start_scan="([^"]+)', block).groups()[0])
                new_iden.charge = int(re.search(b'assumed_charge="([^"]+)', block).groups()[0])
                new_iden.peptide = re.search(b'peptide="([^"]+)', block).groups()[0].decode()
                prev_aa = re.search(b'prev_aa="([^"]+)', block)
                if prev_aa: new_iden.prev_aa = prev_aa.groups()[0].decode()
                next_aa = re.search(b'next_aa="([^"]+)', block)
                if next_aa: new_iden.next_aa = next_aa.groups()[0].decode()

                nterm_mod = re.search(b'mod_nterm_mass="([^"]+)', block)
                if nterm_mod:
                    new_iden.nterm_mod = nterm_mod.groups()[0].decode()
                cterm_mod = re.search(b'mod_cterm_mass="([^"]+)', block)
                if cterm_mod:
                    new_iden.cterm_mod = cterm_mod.groups()[0].decode()
                for mod in re.findall(b'<mod_aminoacid_mass.+>', block):
                    temp_mods_pos.append(int(re.search(b'position="([^"]+)', mod).groups()[0]))
                    temp_mods_mass.append(float(re.search(b'mass="([^"]+)', mod).groups()[0]))
                new_iden.mods_pos = np.array(temp_mods_pos, dtype=ID_MODS_POS_DATA_TYPE)
                new_iden.mods_mass = np.array(temp_mods_mass, dtype=ID_MODS_MASS_DATA_TYPE)
                new_iden.l_offset = start_pos
                new_iden.r_offset = end_pos
            except Exception:
                if ignore_errors:
                    fail_count += 1
                    continue
                else:
                    err_msg = '\nUnable to get the necessary information of the identification located at {}.'\
                              '\nFile name: {}'\
                              '\nScan num: {}    Charge: {}'\
                              '\nProbability: {}    Peptide: {}    Is decoy: {}'\
                              '\nPrev aa: {}    Next aa: {}    Nterm mod: {}    Cterm mod: {}'\
                              '\nMods pos: {}    Mods mass: {}'\
                              .format(start_pos, new_iden.ms_file_name, new_iden.native_id, new_iden.charge,
                                      new_iden.probability, new_iden.peptide, new_iden.is_decoy, new_iden.prev_aa,
                                      new_iden.next_aa, new_iden.nterm_mod, new_iden.cterm_mod, new_iden.mods_pos,
                                      new_iden.mods_mass)
                    with log_lock:
                        logging.error(err_msg)
                    raise SyntaxError(err_msg)

            # add to table
            if ms_file_name not in iden_lut.keys():
                iden_lut[ms_file_name] = {}
            iden_lut[ms_file_name][native_id] = new_iden
            total_iden_count += 1
            valid_iden_count += 1

    if fail_count > 0:
        wrn_msg = '{} identifications were skipped while reading file \"{}\" due to reading errors.'\
                  .format(fail_count, file)
        with log_lock:
            logging.warning(wrn_msg)

    with log_lock:
        logging.debug('{} valid identifications (out of {}) were imported from file: {}'\
                      .format(valid_iden_count, total_iden_count, file))

    return iden_lut, valid_iden_count, total_iden_count


def _get_offset_list(fp, log_lock=Lock()):
    offset_list = []
    mm = mmap.mmap(fp.fileno(), 0, prot=mmap.PROT_READ)
    while True:
        start = mm.find(b'<spectrum_query')
        if start == -1:
            mm.close()
            break
        mm.seek(start)
        end = mm.find(b'</spectrum_query>')
        if end == -1:
            if not ignore_errors:
                err_msg = '\nCould not find the ending tag "</spectrum_query>" after the starting tag.'
                with log_lock:
                    logging.error(err_msg)
            raise(err_msg)
        offset_list.append((start, end))
        mm.seek(end)
    return offset_list


def _refresh_session():
    global ignore_errors, prob_threshold, include_decoy
    ignore_errors = session.flags.ignore_errors
    prob_threshold = session.config.id_prob_threshold.value
    include_decoy = session.config.id_include_decoy.value
    return


_refresh_session()
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
