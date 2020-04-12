# -*- coding: utf-8 -*-
"""
Created on 22:08:57 13/11/2017

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
class IdentificationLUT:

    def __init__(self):
        self.file_path = None
        self.mods_pos_data_type = None
        self.mods_mass_data_type = None
        self.files = None
        self.magic_label = uuid4()
        self.session_label = None
        self.connection = None
        self.cursor = None
        self.is_connected = False
        return

    @staticmethod
    def mount_file():
        file = Path.cwd().joinpath(session.name + FILE_EXTENSION_IDEN_LUT)
        try:
            file = file.resolve()
        except Exception:
            err_msg = '\nCould not resolve identification lookup table file: {}'.format(file)
            logging.error(err_msg)
            raise
        try:
            iden_lut = sqlite3.connect(str(file))
            cur = iden_lut.cursor()
            file_type = cur.execute('SELECT "file_type" FROM "metadata"').fetchone()[0]
            mods_pos_data_type = cur.execute('SELECT "mods_pos_data_type" FROM "metadata"').fetchone()[0]
            mods_mass_data_type = cur.execute('SELECT "mods_mass_data_type" FROM "metadata"').fetchone()[0]
            magic_label = pickle.loads(cur.execute('SELECT "magic_label" FROM "metadata"').fetchone()[0])
            session_label = pickle.loads(cur.execute('SELECT "session_label" FROM "metadata"').fetchone()[0])

            if file_type != HEADER_LABEL_IDEN_LUT:
                err_msg = '\nTarget file is not an identification lookup table file. It is \"{}\".'.format(file_type)
                logging.error(err_msg)
                raise Exception(err_msg)

            if integrity_check and session_label != session.magic_label:
                err_msg = '\nSession label does not match the current session.'\
                          '\nThis identification lookup file was not produced from the current session.'
                logging.error(err_msg)
                raise Exception(err_msg)

            iden_lut.close()
        except Exception:
            err_msg = '\nIncorrect file type or corrupted file.'
            logging.error(err_msg)
            iden_lut.close()
            raise

        iden_lut = IdentificationLUT()
        iden_lut.file_path = file
        iden_lut.mods_pos_data_type = mods_pos_data_type
        iden_lut.mods_mass_data_type = mods_mass_data_type
        iden_lut.magic_label = magic_label
        iden_lut.session_label = session_label
        return iden_lut

    def connect(self):
        if self.is_connected: return
        self.connection = sqlite3.connect(str(self.file_path))
        self.cursor = self.connection.cursor()
        self.cursor.execute('SELECT "name" FROM "sqlite_master" WHERE "type" = "table"')
        self.files = {f[0] for f in self.cursor.fetchall()}
        self.is_connected = True
        return

    def disconnect(self):
        if not self.is_connected: return
        self.connection.close()
        self.cursor = None
        self.connection = None
        self.files = None
        self.is_connected = False
        return

    def get_identification(self, ms_exp_file_name, native_id):
        if ms_exp_file_name not in self.files: return None
        query = self.cursor.execute('SELECT * FROM "{}" WHERE "native_id"="{}"'.format(ms_exp_file_name, native_id)).fetchone()
        if query is None: return None
        identification = Identification()
        identification.peptide = query[1]
        identification.charge = query[2]
        identification.probability = query[3]
        identification.source = query[4]
        identification.is_decoy = query[5]
        identification.prev_aa = query[6]
        identification.next_aa = query[7]
        identification.mods_pos = np.frombuffer(query[8], dtype=self.mods_pos_data_type)
        identification.mods_mass = np.frombuffer(query[9], dtype=self.mods_mass_data_type)
        identification.nterm_mod = query[10]
        identification.cterm_mod = query[11]
        identification.iden_file_id = query[12]
        identification.l_offset = query[13]
        identification.r_offset = query[14]
        return identification

    def __repr__(self):
        return pformat(vars(self))


class Identification:

    def __init__(self):
        self.peptide = None
        self.charge = None
        self.probability = None
        self.source = None
        self.is_decoy = None
        self.prev_aa = None
        self.next_aa = None
        self.mods_pos = None
        self.mods_mass = None
        self.nterm_mod = None
        self.cterm_mod = None
        self.iden_file_id = None
        self.l_offset = None
        self.r_offset = None
        return

    def to_tpp_string(self):
        idx = np.argsort(self.mods_pos)
        mods_pos = self.mods_pos[idx]
        mods_mass = self.mods_mass[idx]
        full_string = self.peptide
        for i, pos in enumerate(mods_pos):
            offset = i * 4
            pos += offset
            full_string = full_string[:pos] + '[{}]' + full_string[pos:]
        full_string = full_string.format(*mods_mass)
        if self.nterm_mod: full_string = 'n[' + str(self.nterm_mod) + ']' + full_string
        if self.cterm_mod: full_string += 'c[' + str(self.cterm_mod) + ']'
        return full_string

    def to_tpp_string_integer(self):
        idx = np.argsort(self.mods_pos)
        mods_pos = self.mods_pos[idx]
        mods_mass = self.mods_mass[idx]
        full_string = self.peptide
        for i, pos in enumerate(mods_pos):
            offset = i * 4
            pos += offset
            full_string = full_string[:pos] + '[{}]' + full_string[pos:]
        full_string = full_string.format(*mods_mass.astype(np.int32))
        if self.nterm_mod: full_string = 'n[' + str(int(self.nterm_mod)) + ']' + full_string
        if self.cterm_mod: full_string += 'c[' + str(int(self.cterm_mod)) + ']'
        return full_string

    def to_string(self, tpp_string=None):
        full_string = tpp_string if tpp_string else self.to_tpp_string()
        prev_aa = self.prev_aa if self.prev_aa else ''
        next_aa = self.next_aa if self.next_aa else ''
        full_string = prev_aa + '.' + full_string + '.' + next_aa + '/' + str(self.charge)
        if self.is_decoy: full_string = 'DECOY_' + full_string
        full_string += ' (' + str(self.probability) + ')'
        return full_string

    def get_raw_text(self):
        with session.iden_files[self.iden_file_id].open('rb') as fp:
            fp.seek(self.l_offset)
            raw_text = fp.read(self.r_offset - self.l_offset)
        return raw_text

    def __repr__(self):
        return pformat(vars(self))
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
