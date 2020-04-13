# -*- coding: utf-8 -*-
"""
Created on 02:28:48 12/09/2017

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
from pprint import pformat

from ClusterSheep.envr.session import get_session
from ClusterSheep.share.misc import extract_header
from ClusterSheep.property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
    iden_lut = session.iden_lut
    integrity_check = session.config.gr_integrity_check.value
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
class Index(object):

    def __init__(self):
        self.precursor_mass = None
        self.precursor_charge = None
        self.file_id = None
        self.offset = None
        self.internal_id = None
        self.native_id = None
        self.length = None
        self.true_precursor_mass = None
        self.magic_label = uuid4()
        self.session_label = None
        return

    def __getitem__(self, item):
        if isinstance(item, int):
            new_entry = Entry()
            new_entry.precursor_mass = self.precursor_mass[item]
            new_entry.precursor_charge = self.precursor_charge[item]
            new_entry.file_id = self.file_id[item]
            new_entry.offset = self.offset[item]
            new_entry.internal_id = self.internal_id[item]
            new_entry.native_id = self.native_id[item]
            new_entry.format = Index.get_format(new_entry.file_id)
            return new_entry
        elif isinstance(item, slice):
            new_index = Index()
            new_index.precursor_mass = self.precursor_mass[item.start:item.stop].copy()
            new_index.precursor_charge = self.precursor_charge[item.start:item.stop].copy()
            new_index.file_id = self.file_id[item.start:item.stop].copy()
            new_index.offset = self.offset[item.start:item.stop].copy()
            new_index.internal_id = self.internal_id[item.start:item.stop].copy()
            new_index.native_id = self.native_id[item.start:item.stop].copy()
            new_index.length = item.stop - item.start
            return new_index
        else:
            err_msg = '\nInvalid index type: {}'.format(type(item))
            logging.error(err_msg)
            raise TypeError(err_msg)

    def __iter__(self):
        for i in range(self.length):
            yield self[i]

    def __len__(self):
        return self.length

    @staticmethod
    def get_format(file_id):
        return ''.join(session.ms_exp_files[file_id].suffixes)

    @staticmethod
    def mount_file(read_only=True):
        file = Path.cwd().joinpath(session.name + FILE_EXTENSION_INTERNAL_INDEX)
        try:
            file = file.resolve()
        except Exception:
            err_msg = '\nCould not resolve internal index file: {}'.format(file)
            logging.error(err_msg)
            raise
        try:
            header = extract_header(file)
            header = pickle.loads(header)
            if header['file_type'] != HEADER_LABEL_INTERNAL_INDEX:
                err_msg = '\nTarget file is not an index file. It is \"{}\".'.format(header['file_type'])
                logging.error(err_msg)
                raise Exception(err_msg)

            session_label = header['session_label']
            if integrity_check and session_label != session.magic_label:
                err_msg = '\nSession label does not match the current session.'\
                          '\nThis index was not produced from the current session.'
                logging.error(err_msg)
                raise Exception(err_msg)

            magic_label = header['magic_label']
            length = header['length']
            true_precursor_mass = header['true_precursor_mass']
            precursor_mass_data_type = header['precursor_mass_data_type']
            precursor_charge_data_type = header['precursor_charge_data_type']
            file_id_data_type = header['file_id_data_type']
            offset_data_type = header['offset_data_type']
            internal_id_data_type = header['internal_id_data_type']
            native_id_data_type = header['native_id_data_type']
            precursor_mass_offset = header['precursor_mass_offset']
            precursor_charge_offset = header['precursor_charge_offset']
            file_id_offset = header['file_id_offset']
            offset_offset = header['offset_offset']
            internal_id_offset = header['internal_id_offset']
            native_id_offset = header['native_id_offset']
        except Exception:
            err_msg = '\nIncorrect file type or corrupted file.'
            logging.error(err_msg)
            raise

        index = Index()
        index.magic_label = magic_label
        index.session_label = session_label
        index.length = length
        index.true_precursor_mass = true_precursor_mass

        mode = 'r' if read_only else 'r+'
        index.precursor_mass = np.memmap(str(file), offset=precursor_mass_offset,
                                         dtype=precursor_mass_data_type, mode=mode, shape=length)
        index.precursor_charge = np.memmap(str(file), offset=precursor_charge_offset,
                                           dtype=precursor_charge_data_type, mode=mode, shape=length)
        index.file_id = np.memmap(str(file), offset=file_id_offset,
                                  dtype=file_id_data_type, mode=mode, shape=length)
        index.offset = np.memmap(str(file), offset=offset_offset,
                                 dtype=offset_data_type, mode=mode, shape=(length, 2))
        index.internal_id = np.memmap(str(file), offset=internal_id_offset,
                                      dtype=internal_id_data_type, mode=mode, shape=length)
        index.native_id = np.memmap(str(file), offset=native_id_offset,
                                    dtype=native_id_data_type, mode=mode, shape=length)
        return index

    def __repr__(self):
        return pformat(vars(self))


class Entry(object):

    def __init__(self):
        self.precursor_mass = None
        self.precursor_charge = None
        self.file_id = None
        self.offset = None
        self.internal_id = None
        self.native_id = None
        self.format = None
        return

    def get_file_path(self):
        return session.ms_exp_files[self.file_id]

    def get_identification(self):
        if session.iden_lut is None:
            return None
        ms_exp_file_name = self.get_file_path().stem
        return session.iden_lut.get_identification(ms_exp_file_name, self.native_id)

    def get_spectrum(self):
        from share.spectrum import Spectrum
        return Spectrum(self)

    def s(self):
        return self.get_spectrum()

    def __repr__(self):
        return pformat(vars(self))
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
