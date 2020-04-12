# -*- coding: utf-8 -*-
"""
Created on 02:42:35 12/09/2017

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

from envr.session import get_session
from share.misc import extract_header
from property import *
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
class RankedSpectra:

    def __init__(self):
        self.precursor_mass = None
        self.mz = None
        self.intensity = None
        self.length = None
        self.num_of_peaks = None
        self.magic_label = uuid4()
        self.session_label = None
        self.index_label = None
        return

    @staticmethod
    def mount_file(read_only=True):
        file = Path.cwd().joinpath(session.name + FILE_EXTENSION_RANKED_SPECTRA)
        try:
            file = file.resolve()
        except Exception:
            err_msg = '\nCould not resolve ranked spectra file: {}'.format(file)
            logging.error(err_msg)
            raise
        try:
            header = extract_header(file)
            header = pickle.loads(header)
            if header['file_type'] != HEADER_LABEL_RANKED_SPECTRA:
                err_msg = '\nTarget file is not a ranked spectra file. It is \"{}\".'.format(header['file_type'])
                logging.error(err_msg)
                raise Exception(err_msg)

            session_label = header['session_label']
            if integrity_check and session_label != session.magic_label:
                err_msg = '\nSession label does not match the current session.'\
                          '\nThese spectra were not produced from the current session.'
                logging.error(err_msg)
                raise Exception(err_msg)

            index_label = header['index_label']
            if integrity_check and index_label != session.internal_index.magic_label:
                err_msg = '\nIndex label does not match the current internal index.'\
                          '\nThese spectra were not produced from the current internal index.'
                logging.error(err_msg)
                raise Exception(err_msg)

            magic_label = header['magic_label']
            length = header['length']
            num_of_peaks = header['num_of_peaks']
            mz_data_type = header['mz_data_type']
            intensity_data_type = header['intensity_data_type']
            mz_offset = header['mz_offset']
            intensity_offset = header['intensity_offset']
        except Exception:
            err_msg = '\nIncorrect file type or corrupted file.'
            logging.error(err_msg)
            raise

        ranked_spectra = RankedSpectra()
        ranked_spectra.magic_label = magic_label
        ranked_spectra.session_label = session_label
        ranked_spectra.index_label = index_label
        ranked_spectra.length = length
        ranked_spectra.num_of_peaks = num_of_peaks
        ranked_spectra.mz = session.internal_index.precursor_mass

        mode = 'r' if read_only else 'r+'
        ranked_spectra.mz = np.memmap(str(file), offset=mz_offset, dtype=mz_data_type,
                                      mode=mode, shape=(length, num_of_peaks))
        ranked_spectra.intensity = np.memmap(str(file), offset=intensity_offset, dtype=intensity_data_type,
                                             mode=mode, shape=(length, num_of_peaks))
        ranked_spectra.precursor_mass = session.internal_index.precursor_mass
        return ranked_spectra

    def __repr__(self):
        return pformat(vars(self))
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
