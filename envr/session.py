# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 06:14:03 2017

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
from pathlib import Path
from datetime import datetime
import pickle
from uuid import uuid4
from pprint import pformat
import os

from envr.configuration import Configuration
from property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
_session = None
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
class Session:

    def __init__(self):
        time = datetime.now()
        self.name = NAME.lower() + '_' + time.isoformat()
        self.creation_time = time
        self.magic_label = uuid4()
        self.flags = None
        self.config = None
        self.ms_exp_files = None
        self.iden_files = None
        self.internal_index = None
        self.ranked_spectra = None
        self.iden_lut = None
        self.clusters = None
        self.console_hist = ''
        self.viewer_hist = ''
        return

    def __str__(self):
        string = ''
        string += '\nSession name: ' + self.name
        string += '\nCreation time: ' + self.creation_time.isoformat()
        string += '\nMagic label: ' + str(self.magic_label)
        string += '\n\nConfiguration:\n' + str(self.config)
        string += '\n\nMS exp files:\n' + str(self.ms_exp_files)
        string += '\n\nIdentification files:\n' + str(self.iden_files)
        return string

    def __repr__(self):
        return pformat(vars(self))

    def mount_internal_index(self):
        from share.internal_index import Index
        if type(self.internal_index) is Index: return
        self.internal_index = Index.mount_file()
        return

    def mount_identification_lut(self):
        from share.identification_lut import IdentificationLUT
        if type(self.iden_lut) is IdentificationLUT: return
        self.iden_lut = IdentificationLUT.mount_file()
        return

    def mount_ranked_spectra(self):
        from share.ranked_spectra import RankedSpectra
        if type(self.ranked_spectra) is RankedSpectra: return
        self.ranked_spectra = RankedSpectra.mount_file()
        return

    def mount_clusters(self):
        from share.clusters import Clusters
        if type(self.clusters) is Clusters: return
        self.clusters = Clusters.mount_file()
        return

    def save_session(self):
        logging.info('Saving session.')
        session_file = Path.cwd().joinpath(_session.name + FILE_EXTENSION_SESSION)
        if session_file.exists():
            session_file.unlink()
            session_file.touch()
        vars_ = dict()
        vars_['version'] = VERSION
        vars_['file_type'] = HEADER_LABEL_SESSION
        vars_['name'] = self.name
        vars_['creation_time'] = self.creation_time
        vars_['magic_label'] = self.magic_label
        vars_['config'] = self.config
        vars_['ms_exp_files'] = self.ms_exp_files
        vars_['iden_files'] = self.iden_files
        vars_['console_hist'] = self.console_hist
        vars_['viewer_hist'] = self.viewer_hist
        with session_file.open('wb') as fp:
            pickle.dump(vars_, fp)
        logging.debug('Session saved.')
        return

    @staticmethod
    def load_session(flags):
        logging.info('Loading session.')
        session_file = Path(flags.load_session).resolve()
        session_file_dir = session_file.parent

        with session_file.open('rb') as fp:
            try:
                vars_ = pickle.load(fp)
                if vars_['file_type'] != HEADER_LABEL_SESSION:
                    err_msg = '\nTarget file is not a session file. It is \"{}\".'.format(vars_['file_type'])
                    logging.error(err_msg)
                    raise Exception(err_msg)

                name = vars_['name']
                creation_time = vars_['creation_time']
                magic_label = vars_['magic_label']
                config = vars_['config']
                ms_exp_files = vars_['ms_exp_files']
                iden_files = vars_['iden_files']
                console_hist = vars_['console_hist']
                viewer_hist = vars_['viewer_hist']
            except Exception:
                err_msg = '\nIncorrect file type or corrupted file.'
                logging.error(err_msg)
                raise

        session = Session()
        session.name = name
        session.creation_time = creation_time
        session.magic_label = magic_label
        session.config = config
        session.ms_exp_files = ms_exp_files
        session.iden_files = iden_files
        session.console_hist = console_hist
        session.viewer_hist = viewer_hist
        session.flags = flags

        if flags.fork:
            logging.info('The forked session was named as \"{}\"'.format(flags.fork))
            if len(list(Path.cwd().glob(flags.fork + '.*'))) > 0:
                if flags.ignore_errors:
                    wrn_msg = 'Working directory contains files with names ' + \
                              'same as the session name \"{}\". '.format(session.name) + \
                              'Original files will be overwritten.'
                    logging.warning(wrn_msg)
                    for f in list(Path.cwd().glob(flags.fork + '.*')):
                        f.unlink()
                else:
                    err_msg = '\nWorking directory contains files with names same as the session name \"{}\".'.format(
                        session.name)
                    logging.error(err_msg)
                    raise FileExistsError(err_msg)

            fmts = (FILE_EXTENSION_IDEN_LUT, FILE_EXTENSION_RANKED_SPECTRA, FILE_EXTENSION_CLUSTERS,
                    FILE_EXTENSION_RAW_CLUSTERS, FILE_EXTENSION_INTERNAL_INDEX, FILE_EXTENSION_SESSION)
            for f in list(session_file_dir.glob(session.name + '.*')):
                if f.is_file and f.suffix in fmts:
                    Path(Path.cwd(), flags.fork + f.suffix).symlink_to(f)
            log_file = session_file_dir.joinpath(session.name + FILE_EXTENSION_LOG)
            if log_file.exists():
                with log_file.open(encoding='utf-8') as lf:
                    with Path.cwd().joinpath(flags.fork + FILE_EXTENSION_LOG).open('w', encoding='utf-8') as nf:
                        nf.write(lf.read())
            session.name = flags.fork
        else:
            if session_file_dir != Path.cwd():
                logging.info('Changes working directory to "{}"'.format(session_file_dir))
                os.chdir(str(session_file_dir))

        logging.debug('Embedded configuration:\n' + str(session.config))
        session.config.parse_config(flags)
        ms_exp_files, iden_files = _gather_paths(flags)
        if flags.rebuild_iden_lut and len(iden_files) != 0:
            session.iden_files = iden_files
        return session

    @staticmethod
    def new_session(flags):
        logging.info('Creating new session.')
        session = Session()
        if flags.name:
            session.name = flags.name
            logging.info('The new session was named as \"{}\"'.format(flags.name))
        else:
            logging.info('No session name provided, using creation time as session name: {}'.format(session.name))
        if len(list(Path.cwd().glob('{}.*'.format(session.name)))) > 0:
            if flags.ignore_errors:
                wrn_msg = 'Working directory contains files with names ' + \
                          'same as the session name \"{}\". '.format(session.name) + \
                          'Original files will be overwritten.'
                logging.warning(wrn_msg)
            else:
                err_msg = '\nWorking directory contains files with names same as the session name \"{}\".'.format(
                    session.name)
                logging.error(err_msg)
                raise FileExistsError(err_msg)
        session.flags = flags
        session.config = Configuration().parse_config(flags)
        session.ms_exp_files, session.iden_files = _gather_paths(flags)
        if len(session.ms_exp_files) == 0:
            err_msg = '\nNo MS experiement file specified.'
            logging.error(err_msg)
            raise Exception(err_msg)
        return session
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def get_session():
    if _session:
        return _session
    else:
        err_msg = '\nSession has not been setup, cannot get session.'
        logging.error(err_msg)
        raise ImportError(err_msg)
    return


def set_session(session):
    global _session
    _session = session
    return


def setup(flags):
    if flags.load_session or flags.print_session:
        logging.info('Session file specified: {}'.format(flags.load_session))
        try:
            if not flags.load_session.endswith(FILE_EXTENSION_SESSION): flags.load_session += FILE_EXTENSION_SESSION
            session_file = Path.cwd().joinpath(flags.load_session).resolve()
        except FileNotFoundError:
            err_msg = '\nCould not find the session file: {}'.format(flags.load_session)
            logging.error(err_msg)
            raise
        if session_file.is_dir():
            err_msg = '\nThe path of the session file is not pointing to valid file.'
            logging.error(err_msg)
            raise IsADirectoryError
        set_session(Session.load_session(flags))

        # logging.debug('Mounting files.')
        # if _session.config.ii_finished.value:
        #     _session.mount_internal_index()
        # if _session.config.rt_finished.value:
        #     if Path.cwd().joinpath(_session.name + FILE_EXTENSION_RANKED_SPECTRA).exists():
        #         _session.mount_internal_index()
        # if _session.config.cg_finished.value:
        #     _session.mount_clusters()
        #     _session.clusters.connect()
        # if _session.config.id_finished.value:
        #     if Path.cwd().joinpath(_session.name + FILE_EXTENSION_IDEN_LUT).exists():
        #         _session.mount_identification_lut()
        #         _session.iden_lut.connect()
        #     else:
        #         wrn_msg = 'Lost identification lookup table.'
        #         logging.warning(wrn_msg)
        #         _session.config.id_finished.value = False
        # logging.debug('Files mounted.')
    else:
        set_session(Session.new_session(flags))
    return


def _gather_paths(flags):
    if flags.ignore_errors:
        logging.info('\"ignore_errors\" was set to True, missing files will be ignored.')
    all_files = flags.files[:]

    if flags.file_list:
        logging.info('File list provided, reading file list: {}.'.format(flags.file_list))
        files_in_list = []
        with flags.file_list.open(encoding='utf-8') as fp:
            for l in fp.readlines():
                stripped = l.strip()
                if stripped != '':
                    files_in_list.append(stripped)
        all_files.extend(files_in_list)

    ms_exp_files = []
    iden_files = []
    for f in all_files:
        try:
            temp = Path(f).resolve()
            logging.debug('Resolved file path: {}'.format(f))
        except FileNotFoundError:
            if flags.ignore_errors:
                wrn_msg = 'Missing file: {}'.format(f)
                logging.warning(wrn_msg)
                continue
            else:
                err_msg = '\nFile could not be found.' +\
                          '\nFile path: {}'.format(f)
                logging.error(err_msg)
            raise

        patterns = {'.mzxml': ms_exp_files, '.mzml': ms_exp_files, '.pep.xml': iden_files}
        if temp.is_dir():
            logging.info('Searching directory \"{}\"'.format(temp))
            for ff in temp.iterdir():
                if ff.is_file():
                    try:
                        lower_case = ''.join(ff.suffixes).lower()
                        patterns[lower_case].append(ff)
                        logging.debug('Found file: {}'.format(ff))
                    except KeyError:
                        pass
        else:
            try:
                lower_case = ''.join(temp.suffixes).lower()
                patterns[lower_case].append(temp)
            except KeyError:
                if flags.ignore_errors:
                    wrn_msg = 'Unknown format: {}'.format(''.join(temp.suffixes))
                    logging.warning(wrn_msg)
                    continue
                else:
                    err_msg = '\nFile with unknown format was specified.' +\
                              '\nThis program only supports *.mzXML, *.mzML and *.pep.xml.' +\
                              '\nFile path: {}'.format(temp)
                    raise Exception(err_msg)

    logging.info('Number of MS experiment files: {}'.format(len(ms_exp_files)))
    logging.info('Number of identification files: {}'.format(len(iden_files)))
    return np.array(ms_exp_files), np.array(iden_files)
# ====END OF CODE====

