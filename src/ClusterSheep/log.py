# -*- coding: utf-8 -*-
"""
Created on 20:40:22 04/11/2017

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
import sys
from datetime import datetime
from pathlib import Path

import ClusterSheep.envr.session
from ClusterSheep.property import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
class Log:

    def __init__(self):
        self.log_file_current = None
        self.file_handler = None
        self.stream_handler = None
        self.file_formatter = None
        self.stream_formatter = None
        return

    def default(self):
        logger = logging.getLogger()

        self.log_file_current= Path.cwd().joinpath(datetime.now().isoformat() + FILE_EXTENSION_LOG)
        self.file_handler = logging.FileHandler(str(self.log_file_current), mode='w')
        self.stream_handler = logging.StreamHandler(sys.stdout)

        if '--dev-mode' in sys.argv:
            self.file_formatter = logging.Formatter(
                '%(levelname)-8s %(asctime)-18s    %(filename)-25s %(funcName)-25s    %(message)s', '%d%b%Y %H:%M:%S')
            self.stream_formatter = logging.Formatter('%(levelname)-8s %(asctime)-18s    %(message)s', '%d%b%Y %H:%M:%S')
            self.file_handler.setLevel(logging.DEBUG)
            self.stream_handler.setLevel(logging.DEBUG)
        else:
            self.file_formatter = logging.Formatter('%(levelname)-8s %(asctime)-18s    %(message)s', '%d%b%Y %H:%M:%S')
            self.stream_formatter = logging.Formatter('%(levelname)-8s %(asctime)-18s    %(message)s', '%d%b%Y %H:%M:%S')
            self.file_handler.setLevel(logging.INFO)
            self.stream_handler.setLevel(logging.INFO)

        self.file_handler.setFormatter(self.file_formatter)
        self.stream_handler.setFormatter(self.stream_formatter)

        logger.addHandler(self.file_handler)
        logger.addHandler(self.stream_handler)

        logger.setLevel(0)

        logging.info('Program started with arguments: {}'.format(' '.join(sys.argv)))
        return

    def update(self):
        session = envr.session.get_session()

        logger = logging.getLogger()
        logger.removeHandler(self.file_handler)
        new_file_name = session.name + FILE_EXTENSION_LOG
        new_file = Path.cwd().joinpath(new_file_name)

        if not session.flags.no_saving:
            if not new_file.exists():
                new_file.touch()
            with new_file.open('a', encoding='utf-8') as nf:
                with self.log_file_current.open(encoding='utf-8') as lf:
                    nf.write(lf.read())
            self.file_handler = logging.FileHandler(new_file_name, mode='a')
            self.file_handler.setLevel(session.config.gr_logging_level.value)
            self.file_handler.setFormatter(self.file_formatter)
        else:
            self.file_handler = logging.NullHandler()

        self.log_file_current.unlink()
        self.log_file_current = new_file
        logger.addHandler(self.file_handler)

        if session.flags.dev_mode:
            self.file_handler.setLevel(logging.DEBUG)
            self.stream_handler.setLevel(logging.DEBUG)
        else:
            self.file_handler.setLevel(session.config.gr_logging_level.value)
            self.stream_handler.setLevel(session.config.gr_printing_level.value)
        logging.info('Logging level was set to File: {}, Stream: {}.'.format(self.file_handler.level,
                                                                             self.stream_handler.level))
        return
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
