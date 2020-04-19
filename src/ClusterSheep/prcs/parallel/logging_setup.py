# -*- coding: utf-8 -*-
"""
Created on 16:14:36 12/09/2017

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
try:
    session = get_session()
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def logging_setup():
    stream_handler = None
    file_handler = None
    for h in logging.getLogger().handlers:
        stream_handler = h if type(h) is logging.StreamHandler else stream_handler
        file_handler = h if type(h) is logging.FileHandler or logging.NullHandler else file_handler

    if session.flags.dev_mode:
        file_formatter = logging.Formatter(
            '%(levelname)-8s %(asctime)-18s    PID:%(process)-6d    %(filename)-20s %(funcName)-20s    %(message)s', '%d%b%Y %H:%M:%S')
        stream_formatter = logging.Formatter('%(levelname)-8s %(asctime)-18s    PID:%(process)-6d    %(message)s', '%d%b%Y %H:%M:%S')
    else:
        file_formatter = logging.Formatter('%(levelname)-8s %(asctime)-18s    PID:%(process)-6d    %(message)s', '%d%b%Y %H:%M:%S')
        stream_formatter = logging.Formatter('%(levelname)-8s %(asctime)-18s    PID:%(process)-6d    %(message)s', '%d%b%Y %H:%M:%S')

    file_handler.setFormatter(file_formatter)
    stream_handler.setFormatter(stream_formatter)
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
