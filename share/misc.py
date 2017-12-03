# -*- coding: utf-8 -*-
"""
Created on 06:48:26 09/09/2017

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
import colorsys
import math

from envr.session import get_session
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
def with_in_range(range_, value):
    return range_[0] <= value <= range_[1]


def extract_header(file):
    with file.open('rb') as fp:
        block = fp.read(1000)
        loc = block.find(b'%end%')
        if loc == -1:
            err_msg = 'Could not find ending tag.'
            logging.error(err_msg)
            raise Exception(err_msg)
        return block[:loc]


def generate_colors(num_of_colors):
    colors = []
    if num_of_colors != 0:
        num_of_colors = math.ceil(num_of_colors/3)*3
        step = 1 / num_of_colors
        for i in range(num_of_colors):
            colors.append(colorsys.hsv_to_rgb(i*step, 1.0, 1.0))
        for i in range(len(colors)):
            colors[i] = colors[i] + (1.0,)
    return colors


def _clean_temporary(temp_storage):
    if session.flags.keep_trash: return
    logging.debug('Cleaning temporary files')
    for f in list(temp_storage.iterdir()): f.unlink()
    temp_storage.rmdir()
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
