# -*- coding: utf-8 -*-
"""
Created on 01:30:32 03/11/2017

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
from share import code
import rlcompleter
import readline
from pathlib import Path
import sys

from envr.session import get_session
from property import *
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
def interactive_console(globals_):
    banner = 'You are now in Python interactive console, presses "Ctrl-D" or executes "exit()" to leave.'

    readline.parse_and_bind("tab: complete")
    readline.parse_and_bind('set show-all-if-ambiguous on')
    readline.parse_and_bind('"\C-r": reverse-search-history')
    readline.parse_and_bind('"\C-s": forward-search-history')
    history_file = Path.cwd().joinpath(session.name + FILE_EXTENSION_CONSOLE_HISTORY)
    with history_file.open('w', encoding='utf-8') as fp:
        fp.write(session.console_hist)
    readline.clear_history()
    readline.read_history_file(str(history_file))
    logging.info('Entering Python interactive console.')
    sys.ps1 = 'python >>> '
    sys.ps2 = 'python ... '
    console = code.InteractiveConsole(locals=globals_)

    pre_run = [
        'import misc',
        'misc.global_vars = globals()',
        'from misc import ipython',
        'import numpy as np',
        'from importlib import reload',
        'from prcs.cluster_viewing import cluster_viewer',
        'from prcs.cluster_viewing import get_graph',
        'iid = session.internal_index',
        'rks = session.ranked_spectra',
        'clu = session.clusters',
        'ide = session.iden_lut'
    ]
    for line in pre_run: console.push(line)
    console.interact(banner=banner)
    readline.write_history_file(str(history_file))
    with history_file.open(encoding='utf-8') as fp:
        session.console_hist = fp.read()
    history_file.unlink()
    readline.clear_history()
    logging.info('Leaving Python interactive console.')
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
