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
import rlcompleter
import readline
from pathlib import Path
import sys

from ClusterSheep.envr.session import get_session
from ClusterSheep.property import *
from ClusterSheep.share import code
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
class Checkpoint:

    def __init__(self, enabled, locals_):
        self.enabled = enabled
        self.locals = locals_
        return

    def __call__(self, next_process=None):
        if self.enabled:
            if next_process:
                print('Next process is "{}".'.format(next_process))
            print('Press Ctrl-D to proceed or input "exit()" to hard terminate the program.')
            interactive_console(self.locals, banner='', pre_run=[], true_exit=True)
        return
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def interactive_console(globals_, banner=None, pre_run=None, true_exit=False):
    if banner is None:
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

    if pre_run is None:
        pre_run = [
            'import ClusterSheep.misc',
            'import sys',
            'import os',
            'sys.path.append(os.getcwd())',
            'ClusterSheep.misc.global_vars = globals()',
            'from ClusterSheep.misc import ipython',
            'import numpy as np',
            'from importlib import reload',
            'from ClusterSheep.prcs.cluster_viewing import cluster_viewer',
            'from ClusterSheep.prcs.cluster_viewing import get_graph',
            'iid = session.internal_index',
            'rks = session.ranked_spectra',
            'clu = session.clusters',
            'ide = session.iden_lut'
        ]

    for line in pre_run: console.push(line)
    console.interact(banner=banner, true_exit=true_exit)
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
