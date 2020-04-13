# -*- coding: utf-8 -*-
"""
Created on 06:39:54 09/09/2017

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

from ClusterSheep.envr.session import get_session
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
global_vars = None
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
def ipython():
    from traitlets.config import Config
    import IPython

    c = Config()
    c.TerminalIPythonApp.display_banner = False
    c.InteractiveShell.confirm_exit = False

    logging.info('Entering IPython.')
    wrn_msg = 'Codes executed inside IPython cannot be recorded.'
    logging.warning(wrn_msg)
    IPython.start_ipython(argv=[], user_ns=global_vars, config=c)
    logging.info('Left IPython.')

    sys.ps1 = 'python >>> '
    sys.ps2 = 'python ... '
    sys.excepthook = sys.__excepthook__
    return


def bpython():
    import bpython
    bpython.embed(locals_=global_vars)
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
