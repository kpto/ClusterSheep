# -*- coding: utf-8 -*-
"""
Created on 03:48:37 02/12/2017

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
from ClusterSheep.property import *
from ClusterSheep.prcs.template import create_template_file
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
HEADER = '''
{name} is a GPU accelerated program for MS2 spectral clustering.


Usage:

    {namel} [--switch-type-option --input-type-option=value] files


Examples:

    {namel} --file-list=/path/to/list --config=/path/to/config --name=run --ignore-errors
    {namel} --config=/path/to/config /path/to/files/
    {namel} --load-session=/path/to/a/session{sessext}
    {namel} --load-session=/path/to/a/session --print-session
    {namel} --list-gpus
    {namel} (shows this help and generates a configuration file template)


'''

OPTIONS = [
    ['--version', ('Print the version of this software.',)],
    ['--name=', ('Give the new session a name.',
                 'If not provided, the new session will be named by date and time.')],
    ['--fork=', ('Fork a session with a new name.',
                 'It is usually used to run clustering on prepared material for multiple times with different configuration.'),],
    ['--file-list=', ('Specify a file listing MS experiment files and identification files.',
                      'Lines starting with # will be skipped.')],
    ['--force', ('Basically try to have a complete run.',
                 'It ignores all non-critical errors and fallback to CPU if GPU failed.')],
    ['--use-cpu', ('Use CPU for clustering instead of GPU.',)],
    ['--ignore-errors', ('Ignores all non-critical errors.',)],
    ['--config=', ('Specify a file listing parameters and assigned values.',
                   'A template named "config_template" has already been created under the working directory.',
                   'Lines starting with # will be skipped.',
                   'Default parameters are used for a new session if no config file is specified.',
                   'If it is used with --load-session, the embedded config will be modified')],
    ['--dev-mode', ('Turns on developer mode.',)],
    ['--preparation-only', ('Stops after internal index building and rank transformation.',)],
    ['--stay-interactive', ('Go into cluster viewer after all processes.',)],
    ['--load-session=', ('Loads a session. File extension can be automatically filled.',)],
    ['--no-saving', ('Runs a full run and exports all clusters to a text file and deletes all other files.',)],
    ['--keep-trash', ('Keeps all intermediate files produced during the process.',)],
    ['--re-process', ('Abandons all files of a session.',)],
    ['--re-cluster', ('Abandons clustering results only.',)],
    ['--rebuild-iden-lut', ('Abandons the current identification lookup table and rebuild a new one.',
                            'User must specify identification files.')],
    ['--list-gpus', ('List all equipped CUDA device and ends program.',)],
    ['--print-session', ('Print the information of a session, must be used with --load-session.',)]
]
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def print_help():
    print(HEADER.format(name=NAME, namel=NAME.lower(), sessext=FILE_EXTENSION_SESSION))
    print('Options:\n')
    for o in OPTIONS:
        print('    {:25}'.format(o[0]), end='')
        for n, i in enumerate(o[1]):
            if n == 0: print(i)
            else: print(' '*29, i, sep='')
        print()
    print()

    try:
        create_template_file()
    except Exception as e:
        print('\033[33mFailed to create a configuration template file on the current working directory: {}\033[0m'
              .format(e.__class__.__name__))

    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
