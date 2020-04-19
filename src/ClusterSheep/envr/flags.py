# -*- coding: utf-8 -*-
"""
Created on 19:21:02 25/10/2017

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
from pathlib import Path
from pprint import pformat
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
class Flags:

    def __init__(self):
        # Running options
        self.name = None
        self.fork = None
        self.files = None
        self.file_list = None
        self.force = False
        self.use_cpu = False
        self.ignore_errors = False
        self.config = None
        self.dev_mode = False
        self.checkpoint = False
        self.preparation_only = False
        self.stay_interactive = False
        self.load_session = None
        self.no_saving = False
        self.keep_trash = False
        self.re_process = False
        self.re_cluster = False
        self.rebuild_iden_lut = False
        # Non-running options
        self.list_gpus = False
        self.test_gpus = False
        self.print_session = False
        return

    def __repr__(self):
        return pformat(vars(self))
    
    def __str__(self):
        return '--name:' + str(self.name) + '\n' + \
               '--fork:' + str(self.fork) + '\n' + \
               '--file-list:' + str(self.file_list) + '\n' + \
               '--force:' + str(self.force) + '\n' + \
               '--use-cpu:' + str(self.use_cpu) + '\n' + \
               '--ignore-errors:' + str(self.ignore_errors) + '\n' + \
               '--config:' + str(self.config) + '\n' + \
               '--dev-mode:' + str(self.dev_mode) + '\n' + \
               '--checkpoint:' + str(self.checkpoint) + '\n' + \
               '--preparation-only:' + str(self.preparation_only) + '\n' + \
               '--stay-interactive:' + str(self.stay_interactive) + '\n' + \
               '--load-session:' + str(self.load_session) + '\n' + \
               '--no-saving:' + str(self.no_saving) + '\n' + \
               '--keep-trash:' + str(self.keep_trash) + '\n' + \
               '--re-process:' + str(self.re_process) + '\n' + \
               '--re-cluster:' + str(self.re_cluster) + '\n' + \
               '--rebuild-iden-lut:' + str(self.rebuild_iden_lut) + '\n' + \
               '--list-gpus:' + str(self.list_gpus) + '\n' + \
               '--test-gpus:' + str(self.test_gpus) + '\n' + \
               '--print-session:' + str(self.print_session)

    @staticmethod
    def _get_path(option, raw, type_):
        try:
            temp = Path(raw).resolve()
        except FileNotFoundError:
            err_msg = '\nPath \"{}\" cannot be resolved (does not exist).'.format(raw)
            logging.error(err_msg)
            raise
        if type_ == 'file':
            if temp.is_dir():
                err_msg = '\nPath \"{}\" is pointing to a directory instead of a file.'.format(raw) + \
                          '\nOption \"{}\" requires a file path.'.format(option)
                logging.error(err_msg)
                raise NotADirectoryError(err_msg)
        if type_ == 'dir':
            if temp.is_file():
                err_msg = '\nPath \"{}\" is pointing to a file instead of a directory.'.format(raw) + \
                          '\nOption \"{}\" requires a directory path.'.format(option)
                logging.error(err_msg)
                raise IsADirectoryError(err_msg)
        return temp

    @staticmethod
    def _redundant_value_error(option, value):
        if value != '':
            err_msg = '\nOption \"{}\" cannot be used with any value. Value \"{}\" was provided.'.format(option, value)
            logging.error(err_msg)
            raise ValueError(err_msg)
        return

    @staticmethod
    def _no_value_error(option, value):
        if value == '':
            err_msg = '\nOption \"{}\" requires the path of the relevant file.'.format(option) + \
                      '\nInput the path in format \"{}=/path/to/the/file\"'.format(option)
            logging.error(err_msg)
            raise ValueError(err_msg)
        return

    @staticmethod
    def _must_exist_error(must_exist, parter):
        heads = [arg.split('=')[0] for arg in sys.argv]
        for i in must_exist:
            if i not in heads:
                err_msg = '\nOption "{}" must be used with "{}".'.format(parter, i)
                logging.error(err_msg)
                raise Exception(err_msg)
        return

    @staticmethod
    def _conflict_error(cannot_exist, foe):
        heads = [arg.split('=')[0] for arg in sys.argv]
        for i in cannot_exist:
            if i in heads:
                err_msg = '\nOption "{}" cannot be used with "{}".'.format(foe, i)
                logging.error(err_msg)
                raise Exception(err_msg)
        return

    def func_name(self, value):
        Flags._no_value_error('--name', value)
        self.name = value
        return

    def func_fork(self, value):
        option_name = '--fork'
        Flags._no_value_error(option_name, value)
        Flags._conflict_error(('--name',), option_name)
        Flags._must_exist_error(('--load-session',), option_name)
        self.fork = value
        return

    def func_file_list(self, value):
        Flags._no_value_error('--file-list', value)
        self.file_list = Flags._get_path('--file-list', value, 'file')
        return

    def func_force(self, value):
        Flags._redundant_value_error('--force', value)
        self.force = True
        self.ignore_errors = True
        return

    def func_use_cpu(self, value):
        option_name = '--use-cpu'
        Flags._redundant_value_error(option_name, value)
        Flags._conflict_error(('--preparation-only',), option_name)
        self.use_cpu = True
        return

    def func_ignore_errors(self, value):
        Flags._redundant_value_error('--ignore-errors', value)
        self.ignore_errors = True
        return

    def func_config(self, value):
        Flags._no_value_error('--config', value)
        self.config = Flags._get_path('--config', value, 'file')
        return

    def func_dev_mode(self, value):
        Flags._redundant_value_error('--dev-mode', value)
        self.dev_mode = True
        return

    def func_checkpoint(self, value):
        Flags._redundant_value_error('--checkpoint', value)
        self.checkpoint = True
        return

    def func_preparation_only(self, value):
        Flags._redundant_value_error('--preparation-only', value)
        self.preparation_only = True
        return

    def func_stay_interactive(self, value):
        Flags._redundant_value_error('--stay-interactive', value)
        self.stay_interactive = True
        return

    def func_load_session(self, value):
        option_name = '--load-session'
        Flags._no_value_error(option_name, value)
        Flags._conflict_error(('--name',), option_name)
        self.load_session = value
        return

    def func_no_saving(self, value):
        option_name = '--no-saving'
        Flags._redundant_value_error(option_name, value)
        Flags._conflict_error(('--fork', '--preparation-only', '--stay-interactive', '--load-session',
                               '--re-process', '--re-cluster', '--rebuild-iden_lut',
                               '--list-gpus', '--test-gpus', '--print-session'), option_name)
        self.no_saving = True
        return

    def func_keep_trash(self, value):
        Flags._redundant_value_error('--keep-trash', value)
        self.keep_trash = True
        return

    def func_re_process(self, value):
        option_name = '--re-process'
        Flags._redundant_value_error(option_name, value)
        Flags._must_exist_error(('--load-session',), option_name)
        self.re_process = True
        return

    def func_re_cluster(self, value):
        option_name = '--re-cluster'
        Flags._redundant_value_error(option_name, value)
        Flags._must_exist_error(('--load-session',), option_name)
        Flags._conflict_error(('--preparation-only',), option_name)
        self.re_cluster = True
        return

    def func_rebuild_iden_lut(self, value):
        option_name = '--rebuild-iden_lut'
        Flags._redundant_value_error(option_name, value)
        Flags._must_exist_error(('--load-session',), option_name)
        Flags._conflict_error(('--preparation-only',), option_name)
        self.rebuild_iden_lut = True
        return

    def func_list_gpus(self, value):
        option_name = '--list-gpus'
        Flags._redundant_value_error(option_name, value)
        Flags._conflict_error(('--name', '--fork', '--file-list', '--force', '--use-cpu', '--ignore-errors',
                               '--config', '-dev-mode', '--checkpoint', '--preparation-only', '--stay-interactive', '--load-session',
                               '--no-saving', '--keep-trash', '--re-process', '--re-cluster', '--rebuild-iden_lut',
                               '--test-gpus', '--print-session'), option_name)
        self.list_gpus = True
        return

    def func_test_gpus(self, value):
        option_name = '--test-gpus'
        Flags._redundant_value_error(option_name, value)
        Flags._conflict_error(('--name', '--fork', '--file-list', '--force', '--use-cpu', '--ignore-errors',
                               '--config', '-dev-mode', '--checkpoint', '--preparation-only', '--stay-interactive', '--load-session',
                               '--no-saving', '--keep-trash', '--re-process', '--re-cluster', '--rebuild-iden_lut',
                               '--list-gpus', '--print-session'), option_name)
        self.test_gpus = True
        return

    def func_print_session(self, value):
        option_name = '--print-session'
        Flags._redundant_value_error(option_name, value)
        Flags._must_exist_error(('--load-session',), option_name)
        Flags._conflict_error(('--name', '--fork', '--file-list', '--force', '--use-cpu', '--ignore-errors',
                               '--config', '-dev-mode', '--checkpoint', '--preparation-only', '--stay-interactive', '--no-saving',
                               '--keep-trash', '--re-process', '--re-cluster', '--rebuild-iden_lut', '--list-gpus',
                               '--test-gpus'), option_name)
        self.print_session = True
        return

    def get_function_map(self):
        function_map = {
            '--name': self.func_name,
            '--fork': self.func_fork,
            '--file-list': self.func_file_list,
            '--force': self.func_force,
            '--use-cpu': self.func_use_cpu,
            '--ignore-errors': self.func_ignore_errors,
            '--config': self.func_config,
            '--dev-mode': self.func_dev_mode,
            '--checkpoint': self.func_checkpoint,
            '--preparation-only': self.func_preparation_only,
            '--stay-interactive': self.func_stay_interactive,
            '--load-session': self.func_load_session,
            '--no-saving': self.func_no_saving,
            '--keep-trash': self.func_keep_trash,
            '--re-process': self.func_re_process,
            '--re-cluster': self.func_re_cluster,
            '--rebuild-iden-lut': self.func_rebuild_iden_lut,
            '--list-gpus': self.func_list_gpus,
            '--test-gpus': self.func_test_gpus,
            '--print-session': self.func_print_session
        }
        return function_map

    def parse_argv(self):
        logging.debug('Parsing arguments.')
        argv = sys.argv[1:]
        function_map = self.get_function_map()
        files = []
        for arg in argv:
            if arg.startswith('--'):
                try:
                    fragments = arg.split('=')
                    option = fragments[0]
                    value = ''.join(fragments[1:])
                    function_map[option](value)
                    logging.debug('Received argument: {}'.format(arg))
                except KeyError:
                    err_msg = '\nUnknown option \"{}\"'.format(fragments)
                    logging.error(err_msg)
                    raise
            else:
                files.append(arg)
        self.files = files
        return self
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
