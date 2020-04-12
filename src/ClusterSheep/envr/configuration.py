# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:17:54 2017

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
from pathlib import Path
import os

from envr.parameter_base import *
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
class Configuration:

    def __init__(self):
        # general parameters
        self.gr_config_raw_text = ParameterBase('', str, immutable=True)
        self.gr_logging_level = ParameterBase(10, int, range_=(0, 50))
        self.gr_printing_level = ParameterBase(20, int, range_=(0, 50))
        self.gr_num_of_threads = ParameterBase(os.cpu_count(), int, range_=(0, INF))
        self.gr_integrity_check = ParameterBase(True, bool, immutable=True)
        # identification import parameters
        self.id_prob_threshold = ParameterBase(0.9, float, range_=(0.0, 1.0))
        self.id_include_decoy = ParameterBase(False, bool)
        self.id_num_of_threads = self.gr_num_of_threads
        self.id_finished = ParameterBase(False, bool)
        # internal index parameters
        self.ii_input_limit = ParameterBase(0, int, range_=(0, INF))
        self.ii_true_precursor_mass = ParameterBase(False, bool)
        self.ii_min_num_peaks = ParameterBase(10, int, range_=(1, INF))
        self.ii_num_of_threads = self.gr_num_of_threads
        self.ii_finished = ParameterBase(False, bool)
        # rank transformation parameters
        self.rt_num_of_peaks = ParameterBase(50, int, range_=(1, 191))
        self.rt_mz_range = ParameterRange((0, 0), int, element_range=(0, INF))
        self.rt_precursor_removal_range = ParameterRange((-18.0, 6.0), float)
        self.rt_bins_per_th = ParameterBase(1, int, range_=(1, INF))
        self.rt_num_of_threads = self.gr_num_of_threads
        self.rt_finished = ParameterBase(False, bool)
        # clustering general parameters
        self.cg_precursor_tolerance = ParameterBase(1.1, float, range_=(0.0, INF))
        self.cg_dot_product_threshold = ParameterBase(0.7, float, range_=(0.0, 1.0))
        self.cg_block_dimensions = ParameterTuple((2048, 2048), int, element_range=(1, INF), size=2)
        self.cg_allocation_size_initial_divisor = ParameterBase(30, int, range_=(1, INF))
        self.cg_finished = ParameterBase(False, bool)
        # clustering gpu parameters
        self.gpu_use_gpu = ParameterBase(True, bool)
        self.gpu_gpus_used = ParameterList([], int, element_range=(0, INF))
        self.gpu_cuda_block_dimensions = ParameterTuple((32, 32), int, element_range=(1, INF), size=2)
        self.gpu_use_fmad = ParameterBase(True, bool)
        self.gpu_processes_per_device = ParameterBase(1, int, range_=(1, INF))
        self.gpu_threads_per_device = ParameterBase(2, int, range_=(1, INF))
        # clustering cpu parameters
        self.cpu_num_of_threads = self.gr_num_of_threads
        # graph making parameters
        self.gm_num_of_threads = self.gr_num_of_threads
        # cluster refinement parameters
        self.cr_outlier_threshold = ParameterList([], int, element_range=(0, INF))
        self.cr_keep_raw = ParameterBase(False, bool)
        self.cr_num_of_threads = self.gr_num_of_threads
        self.cr_finished = ParameterBase(False, bool)
        return
    
    def __str__(self):
        return 'gr_logging_level:' + str(self.gr_logging_level) + '\n' +\
               'gr_printing_level:' + str(self.gr_printing_level) + '\n' + \
               'gr_num_of_threads:' + str(self.gr_num_of_threads) + '\n' + \
               'id_prob_threshold:' + str(self.id_prob_threshold) + '\n' +\
               'id_include_decoy:' + str(self.id_include_decoy) + '\n' + \
               'id_num_of_threads:' + str(self.id_num_of_threads) + '\n' + \
               'id_finished:' + str(self.ii_finished) + '\n' + \
               'ii_input_limit:' + str(self.ii_input_limit) + '\n' +\
               'ii_true_precursor_mass:' + str(self.ii_true_precursor_mass) + '\n' +\
               'ii_min_num_peaks:' + str(self.ii_min_num_peaks) + '\n' + \
               'ii_num_of_threads:' + str(self.ii_num_of_threads) + '\n' + \
               'ii_finished:' + str(self.ii_finished) + '\n' +\
               'rt_num_of_peaks:' + str(self.rt_num_of_peaks) + '\n' + \
               'rt_mz_range:' + str(self.rt_mz_range) + '\n' + \
               'rt_precursor_removal_range:' + str(self.rt_precursor_removal_range) + '\n' +\
               'rt_bins_per_th:' + str(self.rt_bins_per_th) + '\n' + \
               'rt_num_of_threads:' + str(self.rt_num_of_threads) + '\n' + \
               'rt_finished:' + str(self.rt_finished) + '\n' +\
               'cg_precursor_tolerance:' + str(self.cg_precursor_tolerance) + '\n' +\
               'cg_dot_product_threshold:' + str(self.cg_dot_product_threshold) + '\n' +\
               'cg_block_dimensions:' + str(self.cg_block_dimensions) + '\n' +\
               'cg_allocation_size_initial_divisor' + str(self.cg_allocation_size_initial_divisor) + '\n' +\
               'cg_finished:' + str(self.cg_finished) + '\n' +\
               'gpu_use_gpu:' + str(self.gpu_use_gpu) + '\n' +\
               'gpu_gpus_used:' + str(self.gpu_gpus_used) + '\n' + \
               'gpu_cuda_block_dimensions:' + str(self.gpu_cuda_block_dimensions) + '\n' + \
               'gpu_use_fmad:' + str(self.gpu_use_fmad) + '\n' + \
               'gpu_processes_per_device:' + str(self.gpu_processes_per_device) + '\n' + \
               'gpu_threads_per_device:' + str(self.gpu_threads_per_device) + '\n' + \
               'cpu_num_of_threads:' + str(self.cpu_num_of_threads) + '\n' + \
               'gm_num_of_threads:' + str(self.gm_num_of_threads) + '\n' + \
               'cr_outlier_threshold:' + str(self.cr_outlier_threshold) + '\n' +\
               'cr_keep_raw:' + str(self.cr_keep_raw) + '\n' + \
               'cr_num_of_threads:' + str(self.cr_num_of_threads) + '\n' + \
               'cr_finished:' + str(self.cr_finished)

    def func_gr_verbose(self, value):
        temp_bool = ParameterBase(True, bool)
        temp_bool.value = value
        if temp_bool.value is True:
            self.gr_logging_level.value = logging.DEBUG
            self.gr_printing_level.value = logging.INFO
        else:
            self.gr_logging_level.value = logging.INFO
            self.gr_printing_level.value = logging.INFO
        return

    def func_gr_num_of_threads(self, value):
        temp_int = ParameterBase(0, int)
        temp_int.value = value
        num_logical = os.cpu_count()
        if temp_int.value <= 0 or temp_int.value > num_logical:
            self.gr_num_of_threads.value = num_logical
        else:
            self.gr_num_of_threads.value = temp_int.value
        return

    def func_id_prob_threshold(self, value):
        self.id_prob_threshold.value = value
        return

    def func_id_include_decoy(self, value):
        self.id_include_decoy.value = value
        return

    def func_id_num_of_threads(self, value):
        self.id_num_of_threads = self._resolve_num_of_threads(value)
        return

    def func_ii_input_limit(self, value):
        self.ii_input_limit.value = value
        return

    def func_ii_true_precursor_mass(self, value):
        self.ii_true_precursor_mass.value = value
        return

    def func_ii_min_num_peaks(self, value):
        self.ii_min_num_peaks.value = value
        return

    def func_ii_num_of_threads(self, value):
        self.ii_num_of_threads = self._resolve_num_of_threads(value)
        return

    def func_rt_num_of_peaks(self, value):
        self.rt_num_of_peaks.value = value
        return

    def func_rt_mz_range(self, value):
        self.rt_mz_range.value = value
        return

    def func_rt_precursor_removal_range(self, value):
        self.rt_precursor_removal_range.value = value
        return

    def func_rt_bins_per_th(self, value):
        self.rt_bins_per_th.value = value
        return

    def func_rt_num_of_threads(self, value):
        self.rt_num_of_threads = self._resolve_num_of_threads(value)
        return

    def func_cg_precursor_tolerance(self, value):
        self.cg_precursor_tolerance.value = value
        return

    def func_cg_dot_product_threshold(self, value):
        self.cg_dot_product_threshold.value = value
        return

    def func_cg_block_dimensions(self, value):
        self.cg_block_dimensions.value = value
        self._check_block_dimension_validity()
        return

    def func_cg_allocation_size_initial_divisor(self, value):
        self.cg_allocation_size_initial_divisor.value = value
        return

    def func_gpu_use_gpu(self, value):
        self.gpu_use_gpu.value = value
        return

    def func_gpu_gpus_used(self, value):
        self.gpu_gpus_used.value = value
        return

    def func_gpu_cuda_block_dimensions(self, value):
        self.gpu_cuda_block_dimensions.value = value
        self._check_block_dimension_validity()
        return

    def func_gpu_use_fmad(self, value):
        self.gpu_use_fmad.value = value
        return

    def func_gpu_processes_per_device(self, value):
        self.gpu_processes_per_device.value = value
        return

    def func_gpu_threads_per_device(self, value):
        self.gpu_threads_per_device.value = value
        return

    def func_cpu_num_of_threads(self, value):
        self.cpu_num_of_threads.value = self._resolve_num_of_threads(value)
        return

    def func_gm_num_of_threads(self, value):
        self.gm_num_of_threads.value = self._resolve_num_of_threads(value)
        return

    def func_cr_outlier_threshold(self, value):
        self.cr_outlier_threshold.value = value
        return

    def func_cr_keep_raw(self, value):
        self.cr_keep_raw.value = value
        return

    def func_cr_num_of_threads(self, value):
        self.cr_num_of_threads.value = self._resolve_num_of_threads(value)
        return

    def _check_block_dimension_validity(self):
        cuda_bk_h, cuda_bk_w = self.gpu_cuda_block_dimensions.value[0], self.gpu_cuda_block_dimensions.value[1]
        if cuda_bk_h * cuda_bk_w > 1024:
            err_msg = ('\nTotal number of threads inside a CUDA block cannot exceed 1024.'
                       '\nCurrent: {} x {} = {}'.format(cuda_bk_h, cuda_bk_h, cuda_bk_h*cuda_bk_w))
            logging.error(err_msg)
            raise ValueError(err_msg)
        bk_h, bk_w = self.cg_block_dimensions.value[0], self.cg_block_dimensions.value[1]
        if bk_h % cuda_bk_h != 0 or bk_w % cuda_bk_w != 0:
            err_msg = ('\nBlock height/width must be multiple of CUDA block height/width'
                       '\nCurrent: Block height: {}    CUDA block height: {}'
                       '\n         Block width: {}    CUDA block width: {}'.
                       format(bk_h, cuda_bk_h, bk_w, cuda_bk_w))
            logging.error(err_msg)
            raise ValueError(err_msg)
        if bk_h * bk_w > 2**32:
            err_msg = ('\nTotal number of cells inside a block cannot exceed {}.'
                       '\nCurrent: {} x {} = {}'.format(bk_h, bk_w, bk_h*bk_w))
            logging.error(err_msg)
            raise ValueError(err_msg)
        return

    def get_function_map(self):
        function_map = {
            'gr_verbose': self.func_gr_verbose,
            'gr_num_of_threads': self.func_gr_num_of_threads,
            'id_prob_threshold': self.func_id_prob_threshold,
            'id_include_decoy': self.func_id_include_decoy,
            'id_num_of_threads': self.func_id_num_of_threads,
            'ii_input_limit': self.func_ii_input_limit,
            'ii_true_precursor_mass': self.func_ii_true_precursor_mass,
            'ii_min_num_peaks': self.func_ii_min_num_peaks,
            'ii_num_of_threads': self.func_ii_num_of_threads,
            'rt_num_of_peaks': self.func_rt_num_of_peaks,
            'rt_mz_range': self.func_rt_mz_range,
            'rt_precursor_removal_range': self.func_rt_precursor_removal_range,
            'rt_bins_per_th': self.func_rt_bins_per_th,
            'rt_num_of_threads': self.func_rt_num_of_threads,
            'cg_precursor_tolerance': self.func_cg_precursor_tolerance,
            'cg_dot_product_threshold': self.func_cg_dot_product_threshold,
            'cg_block_dimensions': self.func_cg_block_dimensions,
            'cg_allocation_size_initial_divisor': self.func_cg_allocation_size_initial_divisor,
            'gpu_use_gpu': self.func_gpu_use_gpu,
            'gpu_gpus_used': self.func_gpu_gpus_used,
            'gpu_cuda_block_dimensions': self.func_gpu_cuda_block_dimensions,
            'gpu_use_fmad': self.func_gpu_use_fmad,
            'gpu_processes_per_device': self.func_gpu_processes_per_device,
            'gpu_threads_per_device': self.func_gpu_threads_per_device,
            'cpu_num_of_threads':self.func_cpu_num_of_threads,
            'gm_num_of_threads': self.func_gm_num_of_threads,
            'cr_outlier_threshold': self.func_cr_outlier_threshold,
            'cr_keep_raw': self.func_cr_keep_raw,
            'cr_num_of_threads': self.func_cr_num_of_threads
        }
        return function_map

    def parse_config(self, flags):
        if flags.config:
            config_file = Path.cwd().joinpath(flags.config).resolve()
            self._read_file(config_file)
        else:
            logging.info('No config file specified, configuration remains unmodified (or default).')
        return self

    def _resolve_num_of_threads(self, value):
        temp_int = ParameterBase(0, int, range_=(-1, INF))
        temp_int.value = value
        num_logical = os.cpu_count()
        if temp_int.value == -1 or temp_int.value > num_logical:
            temp_int.value = num_logical
        elif temp_int.value == 0:
            return self.gr_num_of_threads.value
        else:
            return temp_int.value

    def _read_file(self, file):
        function_map = self.get_function_map()
        logging.debug('Parsing configuration file: {}'.format(file))
        with file.open(encoding='utf-8') as fp:
            text = fp.read()
            buffer = self.gr_config_raw_text.value
            lines = text.split('\n')
        for n, l in enumerate(lines):
            l = l.strip()
            if l != '' and not l.startswith('#'):
                try:
                    parameter, value = map(str.strip, l.split('='))
                    function_map[parameter](value)
                    buffer += (l + '\n')
                    logging.debug('Parsed line {}: {}'.format(n + 1, l))
                except Exception:
                    err_msg = '\nError occurred while parsing configuration file.' + \
                              '\nFile: {}'.format(file) + \
                              '\nLine {}: {}\n'.format(n + 1, l)
                    logging.error(err_msg)
                    raise
        buffer += '\n\n'
        self.gr_config_raw_text = ParameterBase(buffer, str, immutable=True)
        return
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
# ====END OF CODE====
