# -*- coding: utf-8 -*-
"""
Created on 15:34:27 16/09/2017

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
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
CONFIG_FILE_TEMPLATE = '''### FRAME BEGIN
READ ME BEFORE ANY MODIFICATION:
You can find all settings available for users in here.

All parameter lines follow the format of "ParameterName = DefaultValue".
Each parameter line is followed by lines of description.
The first line states its data type and the allowed value range.
The resting lines explain the parameter.

Uncomment a parameter line to activate it.
# gr_ignore_files_error = False    <- commented line
gr_ignore_files_error = False      <- uncommented line
Input values must follow Python syntax and have correct data type.
Incorrect data types will not be converted implicitly to prevent ambiguity.

In case that you are not a Python user, here is a crash course.
        boolean(bool):  True / False (case sensitive)
         integer(int):  -1 / 0 / 1 / 2 / 3
floating-point(float):  -1.0 / 0.0 / 1.0 / 1.5
          string(str):  'A string' / 'Another string' / 'A valid "string" <- double quotation mark is okay'
                        'Not a valid 'string' <- no single quotation mark within a string'
                tuple:  (-1,) / (-1, 0) / (-1, 0, 1) / (-1.0, 0.0, 1.0) / ('string', 'another string')
                           ^---DO NOT MISS the comma even if there is only one element in a tuple.
                 list:  [-1] / [-1, 0] / [-1, 0, 1] / [-1.0, 0.0, 1.0] / ['string', 'another string']
If you have difficulty of understanding above rules, imitate the style of the default value.

The prefix of a parameter hints what it is related to.
"gr"  > "general"            : It can be referenced by any module.
"id"  > "identification"     : Identification import settings such as probability threshold.
"ii"  > "internal index"     : Index building settings such as spectra input and filtering.
"rt"  > "rank transformation": Rank transformation settings such as spectrum touching.
"cg"  > "clustering general" : Clustering settings that is applicable to both CPU and GPU.
"gpu" > "clustering by GPU"  : GPU specific clustering setting.
"cpu" > "clustering by CPU"  : CPU specific clustering setting.
"gm"  > "graph making"       : Conversion of edge list to graph objects.
"cr"  > "cluster refinement" : Cluster refinement such as refinement threshold.

Parameters are arranged by their importance, usually you only need to modify parameters under "USUAL" section.
Though this program prohibits impossible values, it DOES NOT handle extreme values.
It does not guarantee the stability under extreme settings.

* If you changed the configuration after a program run, you NEED TO re-run the related processes.
  For example, changing "rt"-starting parameters requires re-run of clustering and processes following it.
### FRAME END


USUAL parameters:

gr_num_of_threads = 0
/(int, 0 - INF)
/Global value of number of threads spawned in processes supporting parallel computing.
/Python does not have true multi-threading, this is achieved by forking Python process.
/Though this parameter says "thread", here it actually means "process".
/The same applies on all thread number related parameters, EXCEPT parameters starting with "gpu".
/If you want to tune this parameter of each process separately, go to TECHNICAL section.
/Value of 0 means all logical units.

id_prob_threshold = 0.9
/(float, 0.0 - 1.0)
/Only identifications with iProphet/PeptideProphet probability higher than the set threshold are imported.

id_include_decoy = False
/(bool)
/Whether identifications hitting decoys should be imported.

ii_input_limit = 0
/(int, 0 - INF)
/Limits the number of input spectra, 0 means unlimited.
/It is a file level limiter, the resting spectra within the same MS experiment file are still imported.

ii_true_precursor_mass = False
/(bool)
/Uses actual masses of precursors by multiplying m/z values to estimated charges instead of m/z values.

ii_min_num_peaks = 10
/(int, 1 - INF)
/Spectra with numbers of peaks smaller than the set value are skipped.

rt_mz_range = 0.0,2000.0
The range of the acceptable m/z. Peaks that are out of the range are discarded.
Expects a tuple with a structure of (float x,float y).

rt_precursor_removal_range = (-18.0, 6.0)
/(tuple of float, -INF - INF)
/Removes peaks near precursor m/z possibly caused by unfragmented precursors.
/Peaks within and equal to (precursor m/z + left value) - (precursor m/z + right value) are removed.
/For example, peaks with m/z values within and equal to 482.0 - 506.0 are removed if the precursor m/z is 500.0 .

cg_precursor_tolerance = 1.1
/(float, 0.0 - INF)
/A spectrum pair with precursor difference larger than the set threshold are assumed to be dissimilar and skipped.
/The higher of the value, the longer of the running time.

cg_dot_product_threshold = 0.7
/(float, 0.0 - 1.0)
/An edge is formed only if a spectrum pair have the similarity score larger than the set threshold.

cr_outlier_threshold = []
/In cluster refinement, a spectrum having a betweenness centrality value with m-score higher than the set threshold is identified as a bridge and removed.
/The list can have any length. If cluster refinement is not desired, input an empty list [] (default value).
/Refined clusters can be refined again. The length of the list is the number of multiple refinements.
/For example, [50, 20, 10] means three rounds of refinement, starting from threshold 50, then 20, and lastly 10.
/The longer of the list, the longer of the running time.


ADVANCED parameters:

gr_verbose = False
/(bool)
/Prints debugging level messages on terminal.
/You rarely need it as the log file has already logged them.

rt_num_of_peaks = 50
/(int, 1 - 191)
/The number of peaks retained after rank transformation.
/Since spectra needs to be stored in shared memory inside GPU, the size of the memory limited the number.
/For default CUDA block dimensions, the maximum number is 191 for GPUs equipped with 96KiB shared memory and 127 for 64KiB.
/Though a smaller CUDA block dimensions allows a higher maximum value, the range here is locked.
/If you really want to unlock the limit, go to envr/configuration.py and change the limit setting.
/Increasing the value causes a longer running time.

rt_bins_per_th = 1
/(int, 1 - INF)
/Number of bins within 1 Th.
/For example, value of 1 yields bins of [0.0-0.999..., 1.0-1.999..., 2.0-2.999..., and so on]
/Value of 2 yields bins of [0.0-0.499..., 0.5-0.999..., 1.0-1.499..., 1.5-1.999..., and so on]

gpu_use_gpu = True
/(bool)
/Uses GPU to compute the similarity matrix. CPU is used instead if it is False.

gpu_gpus_used = []
/(list of int, 0 - INF)
/Specify which gpus are going to be used. If you want to use all gpus found, input an empty list (default value).
/If you do not know the IDs of gpus, use option --list-gpus.
/Sometimes you may want to keep the GPU drawing your GUI relaxed to prevent instability of the desktop environment.

cr_keep_raw = False
/(bool)
/Keep the untouched clusters separately after refinement.

cg_use_justin_similarity_func = False
/(bool)
/When true, justin similarity function is used instead of simple dot product function.
/The function is
/   1 / (1 + 
/            EXP(
/                -(
/                    ({cg_dbscan_distance_func_pmass_multiplier} * {precursor mass difference}) +
/                    ({cg_dbscan_distance_func_dp_multiplier} * {dot product}) +
/                    {cg_dbscan_distance_func_constant}
/                 )
/               )
/       )

cg_justin_similarity_func_pmass_multiplier = -0.5275
/(float, -INF - INF)
/Used when justin similarity function is enabled. The multiplier of precursor mass in distance function.
/See description of cg_use_dbscan for the distance function.

cg_justin_similarity_func_dp_multiplier = 4.5572
/(float, -INF - INF)
/Used when justin similarity function is enabled. The multiplier of dot product in distance function.
/See description of cg_use_dbscan for the distance function.

cg_justin_similarity_func_constant = -1.8332
/(float, -INF - INF)
/Used when justin similarity function is enabled. The constant in distance function.
/See description of cg_use_dbscan for the distance function.

cg_use_dbscan = False
/(bool)
/When true, spectra will be clustered using DBSCAN, except that similarity is used instead of distance.

cg_dbscan_min_points = 3
/(int, 1 - INF)
/The minimum points parameter of DBSCAN, ignored if cg_use_dbscan is False.


TECHNICAL parameters:

id_num_of_threads = -1
/(int, -1 - INF)
/Number of threads spawned in identification import process.
/Value of -1 means following the global value, 0 means all logical units

ii_num_of_threads = -1
/(int, -1 - INF)
/Number of threads spawned in index building process.
/Value of -1 means following the global value, 0 means all logical units

rt_num_of_threads = -1
/(int, -1 - INF)
/Number of threads spawned in rank transformation process.
/Value of -1 means following the global value, 0 means all logical units

cg_block_dimensions = (2048, 2048)
/(tuple of int, 1 - INF)
/The dimensions of blocks diving the similarity matrix.
/A larger size reduces division overhead while a smaller size can have a better fitting to the necessary area.
/It also determines the number of spectra uploaded to the computing device.
/The LEFT value is the height and the RIGHT value is the width.

cg_allocation_size_initial_divisor = 30
/(int, 1 - INF)
/It determines the length of the edge list allocated in the device for the first try of computing a block.
/The length is the possible maximum number of edges produced in a block over the value of this parameter.
/For example, with the default value of block dimensions, the possible maximum number of edges produced is 2048 x 2048.
/This number divided by the value of this parameter becomes the initial length of the edge list allocated.
/A larger value reduces the time needed of data transfer between the host and the device but increases the chance of overflow.
/This program utilizes CUDA stream, as long as the transfer time is shorter than  the kernel execution time, it will not be exposed.
/If you see many overflow messages appears during clustering process in the log file, you may try reducing this value.

gpu_cuda_block_dimensions = (32, 32)
/(tuple of int, 1 - INF)
/The dimensions of CUDA blocks dividing a block.
/It determines the space of shared memory in GPU requested by a CUDA block and may possibly affect the performance.
/The LEFT value is the height and the RIGHT value is the width.

gpu_use_fmad = True
/(bool)
/Uses fused multiply–add. It reduces computation of a = a + (b * c) to one step.
/It may improves the performance and the precision.
/Dot products computed may be slightly different from dot products computed by a CPU.
/If you want to verify the computed result from GPU, you can turn it off to get the exact result as CPU.

gpu_processes_per_device = 1
/(int, 1 - INF)
/The number of processes for each CUDA device.
/Since different process cannot attach to the same CUDA context, concurrency achieved by CUDA stream is not possible between processes.
/If the CPU is too slow, Python may be not fast enough to prepare material before the GPU finished kernel execution.
/In this case, you can increase the number of processes to sacrifice the perfect concurrency and get the true "multi-threading".

gpu_threads_per_device = 2
/(int, 1 - INF)
/The number of threads (not Python processes) for each CUDA device.
/CUDA stream is utilized if multiple threads are acquiring the same device.
/This can hide the latency caused by data transfer between the host and the device.
/Be reminded that Python serializes all threads. Higher number means that a thread needs to prepare materials for multiple blocks.
/Since two CUDA streams can already achieve 100% GPU usage, you should not have reason to change this parameter.

cpu_num_of_threads = -1
/(int, -1 - INF)
/Number of threads spawned in clustering process in the case that CPU is used instead of GPUs.
/Value of -1 means following the global value, 0 means all logical units

gm_num_of_threads = -1
/(int, -1 - INF)
/Number of threads spawned in graph making process.
/Value of -1 means following the global value, 0 means all logical units

cr_num_of_threads
/(int, -1 - INF)
/Number of threads spawned in cluster refinement process.
/Value of -1 means following the global value, 0 means all logical units
'''
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def _draw_frame(lines):
    ranges = []
    counter = 0
    begin_found = False
    begin = 0
    while counter < len(lines):
        line = lines[counter]
        if begin_found:
            if line.startswith('### FRAME END'):
                end = counter
                ranges.append((begin, end))
                begin_found = False
        else:
            if line.startswith('### FRAME BEGIN'):
                begin = counter
                begin_found = True
        counter += 1
    for r in ranges:
        longest = max(map(len, lines[r[0]+1:r[1]]))
        lines[r[0]] = lines[r[1]] = '# ' + '+' + '-'*longest + '+'
        for n in range(r[0]+1, r[1]):
            lines[n] = '# ' + '|' + lines[n] + ' '*(longest-len(lines[n])) + '|'
    return lines


def _indentate(lines):
    for i in range(len(lines)):
        if lines[i].startswith('#') or lines[i] == '':
            continue
        elif lines[i].startswith('/'):
            lines[i] = '#   ' + lines[i][1:]
        else:
            lines[i] = '# ' + lines[i]
    return lines


def get_template():
    global CONFIG_FILE_TEMPLATE
    lines = _draw_frame(CONFIG_FILE_TEMPLATE.split('\n'))
    lines = _indentate(lines)
    return '\n'.join(lines)


def create_template_file():
    file = Path.cwd().joinpath('config_template')
    if file.exists():
        file.unlink()
    with file.open('w', encoding='utf-8') as fp:
        fp.write(get_template())
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
