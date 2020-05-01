# -*- coding: utf-8 -*-
"""
Created on 12:31:37 08/11/2017

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
import numpy as np

from ClusterSheep._version import VERSION
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
NAME = 'ClusterSheep'

# size of the header of customized file format for metadata
FILE_HEADER_SIZE = 10240

DATA_TYPE_SIZE = {'uint8': 1,
                  'uint32': 4,
                  'uint64': 8,
                  'int32': 4,
                  'int64': 8,
                  'float32': 4,
                  'float64': 8}

NUMPY_DATA_TYPE = {'uint8': np.uint8,
                   'uint32': np.uint32,
                   'uint64': np.uint64,
                   'int32': np.int32,
                   'int64': np.int64,
                   'float32': np.float32,
                   'float64': np.float64}

# mzXML reading, only read the last portion of the file to find the offset to the index
II_INDEX_OFFSET_SEARCH_AREA = 1024
II_PRECURSOR_MASS_DATA_TYPE = 'float32'
II_PRECURSOR_CHARGE_DATA_TYPE = 'uint8'
II_FILE_ID_DATA_TYPE = 'uint32'
II_OFFSET_DATA_TYPE = 'uint64'
II_INTERNAL_ID_DATA_TYPE = 'uint64'
II_NATIVE_ID_DATA_TYPE = 'uint32'

ID_MODS_POS_DATA_TYPE = 'uint32'
ID_MODS_MASS_DATA_TYPE = 'float64'

RT_MZ_DATA_TYPE = 'int32'
RT_INTENSITY_DATA_TYPE = 'float32'
RT_INPUT_HIGH_PRECISION = False

CG_EDGE_DATA_TYPE = 'uint64'
CG_DOT_PRODUCT_DATA_TYPE = 'float32'
CG_OVERFLOWED_DATA_TYPE = 'bool'

CG_PRECURSOR_MASS_DATA_TYPE = II_PRECURSOR_MASS_DATA_TYPE
CG_MZ_DATA_TYPE = RT_MZ_DATA_TYPE
CG_INTENSITY_DATA_TYPE = RT_INTENSITY_DATA_TYPE
CG_PRECURSOR_TOLERANCE_DATA_TYPE = CG_PRECURSOR_MASS_DATA_TYPE
CG_BLOCK_DIMENSIONS_DATA_TYPE = 'uint32'
CG_OFFSET_DATA_TYPE = CG_EDGE_DATA_TYPE
CG_ALLOCATION_SIZE_DATA_TYPE = 'uint32'
CG_COUNTER_DATA_TYPE = 'uint32'

FG_PARENT_DATA_TYPE = CG_EDGE_DATA_TYPE
FG_RANK_DATA_TYPE = 'uint32'
FG_DATA_COPY_CHUNK_SIZE = 536870912  # in byte

GM_BUFFER_SIZE = 104857600  # in byte

C_DATA_TYPE = {'int32': 'int32_t',
               'uint32': 'uint32_t',
               'uint64': 'uint64_t',
               'float32': 'float',
               'bool': 'bool'}

FILE_EXTENSION_SESSION = '.cssess'
FILE_EXTENSION_LOG = '.cslogg'
FILE_EXTENSION_INTERNAL_INDEX = '.csindx'
FILE_EXTENSION_IDEN_LUT = '.csiden'
FILE_EXTENSION_RANKED_SPECTRA = '.csrksp'
FILE_EXTENSION_CLUSTERS = '.csclut'
FILE_EXTENSION_RAW_CLUSTERS = '.csrawc'
FILE_EXTENSION_TEMPORARY = '.cstemp'
FILE_EXTENSION_CONSOLE_HISTORY = '.cschis'
FILE_EXTENSION_VIEWER_HISTORY = '.csvhis'

HEADER_LABEL_SESSION = 'cs_session'
HEADER_LABEL_LOG = ''
HEADER_LABEL_INTERNAL_INDEX = 'cs_internal_index'
HEADER_LABEL_IDEN_LUT = 'cs_identification_lut'
HEADER_LABEL_RANKED_SPECTRA = 'cs_ranked_spectra'
HEADER_LABEL_CLUSTERS = 'cs_clusters'
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
