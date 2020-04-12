# -*- coding: utf-8 -*-
"""
Created on 15:37:51 23/11/2017

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
cimport numpy as np
import numpy as np
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
ctypedef np.int32_t RT_MZ_DATA_TYPE
ctypedef np.float32_t RT_INTENSITY_DATA_TYPE
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def _binning(np.float32_t[:] mz, np.float32_t[:] intensity,
             RT_MZ_DATA_TYPE[:] new_mz, RT_INTENSITY_DATA_TYPE[:] new_intensity):

    cdef RT_MZ_DATA_TYPE previous = -1
    cdef RT_MZ_DATA_TYPE mz_int
    cdef np.uint32_t cursor = -1
    cdef np.uint32_t new_length
    cdef np.uint32_t length = mz.shape[0]
    cdef np.uint32_t i

    for i in range(length):
        mz_int = int(mz[i])
        if mz_int == previous:
            new_intensity[cursor] += intensity[i]
        else:
            cursor += 1
            new_mz[cursor] = previous = mz_int
            new_intensity[cursor] = intensity[i]
    new_length = cursor + 1
    return new_length
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
