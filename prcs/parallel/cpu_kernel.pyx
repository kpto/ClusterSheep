# -*- coding: utf-8 -*-
"""
Created on 16:04:53 23/11/2017

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
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
ctypedef np.float32_t CG_PRECURSOR_MASS_DATA_TYPE
ctypedef np.int32_t CG_MZ_DATA_TYPE
ctypedef np.float32_t CG_INTENSITY_DATA_TYPE
ctypedef np.uint32_t CG_BLOCK_DIMENSIONS_DATA_TYPE
ctypedef np.uint64_t CG_OFFSET_DATA_TYPE
ctypedef np.uint64_t CG_EDGE_DATA_TYPE
ctypedef np.float32_t CG_DOT_PRODUCT_DATA_TYPE
ctypedef np.uint32_t CG_COUNTER_DATA_TYPE
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def _cpu_kernel(CG_PRECURSOR_MASS_DATA_TYPE[:] precursor_mass,
                CG_MZ_DATA_TYPE[:,:] mz,
                CG_INTENSITY_DATA_TYPE[:,:] intensity,
                CG_BLOCK_DIMENSIONS_DATA_TYPE[:] block_dimensions,
                CG_OFFSET_DATA_TYPE[:] offset,
                CG_PRECURSOR_MASS_DATA_TYPE precursor_tolerance,
                CG_DOT_PRODUCT_DATA_TYPE dot_product_threshold,
                np.uint32_t num_of_peaks,
                CG_EDGE_DATA_TYPE[:,:] edge,
                CG_DOT_PRODUCT_DATA_TYPE[:] dot_product):

    cdef CG_COUNTER_DATA_TYPE count = 0
    cdef CG_BLOCK_DIMENSIONS_DATA_TYPE location_y, location_x,
    cdef CG_EDGE_DATA_TYPE global_location_y, global_location_x
    cdef CG_MZ_DATA_TYPE[:] temp_mz_y, temp_mz_x
    cdef CG_INTENSITY_DATA_TYPE[:] temp_intensity_y, temp_intensity_x
    cdef CG_DOT_PRODUCT_DATA_TYPE dp
    cdef CG_COUNTER_DATA_TYPE y_ptr, x_ptr

    for location_y in range(block_dimensions[0]):
        for location_x in range(block_dimensions[1]):
            if abs(precursor_mass[location_y] - precursor_mass[block_dimensions[0] + location_x]) <= precursor_tolerance:
                temp_mz_y = mz[location_y]
                temp_mz_x = mz[block_dimensions[0] + location_x]
                temp_intensity_y = intensity[location_y]
                temp_intensity_x = intensity[block_dimensions[0] + location_x]

                dp = 0.0
                y_ptr = 0
                x_ptr = 0

                while y_ptr < num_of_peaks and x_ptr < num_of_peaks:
                    if temp_mz_y[y_ptr] == temp_mz_x[x_ptr]:
                        dp += temp_intensity_y[y_ptr] * temp_intensity_x[x_ptr]
                        y_ptr += 1
                        x_ptr += 1
                    elif temp_mz_y[y_ptr] < temp_mz_x[x_ptr]:
                        y_ptr += 1
                    else:
                        x_ptr += 1

                if dp > dot_product_threshold:
                    global_location_y = location_y + offset[0]
                    global_location_x = location_x + offset[1]
                    if global_location_x > global_location_y:
                        dot_product[count] = dp
                        edge[count][0] = global_location_y
                        edge[count][1] = global_location_x
                        count += 1
    return count
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
