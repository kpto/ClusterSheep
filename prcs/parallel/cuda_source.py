# -*- coding: utf-8 -*-
"""
Created on 12:49:33 17/11/2017

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


template = '''
#include <stdint.h>

__global__ void compute_dot_product(
        /* input */
        {pmass_dtype} precursor_mass[],
        {mz_dtype} mz[][{num_of_peaks}],
        {intensity_dtype} intensity[][{num_of_peaks}],
        {bkdim_dtype} block_dimensions[],
        {offset_dtype} offset[2],
        {alloc_size_dtype} *allocation_size,
        /* output */
        {counter_dtype} *counter,
        {edge_dtype} edge[][2],
        {dot_product_dtype} dot_product[],
        {overflowed_dtype} *overflowed
        ) {{
    const {bkdim_dtype} local_location_x = blockDim.x * blockIdx.x + threadIdx.x;
    const {bkdim_dtype} local_location_y = blockDim.y * blockIdx.y + threadIdx.y;

    const {pmass_dtype} *y_precursor_mass = precursor_mass;
    const {pmass_dtype} *x_precursor_mass = precursor_mass + block_dimensions[0];
    const {mz_dtype} (*y_mz)[{num_of_peaks}] = mz;
    const {mz_dtype} (*x_mz)[{num_of_peaks}] = mz + block_dimensions[0];
    const {intensity_dtype} (*y_intensity)[{num_of_peaks}] = intensity;
    const {intensity_dtype} (*x_intensity)[{num_of_peaks}] = intensity + block_dimensions[0];
    
    __shared__ {pmass_dtype} s_precursor_mass_y[{cuda_block_height}];
    __shared__ {pmass_dtype} s_precursor_mass_x[{cuda_block_width}];
    __shared__ {mz_dtype} s_mz_y[{cuda_block_height}][{num_of_peaks}];
    __shared__ {mz_dtype} s_mz_x[{cuda_block_width}][{num_of_peaks}];
    __shared__ {intensity_dtype} s_intensity_y[{cuda_block_height}][{num_of_peaks}];
    __shared__ {intensity_dtype} s_intensity_x[{cuda_block_width}][{num_of_peaks}];

    if(threadIdx.y == 0) {{
        s_precursor_mass_x[threadIdx.x] = x_precursor_mass[local_location_x];
        for(uint32_t i=0; i<{num_of_peaks}; ++i) {{
            s_mz_x[threadIdx.x][i] = x_mz[local_location_x][i];
            s_intensity_x[threadIdx.x][i] = x_intensity[local_location_x][i];
        }}
    }}

    if(threadIdx.x == 0) {{
        s_precursor_mass_y[threadIdx.y] = y_precursor_mass[local_location_y];
        for(uint32_t i=0; i<{num_of_peaks}; ++i) {{
            s_mz_y[threadIdx.y][i] = y_mz[local_location_y][i];
            s_intensity_y[threadIdx.y][i] = y_intensity[local_location_y][i];
        }}
    }}

    __syncthreads();

    if((local_location_y < block_dimensions[0] && local_location_x < block_dimensions[1]) &&
       (abs(s_precursor_mass_y[threadIdx.y]-s_precursor_mass_x[threadIdx.x]) <= {pmass_tol})) {{
        {dot_product_dtype} dp = 0;
        uint32_t y_ptr = 0;
        uint32_t x_ptr = 0;

        while(y_ptr < {num_of_peaks} && x_ptr < {num_of_peaks}) {{
            if(s_mz_y[threadIdx.y][y_ptr] == s_mz_x[threadIdx.x][x_ptr]) {{
                dp += s_intensity_y[threadIdx.y][y_ptr] * s_intensity_x[threadIdx.x][x_ptr];
                ++y_ptr;
                ++x_ptr;
            }} else if (s_mz_y[threadIdx.y][y_ptr] < s_mz_x[threadIdx.x][x_ptr]) {{
                ++y_ptr;
            }} else {{
                ++x_ptr;
            }}
        }}

        if(dp > {dp_tol}) {{
            {edge_dtype} global_location_y = local_location_y + offset[0];
            {edge_dtype} global_location_x = local_location_x + offset[1];
            if(global_location_x > global_location_y) {{
                uint32_t index = atomicAdd(counter, 1);
                if(index >= *allocation_size) {{
                    *overflowed = true;
                }}
                else {{
                    dot_product[index] = dp;
                    edge[index][0] = global_location_y;
                    edge[index][1] = global_location_x;
                }}
            }}
        }}
    }}
}}
'''
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def get_source_code():
    precursor_tolerance = session.config.cg_precursor_tolerance.value
    dot_product_threshold = session.config.cg_dot_product_threshold.value
    cuda_block_dimensions = session.config.gpu_cuda_block_dimensions.value
    cuda_block_height, cuda_block_width = cuda_block_dimensions
    num_of_peaks = session.config.rt_num_of_peaks.value
    info = {'pmass_dtype': C_DATA_TYPE[CG_PRECURSOR_MASS_DATA_TYPE],
            'mz_dtype': C_DATA_TYPE[CG_MZ_DATA_TYPE],
            'intensity_dtype': C_DATA_TYPE[CG_INTENSITY_DATA_TYPE],
            'ptol_dtype': C_DATA_TYPE[CG_PRECURSOR_TOLERANCE_DATA_TYPE],
            'bkdim_dtype': C_DATA_TYPE[CG_BLOCK_DIMENSIONS_DATA_TYPE],
            'offset_dtype': C_DATA_TYPE[CG_OFFSET_DATA_TYPE],
            'alloc_size_dtype': C_DATA_TYPE[CG_ALLOCATION_SIZE_DATA_TYPE],
            'counter_dtype': C_DATA_TYPE[CG_COUNTER_DATA_TYPE],
            'edge_dtype': C_DATA_TYPE[CG_EDGE_DATA_TYPE],
            'dot_product_dtype': C_DATA_TYPE[CG_DOT_PRODUCT_DATA_TYPE],
            'overflowed_dtype': C_DATA_TYPE[CG_OVERFLOWED_DATA_TYPE],
            'pmass_tol': precursor_tolerance,
            'dp_tol': dot_product_threshold,
            'cuda_block_height': cuda_block_height,
            'cuda_block_width': cuda_block_width,
            'num_of_peaks': num_of_peaks
    }
    return template.format(**info)
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
