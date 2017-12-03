# -*- coding: utf-8 -*-
"""
Created on 08:35:04 20/11/2017

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
from distutils.core import setup
from Cython.Build import cythonize
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
setup(
    ext_modules=cythonize(["prcs/parallel/find_cluster.pyx",
                           "prcs/parallel/binning.pyx",
                           "prcs/parallel/cpu_kernel.pyx"])
)
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
