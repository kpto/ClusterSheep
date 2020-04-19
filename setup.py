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
import os
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
here = os.path.abspath(os.path.dirname(__file__));

extensions = [Extension('ClusterSheep.prcs.parallel.find_cluster', [os.path.join(here, 'src/ClusterSheep/prcs/parallel/find_cluster.pyx')]),
              Extension('ClusterSheep.prcs.parallel.binning', [os.path.join(here, 'src/ClusterSheep/prcs/parallel/find_cluster.pyx')]),
              Extension('ClusterSheep.prcs.parallel.cpu_kernel', [os.path.join(here, 'src/ClusterSheep/prcs/parallel/find_cluster.pyx')])]
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def get_version(path):
    with open(path) as fp:
        for line in fp.readlines():
            if line.startswith('VERSION'):
                return line.split(line.split('=')[1].strip().strip("'"))[1]


setup(
    name='ClusterSheep',
    version=get_version(os.path.join(here, 'src/ClusterSheep/property.py')),
    author='Paul TO',
    author_email='kpto@connect.ust.hk',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    ext_modules=cythonize(extensions),
    entry_points={'console_scripts': ['clustersheep=ClusterSheep.main:main']}
)
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
