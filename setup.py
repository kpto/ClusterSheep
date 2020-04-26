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
import os.path as path
import sys
from setuptools import setup, find_packages, Extension, Command

import numpy
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
here = os.path.abspath(os.path.dirname(__file__))

has_c_files = '.c' in [path.splitext(f)[1] for f in os.listdir(path.join(path.dirname(path.realpath(__file__)), 'lib'))]
USE_CYTHON = (not has_c_files) or os.getenv('USE_CYTHON') == 'true'

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [
    Extension('ClusterSheep.prcs.parallel.find_cluster', [os.path.join(here, 'lib/find_cluster' + ext)]),
    Extension('ClusterSheep.prcs.parallel.binning', [os.path.join(here, 'lib/binning' + ext)]),
    Extension('ClusterSheep.prcs.parallel.cpu_kernel', [os.path.join(here, 'lib/cpu_kernel' + ext)])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, force=True)
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def get_version(path):
    with open(path) as fp:
        for line in fp.readlines():
            if line.startswith('VERSION'):
                return line.split('=')[1].strip().strip("'")


setup(
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    ext_modules=extensions,
    include_dirs=[numpy.get_include()],
    entry_points={'console_scripts': ['clustersheep=ClusterSheep.main:main']},
    name='ClusterSheep',
    version=get_version(os.path.join(here, 'src/ClusterSheep/property.py')),
    author='Paul TO',
    author_email='kpto@connect.ust.hk',
    url='https://github.com/kpto/ClusterSheep',
    description='CUDA accelerated MS2 spectral clustering.',
    long_description=open('README.rst', encoding='utf-8').read(),
    long_description_content_type='text/x-rst',
    license='LGPL',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: C',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
