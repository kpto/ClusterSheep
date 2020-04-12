# -*- coding: utf-8 -*-
"""
Created on 12:52:19 01/12/2017

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
import pycuda.driver as drv
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def list_gpus():
    drv.init()
    count = drv.Device.count()
    for id_ in range(count):
        d = drv.Device(id_)
        cc = d.compute_capability()
        print('Device id: {}    Device name: {}    Memory: {} MiB    Compute capability: {}'
              .format(id_, d.name(), int(d.total_memory() / (1024 ** 2)), '.'.join(map(str, cc))),
              end='')
        if cc[0] >= 2: print('    *valid')
        else: print()
    return
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
