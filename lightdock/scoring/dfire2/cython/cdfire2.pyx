#cython: language_level=3

import scipy.spatial
cimport numpy as np
from libc.math cimport sqrt

from cython cimport boundscheck, wraparound, nonecheck

@boundscheck(False)
@wraparound(False)
@nonecheck(False)
cpdef calculate_dfire2(int[:] res_index,
                       int[:] atom_index,
                       np.ndarray[np.float64_t, ndim=2] coordinates,
                       np.ndarray[np.float64_t, ndim=3] potentials,
                       np.uint32_t n, np.float64_t interface_cutoff):
    cdef np.float64_t energy = 0.
    cdef size_t i, j
    cdef double d = 0.
    cdef unsigned int b = 0
    cdef np.ndarray[np.float64_t, ndim=2] dist_matrix = scipy.spatial.distance.cdist(coordinates, coordinates)*2
    interface_receptor = []
    interface_ligand = []

    for i in range(n):
        for j in range(i+1, n):
            if res_index[i] != res_index[j]:
                d = dist_matrix.item(i,j)
                if d <= interface_cutoff:
                  interface_receptor.append(i)
                  interface_ligand.append(j)
                b = int(d)
                if b < 30:
                    energy += potentials.item(atom_index[i], atom_index[j], b)
                    #print "%d %d %4.6f %d %f" % (res_index[i], res_index[j], d, b, potentials.item(atom_index[i], atom_index[j], b))
    return energy/100., set(interface_receptor), set(interface_ligand)
