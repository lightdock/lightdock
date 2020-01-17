#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

from libc.math cimport sqrt, pow
import numpy as np
import scipy.spatial

cdef:
    unsigned int num_bins = 3
    unsigned int num_bins_spline = 5
    unsigned int num_atom_types = 32
    int r = -1
    float max_distance = 8.
    float min_distance = 2.

    r1Spline = [2.0, 3.0, 4.0, 4.5, 6.0]
    r2Spline = [3.0, 4.0, 4.5, 6.0, 8.0]

    r1 = [2.0, 3.5, 5.0]
    r2 = [3.5, 5.0, 8.0]

    float dist = 0.
    unsigned int itype = 0
    unsigned int jtype = 0

cpdef get_distance_to_bin(float dist):
    if dist >= max_distance or dist < min_distance:
        return -1
    i = num_bins_spline - 1
    while dist < r1Spline[i] or dist >= r2Spline[i]:
        i -= 1
    return i

cpdef calculate_pisa(receptor, receptor_coordinates, ligand, ligand_coordinates, pisa_energy, interface_cutoff):
    num_contacts = np.zeros((num_atom_types, num_atom_types, num_bins), dtype=np.float)
    energy = 0.
    dist_matrix = scipy.spatial.distance.cdist(receptor_coordinates, ligand_coordinates)
    atom_indexes = np.where((dist_matrix >= min_distance) & (dist_matrix <= max_distance))
    interface_receptor = []
    interface_ligand = []

    for i,j in zip(atom_indexes[0], atom_indexes[1]):
        dist = dist_matrix[i][j]
        if dist <= interface_cutoff:
            interface_receptor.append(i)
            interface_ligand.append(j)
        r = get_distance_to_bin(dist)
        itype = receptor.objects[i].pisa_type
        jtype = ligand.objects[j].pisa_type
        # swap types, so that itype is always <= jtype
        if jtype < itype:
            jtype, itype = itype, jtype

        if r == 0:
            num_contacts[itype-1][jtype-1][0] += 1.
        elif r == 1:
            num_contacts[itype-1][jtype-1][0] += (4. - dist)
            num_contacts[itype-1][jtype-1][1] += (dist - 3.)
        elif r == 2:
            num_contacts[itype-1][jtype-1][1] += 1.
        elif r == 3:
            num_contacts[itype-1][jtype-1][1] += (4. - (dist/3.))
            num_contacts[itype-1][jtype-1][2] += (dist/3.)-3.
        elif r == 4:
            num_contacts[itype-1][jtype-1][2] += 1.
        else:
            pass

    for i in range(num_atom_types):
        for j in range(num_atom_types):
            for r in range(num_bins):
                energy += num_contacts[i][j][r] * pisa_energy[i][j][r]

    return energy * -1., set(interface_receptor), set(interface_ligand)
