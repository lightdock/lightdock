#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

import scipy.spatial
import numpy as np

cpdef calculate_ddna(receptor, receptor_coordinates, ligand, ligand_coordinates, 
                     ddna_potentials, ddna_map, interface_cutoff=3.9):

    # Calculate all atom distances
    dist_matrix = scipy.spatial.distance.cdist(receptor_coordinates, ligand_coordinates)
    # Keep only atom indexes below 12A
    atom_indexes = np.where(dist_matrix <= 12.)

    energy = 0.
    U = 0.
    cdef unsigned int d
    interface_receptor = []
    interface_ligand = []
    
    for i,j in zip(atom_indexes[0], atom_indexes[1]):
        # convert numpy.float64 to int
        d = dist_matrix[i][j]
        # Calculate interface
        if d <= interface_cutoff:
            interface_receptor.append(i)
            interface_ligand.append(j)

        rec_atom = receptor.objects[i]
        lig_atom = ligand.objects[j]

        jj = ddna_map[int(d * 2.0)]
        if jj > 0 and jj <= 20:
            U = ddna_potentials[jj * 20*20 + rec_atom * 20 + lig_atom];
            if U < -5.0:
                U = 0.0
            energy += U
    
    # Convert and change energy sign
    return (energy * 0.0021297 - 5.4738) * -1., set(interface_receptor), set(interface_ligand)
