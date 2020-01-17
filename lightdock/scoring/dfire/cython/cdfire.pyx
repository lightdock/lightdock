#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

import scipy.spatial
import numpy as np

cpdef calculate_dfire(receptor, receptor_coordinates, ligand, ligand_coordinates, 
                        dfire_dist_to_bins, dfire_energy, interface_cutoff=3.9):
    dist_matrix = scipy.spatial.distance.cdist(receptor_coordinates, ligand_coordinates)
    atom_indexes = np.where(dist_matrix <= 15.)
    dist_matrix *= 2.0
    dist_matrix -= 1.0
    energy = 0.
    cdef unsigned int d
    interface_receptor = []
    interface_ligand = []
    
    for i,j in zip(atom_indexes[0], atom_indexes[1]):
        rec_atom = receptor.objects[i]
        lig_atom = ligand.objects[j]
        rnuma = rec_atom.dfire_residue_index
        anuma = rec_atom.atom_index
        rnumb = lig_atom.dfire_residue_index
        anumb = lig_atom.atom_index
        # convert numpy.float64 to int
        d = dist_matrix[i][j]
        if d <= interface_cutoff:
            interface_receptor.append(i)
            interface_ligand.append(j)
        dfire_bin = dfire_dist_to_bins[d]-1
        energy += dfire_energy[rnuma][anuma][rnumb][anumb][dfire_bin]
    
    # Convert and change energy sign
    return (energy * 0.0157 - 4.7) * -1., set(interface_receptor), set(interface_ligand)
