"""Side-chain prediction"""

import scipy.spatial
import numpy as np
import Bio.PDB


chi_atoms = {"ALA": {},
             "ARG": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "CD"],
                     "x3": ["CB", "CG", "CD", "NE"], "x4": ["CG", "CD", "NE", "CZ"]},
             "ASN": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "OD1"]},
             "ASP": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "OD1"]},
             "CYS": {"x1": ["N", "CA", "CB", "SG"]},
             "GLU": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "CD"],
                     "x3": ["CB", "CG", "CD", "OE1"]},
             "GLN": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "CD"],
                     "x3": ["CB", "CG", "CD", "OE1"]},
             "GLY": {},
             "HIS": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "CD2"]},
             "ILE": {"x1": ["N", "CA", "CB", "CG1"], "x2": ["CA", "CB", "CG1", "CD1"]},
             "LEU": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "CD1"]},
             "LYS": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "CD"],
                     "x3": ["CB", "CG", "CD", "CE"], "x4": ["CG", "CD", "CE", "NZ"]},
             "MET": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "SD"],
                     "x3": ["CB", "CG", "SD", "CE"]},
             "PHE": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "CD1"]},
             "PRO": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "CD"],
                     "x3": ["CB", "CG", "CD", "N"], "x4": ["CG", "CD", "N", "CA"]},
             "SER": {"x1": ["N", "CA", "CB", "OG"]},
             "THR": {"x1": ["N", "CA", "CB", "OG1"]},
             "TRP": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "CD1"]},
             "TYR": {"x1": ["N", "CA", "CB", "CG"], "x2": ["CA", "CB", "CG", "CD1"]},
             "VAL": {"x1": ["N", "CA", "CB", "CG1"]}}


hard_sphere_radii = {'C': 1.6, 'O': 1.3, 'N': 1.3, 'S': 1.7}


def calculate_chi_angles(residue):
    """Calculate the chi angles for a given residue"""
    angles = chi_atoms[residue.name]
    chi = {"x1": None, "x2": None, "x3": None, "x4": None}
    for x in list(angles.keys()):
        xa1 = residue.get_atom(angles[x][0])
        xa2 = residue.get_atom(angles[x][1])
        xa3 = residue.get_atom(angles[x][2])
        xa4 = residue.get_atom(angles[x][3])
        if not xa1 or not xa2 or not xa3 or not xa4:
            continue
        xv1 = Bio.PDB.Vector(xa1.get_coordinates())
        xv2 = Bio.PDB.Vector(xa2.get_coordinates())
        xv3 = Bio.PDB.Vector(xa3.get_coordinates())
        xv4 = Bio.PDB.Vector(xa4.get_coordinates())
        chi[x] = Bio.PDB.calc_dihedral(xv1, xv2, xv3, xv4)
    return chi


def calculate_phi_psi(residue, predecessor):
    """This method calculates the PSI and PHI angles for the proteins backbone"""
    n = Bio.PDB.Vector(residue.get_atom('N').get_coordinates())
    ca = Bio.PDB.Vector(residue.get_atom('CA').get_coordinates())
    c = Bio.PDB.Vector(residue.get_atom('C').get_coordinates())
    n_pd = Bio.PDB.Vector(predecessor.get_atom('N').get_coordinates())
    ca_pd = Bio.PDB.Vector(predecessor.get_atom('CA').get_coordinates())
    c_pd = Bio.PDB.Vector(predecessor.get_atom('C').get_coordinates())

    psi = Bio.PDB.calc_dihedral(n_pd, ca_pd, c_pd, n)

    phi = Bio.PDB.calc_dihedral(c_pd, n, ca, c)

    return phi, psi


def get_interface_residues(receptor, ligand, min_cutoff=0., max_cutoff=10.):
    """Gets the interface residues depending on the minimum and maximum cutoff distances"""
    # TODO: use not hardcoded structure
    dist_matrix = scipy.spatial.distance.cdist(receptor.coordinates[0], ligand.coordinates[0])
    atom_indexes = np.where((dist_matrix >= min_cutoff) & (dist_matrix <= max_cutoff))
    receptor_residues = []
    ligand_residues = []
    for i, j in zip(atom_indexes[0], atom_indexes[1]):
        rec_atom = receptor.objects[i]
        lig_atom = ligand.objects[j]
        receptor_residues.append(rec_atom.residue_index)
        ligand_residues.append(lig_atom.residue_index)
    return list(set(receptor_residues)), list(set(ligand_residues)), dist_matrix, atom_indexes


def steric_energy(interatom_distance, atom1_element, atom2_element):
    """Calculates a simple steric energetic term taken from SCWRL3 paper "A graph-theory
    algorithm for rapid protein side-chain prediction" Protein Sci Sep 2003; 12(9) 2001-2014
    """
    sum_radii = float(hard_sphere_radii[atom1_element] + hard_sphere_radii[atom2_element])
    if interatom_distance > sum_radii:
        return 0.
    if interatom_distance < (0.8254 * sum_radii):
        return 10.
    return 57.273 * (1 - interatom_distance / sum_radii) * 0.8254 * sum_radii
