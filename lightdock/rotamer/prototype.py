import sys
import math

from lightdock.scoring.dfire.driver import DFIREAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.rotamer.predictor import get_interface_residues, calculate_chi_angles, steric_energy
from lightdock.rotamer.library import InterfaceSurfaceLibrary


def usage():
    print("Usage: %s receptor.pdb ligand.pdb" % sys.argv[0])


if __name__ == "__main__":
    if len(sys.argv[1:]) != 2:
        usage()
        raise SystemExit('Wrong command line')

    atoms, residues, chains = parse_complex_from_file(sys.argv[1])
    receptor = Complex(chains, atoms)
    atoms, residues, chains = parse_complex_from_file(sys.argv[2])
    ligand = Complex(chains, atoms)
    adapter = DFIREAdapter(receptor, ligand)

    # 5A of cutoff to be in the same line as the ASA value used in the original paper
    receptor_residue_indexes, ligand_residue_indexes, dist_matrix, atom_indexes = get_interface_residues(adapter.receptor_model,
                                                                                                         adapter.ligand_model,
                                                                                                         max_cutoff=5.)
    rotamer_library = InterfaceSurfaceLibrary()

    print("%d receptor interface residues" % len(receptor_residue_indexes))
    receptor_residues = []
    for residue_index in receptor_residue_indexes:
        residue = receptor.residues[residue_index]
        chi_angles = calculate_chi_angles(residue)
        for chi, angle in chi_angles.items():
            try:
                chi_angles[chi] = math.degrees(angle)
            except TypeError:
                pass
        index_rotamer, rotamer = rotamer_library.get_closest_rotamer(residue, chi_angles)
        if chi_angles['x1']:
            print(residue.name, chi_angles, rotamer)
        receptor_residues.append(residue)

        # Code to move side-chain to closest rotamer and check if the angles correspond
        residue.mutate_side_chain(rotamer)
        chi_angles = calculate_chi_angles(residue)
        for chi, angle in chi_angles.items():
            try:
                chi_angles[chi] = math.degrees(angle)
            except TypeError:
                pass
        print(chi_angles)
        print()

    print("%d ligand interface residues" % len(ligand_residue_indexes))
    ligand_residues = []
    for residue_index in ligand_residue_indexes:
        residue = ligand.residues[residue_index]
        chi_angles = calculate_chi_angles(residue)
        for chi, angle in chi_angles.items():
            try:
                chi_angles[chi] = math.degrees(angle)
            except TypeError:
                pass
        index_rotamer, rotamer = rotamer_library.get_closest_rotamer(residue, chi_angles)
        if chi_angles['x1']:
            print(residue.name, chi_angles, rotamer)
        ligand_residues.append(residue)

        # Code to move side-chain to closest rotamer and check if the angles correspond
        residue.mutate_side_chain(rotamer)
        chi_angles = calculate_chi_angles(residue)
        for chi, angle in chi_angles.items():
            try:
                chi_angles[chi] = math.degrees(angle)
            except TypeError:
                pass
        print(chi_angles)
        print()
    
    # Calculate steric energy
    energy = 0
    for residue_index1 in receptor_residue_indexes:
        residue1 = receptor.residues[residue_index1]
        for residue_index2 in ligand_residue_indexes:
            residue2 = ligand.residues[residue_index2]
            for atom1 in residue1.sidechain:
                for atom2 in residue2.sidechain:
                    energy += steric_energy(dist_matrix[atom1.index][atom2.index], atom1.element, atom2.element)
    print(energy)
