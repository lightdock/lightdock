#!/usr/bin/env python3

"""Calculates the new PDB structure of an input PDB structure using ANM"""

import sys
from prody import parsePDB, ANM, extendModel, confProDy
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex


def usage():
    print("Usage: %s PDB_file n_modes factor" % sys.argv[0])


if __name__ == "__main__":
    confProDy(verbosity='info')

    try:
        pdb_structure = sys.argv[1]
        n_modes = int(sys.argv[2])
        factor = float(sys.argv[3])
    except:
        usage()
        raise SystemExit("Wrong command line")

    protein = parsePDB(pdb_structure)
    ca_atoms = protein.select('name CA')
    protein_anm = ANM('protein ca')
    protein_anm.buildHessian(ca_atoms)
    protein_anm.calcModes(n_modes=n_modes)
    print('Normal modes calculated')

    atoms, residues, chains = parse_complex_from_file(pdb_structure)
    lightdock_structures = [{'atoms': atoms, 'residues': residues, 'chains': chains, 'file_name': pdb_structure}]
    lightdock_structure = Complex.from_structures(lightdock_structures)
    print('Structure read by lightdock')

    num_atoms_prody = len(protein.protein)
    num_atoms_lightdock = len(atoms)

    if num_atoms_prody != num_atoms_lightdock:
        raise SystemExit("Number of atoms is different")

    protein_anm_ext, protein_all = extendModel(protein_anm, ca_atoms, protein, norm=True)

    modes = []
    for i in range(n_modes):
        nm = protein_anm_ext.getEigvecs()[:, i].reshape((num_atoms_lightdock, 3))
        modes.append(nm)

    coordinates = lightdock_structure.atom_coordinates[0].coordinates
    for i in range(n_modes):
        lightdock_structure.atom_coordinates[0].coordinates += modes[i] * factor

    output_file = 'anm_' + pdb_structure
    write_pdb_to_file(lightdock_structure, output_file, lightdock_structure[0])

    print('Structure written to %s' % output_file)
    print('Done.')
