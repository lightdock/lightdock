#!/usr/bin/env python3

"""Calculates the diameter of a given PDB structure"""

import argparse
from scipy import spatial
import numpy as np
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.util.logger import LoggingManager


log = LoggingManager.get_logger('diameter')


def parse_command_line():
    parser = argparse.ArgumentParser(prog='calculate_diameter')
    parser.add_argument("pdb", help="PDB file for structure to calculate maximum diameter")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_command_line()

    atoms, residues, chains = parse_complex_from_file(args.pdb)
    structure = Complex(chains, atoms, structure_file_name=args.pdb)
    distances_matrix = spatial.distance.squareform(spatial.distance.pdist(structure.representative()))
    ligand_max_diameter = np.max(distances_matrix)

    print(ligand_max_diameter)
