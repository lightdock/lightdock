#!/usr/bin/env python

"""Calculates the density of surface points"""

import argparse

from lightdock.constants import MIN_SURFACE_DENSITY
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.util.logger import LoggingManager
from scipy import spatial
import numpy as np


log = LoggingManager.get_logger('surface_density')


def parse_command_line():
    """
    Parses command line arguments
    """
    parser = argparse.ArgumentParser(prog='surface_density')
    parser.add_argument("pdb1", help="PDB file for receptor structure")
    parser.add_argument("pdb2", help="PDB file for ligand structure")
    parser.add_argument("points", type=int, default=400, help="The number of points on the surface")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_command_line()

    # Read receptor and calculate max radius
    atoms, residues, chains = parse_complex_from_file(args.pdb1)
    structure = Complex(chains, atoms, structure_file_name=args.pdb1)
    distances_matrix = spatial.distance.squareform(spatial.distance.pdist(structure.representative()))
    radius1 = np.max(distances_matrix)/2.

    # Read ligand and calculate max radius
    atoms, residues, chains = parse_complex_from_file(args.pdb2)
    structure = Complex(chains, atoms, structure_file_name=args.pdb2)
    distances_matrix = spatial.distance.squareform(spatial.distance.pdist(structure.representative()))
    radius2 = np.max(distances_matrix)/2.

    # Calculate the area of the sphere of radius (Rl + Rr)
    density_area = (4*np.pi*(radius1+radius2)**2)/args.points

    if density_area > MIN_SURFACE_DENSITY:
        log.warning("Surface density is below recommended, please increase the number of points on the surface.")

    print(';'.join([str(x) for x in [radius1, radius2, density_area]]))

