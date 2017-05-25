#!/usr/bin/env python

"""Calculates the minimum ellipsoid containing a given structure in order to free memory for
the main simulation.

This step is recommended before starting the simulation.
"""

import argparse
import numpy as np
from lightdock.constants import DEFAULT_LIGHTDOCK_PREFIX, DEFAULT_ELLIPSOID_DATA_EXTENSION
from lightdock.mathutil.ellipsoid import MinimumVolumeEllipsoid
from lightdock.util.logger import LoggingManager
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


log = LoggingManager.get_logger('lightdock_setup')


def parse_arguments():
    arg_parser = argparse.ArgumentParser(prog='lightdock_setup')
    arg_parser.add_argument("receptor_pdb_file", help="receptor structure in PDB format",
                            type=str, metavar="receptor_pdb_file")
    arg_parser.add_argument("ligand_pdb_file", help="ligand structure in PDB format",
                            type=str, metavar="ligand_pdb_file")
    arguments = arg_parser.parse_args()
    return arguments


if __name__ == "__main__":

    args = parse_arguments()

    log.info("Reading %s receptor PDB file..." % args.receptor_pdb_file)
    atoms, residues, chains = parse_complex_from_file(args.receptor_pdb_file)
    receptor = Complex(chains, atoms, residues, structure_file_name=args.receptor_pdb_file)
    log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    log.info("Reading %s ligand PDB file..." % args.ligand_pdb_file)
    atoms, residues, chains = parse_complex_from_file(args.ligand_pdb_file)
    ligand = Complex(chains, atoms, residues, structure_file_name=args.ligand_pdb_file)
    log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    # Move structures to origin
    receptor.move_to_origin()
    ligand.move_to_origin()

    log.info("Calculating reference points for receptor %s..." % args.receptor_pdb_file)
    rec_ellipsoid = MinimumVolumeEllipsoid(receptor.representative().coordinates)
    ellipsoid_data_file = "%s%s" % (DEFAULT_LIGHTDOCK_PREFIX % receptor.structure_file_names[0],
                                    DEFAULT_ELLIPSOID_DATA_EXTENSION)
    np.save(ellipsoid_data_file, np.array([rec_ellipsoid.center.copy()]))
    log.info("Done.")

    log.info("Calculating reference points for ligand %s..." % args.ligand_pdb_file)
    lig_ellipsoid = MinimumVolumeEllipsoid(ligand.representative().coordinates)
    ellipsoid_data_file = "%s%s" % (DEFAULT_LIGHTDOCK_PREFIX % ligand.structure_file_names[0],
                                    DEFAULT_ELLIPSOID_DATA_EXTENSION)
    np.save(ellipsoid_data_file, np.array([lig_ellipsoid.center.copy()]))
    log.info("Done.")
