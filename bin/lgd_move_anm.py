#!/usr/bin/env python3

"""Generates a set of n_confs conformations of an input pdb_file PDB structure
using the first non-trivial ANM n_modes in a given rmsd interval"""

import os
import numpy as np
import argparse
from prody import parsePDB, ANM, extendModel, confProDy, sampleModes, writePDB
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.util.logger import LoggingManager
from lightdock.constants import DEFAULT_EXTENT_MU, DEFAULT_EXTENT_SIGMA

log = LoggingManager.get_logger("lgd_move_anm")

DEFAULT_RMSD = 1.5
RANDOM_SEED = 324324


def valid_file(file_name):
    """Checks if it is a valid file"""
    if not os.path.exists(file_name):
        raise argparse.ArgumentTypeError("The file does not exist")
    return file_name


def parse_command_line():
    parser = argparse.ArgumentParser(prog="lgd_move_anm")
    parser.add_argument(
        "pdb_file",
        help="Input PDB structure",
        type=valid_file,
        metavar="pdb_file",
    )
    parser.add_argument(
        "n_modes", help="Number of ANM modes", type=int, metavar="n_modes"
    )
    parser.add_argument(
        "n_confs", help="Number of conformations to generate", type=int, metavar="n_confs"
    )
    parser.add_argument(
        "--rmsd",
        "-rmsd",
        "-r",
        help="RMSD interval", type=float, metavar="rmsd", default=DEFAULT_RMSD
    )
    parser.add_argument(
        "--ensemble",
        "-ensemble",
        "-e",
        help="Generate ensemble using native ProDy too", dest="ensemble",
        action="store_true", default=False
    )

    return parser.parse_args()


if __name__ == "__main__":
    confProDy(verbosity="info")

    # Parse arguments
    args = parse_command_line()

    molecule = parsePDB(args.pdb_file)
    atom_nodes = molecule.select("name CA")
    if not atom_nodes:
        atom_nodes = molecule.select("name P")
        if not atom_nodes:
            atom_nodes = molecule.select("name BB")

    if not atom_nodes:
        raise SystemExit("Error: Cannot find nodes for calculating ANM")

    molecule_anm = ANM("molecule nodes")
    molecule_anm.buildHessian(atom_nodes)
    molecule_anm.calcModes(n_modes=args.n_modes)
    log.info("ANM calculated")

    atoms, residues, chains = parse_complex_from_file(args.pdb_file)
    lightdock_structures = [
        {
            "atoms": atoms,
            "residues": residues,
            "chains": chains,
            "file_name": args.pdb_file,
        }
    ]
    lightdock_structure = Complex.from_structures(lightdock_structures)
    log.info("Structure parsed by LightDock")

    molecule_anm_ext, molecule_all = extendModel(
        molecule_anm, atom_nodes, molecule, norm=True
    )

    num_atoms_prody = molecule_anm_ext.numAtoms()
    num_atoms_lightdock = len(atoms)

    if num_atoms_prody != num_atoms_lightdock:
        log.error(f"Number of atoms is different between ProDy ({num_atoms_prody}) and LightDock ({num_atoms_lightdock})")
        raise SystemExit

    modes = molecule_anm_ext
    n_modes = len(modes)
    variances = modes.getVariances()
    magnitudes = np.array([abs(mode) for mode in modes])

    if np.any(variances == 0):
        log.error("Error: One or more modes has zero variance")
        raise SystemExit

    np.random.seed(RANDOM_SEED)
    randn = np.random.normal(DEFAULT_EXTENT_MU, DEFAULT_EXTENT_SIGMA, size=(args.n_confs, n_modes))
    coef = ((randn ** 2 * variances).sum(1) ** 0.5).mean()
    scale = num_atoms_prody**0.5 * args.rmsd / coef
    scale = scale / magnitudes * variances ** 0.5

    array = modes._getArray()
    coordinates = lightdock_structure.atom_coordinates[0].coordinates
    for i in range(args.n_confs):
        conf = (array * scale * randn[i]).sum(1).reshape((num_atoms_prody, 3))
        lightdock_structure.atom_coordinates[0].coordinates = coordinates + conf
        output_file = f"anm_{i+1}_{args.pdb_file}"
        write_pdb_to_file(lightdock_structure, output_file, lightdock_structure[0])
        log.info(f"Conformation {i+1} written to [{output_file}]")

    if args.ensemble:
        np.random.seed(RANDOM_SEED)
        ensemble = sampleModes(modes, molecule, n_confs=args.n_confs, rmsd=args.rmsd)
        backbone = molecule.copy()
        backbone.addCoordset(ensemble)
        ensemble_file_name = "ensemble.pdb"
        writePDB(ensemble_file_name, backbone)
        log.info(f"Ensemble file [{ensemble_file_name}] generated by ProDy")

    log.info("Done.")
