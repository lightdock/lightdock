#!/usr/bin/env python3

"""Generates the top N structures in PDB format given a ranking file"""

import argparse
import os
import numpy as np
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.util.analysis import read_ranking_file
from lightdock.util.logger import LoggingManager
from lightdock.constants import (
    DEFAULT_NMODES_REC,
    DEFAULT_NMODES_LIG,
    DEFAULT_REC_NM_FILE,
    DEFAULT_LIG_NM_FILE,
)
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.structure.nm import read_nmodes
from lightdock.util.parser import (
    valid_file,
    valid_integer_number,
    get_lightdock_structures,
)
from lightdock.prep.simulation import get_setup_from_file


log = LoggingManager.get_logger("lightdock_top")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="conformer_conformations")
    # Receptor
    parser.add_argument(
        "receptor_structures",
        help="receptor structures: PDB file or list of PDB files",
        type=valid_file,
        metavar="receptor_structure",
    )
    # Ligand
    parser.add_argument(
        "ligand_structures",
        help="ligand structures: PDB file or list of PDB files",
        type=valid_file,
        metavar="ligand_structure",
    )
    # Ranking file
    parser.add_argument(
        "lightdock_ranking_file",
        help="LightDock ranking file",
        type=valid_file,
        metavar="lightdock_ranking_file",
    )
    # Number of structures to generate
    parser.add_argument(
        "top",
        help="number of structures to generate",
        type=valid_integer_number,
        metavar="top",
    )
    # Optional, setup file
    parser.add_argument(
        "--setup",
        "-setup",
        "-s",
        help="Simulation setup file",
        dest="setup_file",
        metavar="setup_file",
        type=valid_file,
        default=None,
    )

    args = parser.parse_args()

    # Load setup configuration if provided
    setup = get_setup_from_file(args.setup_file) if args.setup_file else None

    num_anm_rec = DEFAULT_NMODES_REC
    num_anm_lig = DEFAULT_NMODES_LIG
    if setup and setup["use_anm"]:
        num_anm_rec = setup["anm_rec"]
        num_anm_lig = setup["anm_lig"]

    # Receptor
    structures = []
    for structure in get_lightdock_structures(args.receptor_structures):
        log.info("Reading %s receptor PDB file..." % structure)
        atoms, residues, chains = parse_complex_from_file(structure)
        structures.append(
            {
                "atoms": atoms,
                "residues": residues,
                "chains": chains,
                "file_name": structure,
            }
        )
        log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))
    receptor = Complex.from_structures(structures)

    # Ligand
    structures = []
    for structure in get_lightdock_structures(args.ligand_structures):
        log.info("Reading %s ligand PDB file..." % structure)
        atoms, residues, chains = parse_complex_from_file(structure)
        structures.append(
            {
                "atoms": atoms,
                "residues": residues,
                "chains": chains,
                "file_name": structure,
            }
        )
        log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))
    ligand = Complex.from_structures(structures)

    # Read ranking file
    predictions = read_ranking_file(args.lightdock_ranking_file)

    # Destination path is the same as the lightdock output
    destination_path = os.path.dirname(args.lightdock_ranking_file)

    # If normal modes used, need to read them
    nmodes_rec = nmodes_lig = None
    nm_path = os.path.abspath(os.path.dirname(args.receptor_structures))
    # Check NM file for receptor
    nm_rec_file = os.path.join(nm_path, DEFAULT_REC_NM_FILE + ".npy")
    if os.path.exists(nm_rec_file):
        nmodes_rec = read_nmodes(nm_rec_file)
    # Check NM file for ligand
    nm_lig_file = os.path.join(nm_path, DEFAULT_LIG_NM_FILE + ".npy")
    if os.path.exists(nm_lig_file):
        nmodes_lig = read_nmodes(nm_lig_file)

    for i, glowworm in enumerate(predictions):
        if i < args.top:
            receptor_pose = receptor.atom_coordinates[glowworm.receptor_id].clone()
            ligand_pose = ligand.atom_coordinates[glowworm.ligand_id].clone()
            # Use normal modes if provided:
            if nmodes_rec is not None and nmodes_rec.any():
                for nm in range(num_anm_rec):
                    rec_extent = np.array(
                        [float(x) for x in glowworm.pose[7 : 7 + num_anm_rec]]
                    )
                    receptor_pose.coordinates += nmodes_rec[nm] * rec_extent[nm]
            if nmodes_lig is not None and nmodes_lig.any():
                for nm in range(num_anm_lig):
                    lig_extent = np.array(
                        [float(x) for x in glowworm.pose[-num_anm_lig:]]
                    )
                    ligand_pose.coordinates += nmodes_lig[nm] * lig_extent[nm]

            # We rotate first, ligand it's at initial position
            rotation = Quaternion(
                glowworm.pose[3], glowworm.pose[4], glowworm.pose[5], glowworm.pose[6]
            )
            ligand_pose.rotate(rotation)
            ligand_pose.translate(
                [glowworm.pose[0], glowworm.pose[1], glowworm.pose[2]]
            )

            write_pdb_to_file(
                receptor,
                os.path.join(destination_path, "top_%s.pdb" % str(i + 1)),
                receptor_pose,
            )
            write_pdb_to_file(
                ligand,
                os.path.join(destination_path, "top_%s.pdb" % str(i + 1)),
                ligand_pose,
            )
    log.info("Generated %d conformations" % args.top)
