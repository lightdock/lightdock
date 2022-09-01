#!/usr/bin/env python3

"""Generates the PDB structures given a LightDock swarm results file"""

import argparse
import os
import numpy as np
from pathlib import Path
from lightdock.util.logger import LoggingManager
from lightdock.constants import (
    DEFAULT_NMODES_REC,
    DEFAULT_NMODES_LIG,
    DEFAULT_REC_NM_FILE,
    DEFAULT_LIG_NM_FILE,
)
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.structure.nm import read_nmodes
from lightdock.prep.simulation import get_setup_from_file
from lightdock.util.parser import (
    valid_file,
    valid_integer_number,
    get_lightdock_structures,
)


log = LoggingManager.get_logger("generate_conformations")


def parse_output_file(lightdock_output, num_anm_rec, num_anm_lig):
    translations = []
    rotations = []
    receptor_ids = []
    ligand_ids = []
    rec_extents = []
    lig_extents = []

    data_file = open(lightdock_output)
    lines = data_file.readlines()
    data_file.close()

    counter = 0
    for line in lines:
        if line[0] == "(":
            counter += 1
            last = line.index(")")
            coord = line[1:last].split(",")
            translations.append([float(coord[0]), float(coord[1]), float(coord[2])])
            rotations.append(
                Quaternion(
                    float(coord[3]), float(coord[4]), float(coord[5]), float(coord[6])
                )
            )
            if len(coord) > 7:
                rec_extents.append(
                    np.array([float(x) for x in coord[7 : 7 + num_anm_rec]])
                )
                lig_extents.append(np.array([float(x) for x in coord[-num_anm_lig:]]))
            raw_data = line[last + 1 :].split()
            receptor_id = int(raw_data[0])
            ligand_id = int(raw_data[1])
            receptor_ids.append(receptor_id)
            ligand_ids.append(ligand_id)
    log.info("Read %s coordinate lines" % counter)
    return translations, rotations, receptor_ids, ligand_ids, rec_extents, lig_extents


def parse_initial_file(lightdock_output, num_anm_rec, num_anm_lig):
    translations = []
    rotations = []
    receptor_ids = []
    ligand_ids = []
    rec_extents = []
    lig_extents = []

    data_file = open(lightdock_output)
    lines = data_file.readlines()
    data_file.close()

    counter = 0
    for line in lines:
        if line[0] != "#":
            counter += 1
            line = line.rstrip(os.linesep)
            coord = line.split()
            translations.append([float(coord[0]), float(coord[1]), float(coord[2])])
            rotations.append(
                Quaternion(
                    float(coord[3]), float(coord[4]), float(coord[5]), float(coord[6])
                )
            )
            if len(coord) > 7:
                rec_extents.append(
                    np.array([float(x) for x in coord[7 : 7 + num_anm_rec]])
                )
                lig_extents.append(np.array([float(x) for x in coord[-num_anm_lig:]]))

            receptor_ids.append(0)
            ligand_ids.append(0)

    log.info("Read %s coordinate lines" % counter)
    return translations, rotations, receptor_ids, ligand_ids, rec_extents, lig_extents


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
    # Lightdock output file
    parser.add_argument(
        "lightdock_output",
        help="lightdock output file",
        type=valid_file,
        metavar="lightdock_output",
    )
    # Number of glowworms
    parser.add_argument(
        "glowworms", help="number of glowworms", type=valid_integer_number
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

    # Output file
    output_file_path = Path(args.lightdock_output)
    if output_file_path.suffix == ".dat":
        (
            translations,
            rotations,
            receptor_ids,
            ligand_ids,
            rec_extents,
            lig_extents,
        ) = parse_initial_file(args.lightdock_output, num_anm_rec, num_anm_lig)
    else:
        (
            translations,
            rotations,
            receptor_ids,
            ligand_ids,
            rec_extents,
            lig_extents,
        ) = parse_output_file(args.lightdock_output, num_anm_rec, num_anm_lig)

    found_conformations = len(translations)
    num_conformations = args.glowworms
    if num_conformations > found_conformations:
        log.warning(
            "Number of conformations is bigger than found solutions (%s > %s)"
            % (num_conformations, found_conformations)
        )
        log.warning("Clipping number of conformations to %s" % found_conformations)
        num_conformations = found_conformations

    # Destination path is the same as the lightdock output
    destination_path = os.path.dirname(args.lightdock_output)

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

    for i in range(num_conformations):
        receptor_pose = receptor.atom_coordinates[receptor_ids[i]].clone()
        ligand_pose = ligand.atom_coordinates[ligand_ids[i]].clone()

        # Use normal modes if provided:
        if nmodes_rec is not None and nmodes_rec.any():
            try:
                for nm in range(num_anm_rec):
                    receptor_pose.coordinates[receptor.nm_mask, :] += (
                        nmodes_rec[nm] * rec_extents[i][nm]
                    )
            except ValueError:
                log.error("Problem found on calculating ANM for receptor:")
                log.error(
                    "Number of atom coordinates is: %s"
                    % str(receptor_pose.coordinates.shape)
                )
                log.error("Number of ANM is: %s" % str(nmodes_rec.shape))
                raise SystemExit
            except IndexError:
                log.error("Problem found on calculating ANM for receptor:")
                log.error(
                    "If you have used anm_rec different than default, please use --setup"
                )
                raise SystemExit
        if nmodes_lig is not None and nmodes_lig.any():
            try:
                for nm in range(num_anm_lig):
                    ligand_pose.coordinates[ligand.nm_mask, :] += (
                        nmodes_lig[nm] * lig_extents[i][nm]
                    )
            except ValueError:
                log.error("Problem found on calculating ANM for ligand:")
                log.error(
                    "Number of atom coordinates is: %s"
                    % str(receptor_pose.coordinates.shape)
                )
                log.error("Number of ANM is: %s" % str(nmodes_rec.shape))
                raise SystemExit
            except IndexError:
                log.error("Problem found on calculating ANM for ligand:")
                log.error(
                    "If you have used anm_lig different than default, please use --setup"
                )
                raise SystemExit

        # We rotate first, ligand it's at initial position
        ligand_pose.rotate(rotations[i])
        ligand_pose.translate(translations[i])

        write_pdb_to_file(
            receptor,
            os.path.join(destination_path, "lightdock_%s.pdb" % i),
            receptor_pose,
        )
        write_pdb_to_file(
            ligand, os.path.join(destination_path, "lightdock_%s.pdb" % i), ligand_pose
        )
    log.info("Generated %d conformations" % num_conformations)
