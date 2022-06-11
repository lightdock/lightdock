#!/usr/bin/env python3

"""Maps the contacts during simulation on the receptor surface"""

import argparse
import os
import numpy as np
from scipy.spatial import distance
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.util.logger import LoggingManager
from lightdock.constants import (
    DEFAULT_ATOMIC_CONTACT,
    DEFAULT_REC_NM_FILE,
    DEFAULT_LIG_NM_FILE,
)
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.structure.nm import read_nmodes
from lightdock.util.parser import (
    valid_file,
    get_lightdock_structures,
)
from lightdock.prep.simulation import get_setup_from_file


log = LoggingManager.get_logger("lgd_map_contacts")


def calculate_pose(
    translation,
    rotation,
    rec_extent,
    lig_extent,
    receptor,
    num_anm_rec,
    nmodes_rec,
    ligand,
    num_anm_lig,
    nmodes_lig
):
    """Calculates a docking pose encoded by a glowworm"""
    receptor_pose = receptor.atom_coordinates[0].clone()
    ligand_pose = ligand.atom_coordinates[0].clone()

    # Use normal modes if provided:
    if nmodes_rec is not None and nmodes_rec.any():
        for nm in range(num_anm_rec):
            receptor_pose.coordinates += nmodes_rec[nm] * rec_extent[nm]

    if nmodes_lig is not None and nmodes_lig.any():
        for nm in range(num_anm_lig):
            ligand_pose.coordinates += nmodes_lig[nm] * lig_extent[nm]

    # We rotate first, ligand it's at initial position
    ligand_pose.rotate(rotation)
    ligand_pose.translate(translation)

    return receptor_pose, ligand_pose


def parse_output_file(lightdock_output, num_anm_rec, num_anm_lig):
    """Parses a GSO LightDock output"""
    translations = []
    rotations = []
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

    log.info("Read %s coordinate lines" % counter)
    return translations, rotations, rec_extents, lig_extents


def calculate_contacts(rec_atoms_per_residue, rec_pose, lig_pose):
    dist_matrix = distance.cdist(rec_pose, lig_pose)
    atom_indexes = np.where(dist_matrix <= DEFAULT_ATOMIC_CONTACT)

    new_contacts = [rec_atoms_per_residue[i] for i in atom_indexes[0]]

    return list(set(new_contacts))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="lgd_map_contacts")
    # Receptor PDB
    parser.add_argument(
        "receptor_structure",
        help="PDB file of the receptor structure",
        type=valid_file,
        metavar="receptor_structure",
    )
    # Ligand PDB
    parser.add_argument(
        "ligand_structure",
        help="PDB file of the ligand structure",
        type=valid_file,
        metavar="ligand_structure",
    )
    # Setup file
    parser.add_argument(
        "setup_file",
        help="Simulation setup file",
        metavar="setup_file",
        type=valid_file,
    )
    # Output PDB file
    parser.add_argument(
        "output_pdb_file",
        help="Output PDB file",
        metavar="output_pdb_file"
    )
    # GSO files
    parser.add_argument(
        "gso_files",
        help="List of LightDock output files",
        nargs='+', default=[],
        metavar="gso_files",
    )

    args = parser.parse_args()

    if os.path.exists(args.output_pdb_file):
        log.error(f"File {args.output_pdb_file} already exists")
        raise SystemExit

    # Load setup configuration
    setup = get_setup_from_file(args.setup_file)

    num_anm_rec = 0
    num_anm_lig = 0
    if setup["use_anm"]:
        num_anm_rec = setup["anm_rec"]
        num_anm_lig = setup["anm_lig"]

    simulation_path = os.path.abspath(os.path.dirname(args.setup_file))

    # Receptor
    structures = []
    for structure in get_lightdock_structures(args.receptor_structure):
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
    for structure in get_lightdock_structures(args.ligand_structure):
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

    # If normal modes used, need to read them
    nmodes_rec = nmodes_lig = None
    nm_path = os.path.abspath(os.path.dirname(args.receptor_structure))
    # Check NM file for receptor
    nm_rec_file = os.path.join(nm_path, DEFAULT_REC_NM_FILE + ".npy")
    if os.path.exists(nm_rec_file):
        nmodes_rec = read_nmodes(nm_rec_file)
    # Check NM file for ligand
    nm_lig_file = os.path.join(nm_path, DEFAULT_LIG_NM_FILE + ".npy")
    if os.path.exists(nm_lig_file):
        nmodes_lig = read_nmodes(nm_lig_file)

    # Map for each atom at which residue belongs in the same order as atoms list
    rec_atoms_per_residue = [f"{a.chain_id}.{a.residue_name}.{a.residue_number}" for a in receptor.atoms]
    residue_freqs = {}
    for residue in list(set(rec_atoms_per_residue)):
        residue_freqs[residue] = 0

    # Parse GSO files
    for gso_file in args.gso_files:
        log.info(f"Reading {gso_file}")
        translations, rotations, rec_extents, lig_extents = parse_output_file(gso_file, num_anm_rec, num_anm_lig)

        # Generate pose for each glowworm found:
        num_glowworms = len(translations)
        for i in range(num_glowworms):
            print('.', end='', flush=True)
            rec_pose, lig_pose = calculate_pose(translations[i], rotations[i],
                                                rec_extents[i], lig_extents[i],
                                                receptor, num_anm_rec, nmodes_rec,
                                                ligand, num_anm_lig, nmodes_lig)
            new_contacts = calculate_contacts(rec_atoms_per_residue, rec_pose, lig_pose)
            for residue in new_contacts:
                residue_freqs[residue] += 1
        print()

    log.info("Residue frequencies:")
    print(residue_freqs)

    # Normalize frequencies
    freqs = np.array(list(residue_freqs.values()))
    norm = np.linalg.norm(freqs)
    residue_freqs_norm = {k: v / norm * 100.0 for k, v in residue_freqs.items()}
    max_freq = np.max(np.array(list(residue_freqs_norm.values())))
    print(f"spectrum b, white_blue, minimum=0, maximum={max_freq}")

    for atom in receptor.atoms:
        residue_id = f"{atom.chain_id}.{atom.residue_name}.{atom.residue_number}"
        atom.b_factor = residue_freqs_norm[residue_id]

    write_pdb_to_file(receptor, args.output_pdb_file)
