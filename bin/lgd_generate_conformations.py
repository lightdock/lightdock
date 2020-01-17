#!/usr/bin/env python3

"""Generates the PDB structures given a LightDock swarm results file"""

import argparse
import os
import numpy as np
from lightdock.error.lightdock_errors import LightDockError
from lightdock.util.logger import LoggingManager
from lightdock.constants import DEFAULT_LIST_EXTENSION, DEFAULT_LIGHTDOCK_PREFIX, \
    DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG, DEFAULT_REC_NM_FILE, DEFAULT_LIG_NM_FILE
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.structure.nm import read_nmodes
from lightdock.prep.simulation import get_setup_from_file
from lightdock.util.parser import CommandLineParser, get_lightdock_structures


log = LoggingManager.get_logger('generate_conformations')


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
        if line[0] == '(':
            counter += 1
            last = line.index(')')
            coord = line[1:last].split(',')
            translations.append([float(coord[0]), float(coord[1]), float(coord[2])])
            rotations.append(Quaternion(float(coord[3]), float(coord[4]), float(coord[5]), float(coord[6])))
            if len(coord) > 7:
                rec_extents.append(np.array([float(x) for x in coord[7:7+num_anm_rec]]))
                lig_extents.append(np.array([float(x) for x in coord[-num_anm_lig:]]))
            raw_data = line[last+1:].split()
            receptor_id = int(raw_data[0])
            ligand_id = int(raw_data[1])
            receptor_ids.append(receptor_id)
            ligand_ids.append(ligand_id)
    log.info("Read %s coordinate lines" % counter)
    return translations, rotations, receptor_ids, ligand_ids, rec_extents, lig_extents


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="conformer_conformations")
    # Receptor
    parser.add_argument("receptor_structures", help="receptor structures: PDB file or list of PDB files",
                        type=CommandLineParser.valid_file, metavar="receptor_structure")
    # Ligand
    parser.add_argument("ligand_structures", help="ligand structures: PDB file or list of PDB files",
                        type=CommandLineParser.valid_file, metavar="ligand_structure")
    # Lightdock output file
    parser.add_argument("lightdock_output", help="lightdock output file",
                        type=CommandLineParser.valid_file, metavar="lightdock_output")
    # Number of glowworms
    parser.add_argument("glowworms", help="number of glowworms", type=CommandLineParser.valid_integer_number)
    # Optional, setup file
    parser.add_argument("--setup", "-setup", "-s", help="Simulation setup file",
                            dest="setup_file", metavar="setup_file", type=CommandLineParser.valid_file, 
                            default=None)

    args = parser.parse_args()

    # Load setup configuration if provided
    setup = get_setup_from_file(args.setup_file) if args.setup_file else None

    num_anm_rec = DEFAULT_NMODES_REC
    num_anm_lig = DEFAULT_NMODES_LIG
    if setup and setup['use_anm']:
        num_anm_rec = setup['anm_rec']
        num_anm_lig = setup['anm_lig']

    # Receptor
    structures = []
    for structure in get_lightdock_structures(args.receptor_structures):
        log.info("Reading %s receptor PDB file..." % structure)
        atoms, residues, chains = parse_complex_from_file(structure)
        structures.append({'atoms': atoms, 'residues': residues, 'chains': chains, 'file_name': structure})
        log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))
    receptor = Complex.from_structures(structures)

    # Ligand
    structures = []
    for structure in get_lightdock_structures(args.ligand_structures):
        log.info("Reading %s ligand PDB file..." % structure)
        atoms, residues, chains = parse_complex_from_file(structure)
        structures.append({'atoms': atoms, 'residues': residues, 'chains': chains, 'file_name': structure})
        log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))
    ligand = Complex.from_structures(structures)

    # Output file
    translations, rotations, receptor_ids, ligand_ids, \
        rec_extents, lig_extents = parse_output_file(args.lightdock_output, num_anm_rec, num_anm_lig)

    found_conformations = len(translations)
    num_conformations = args.glowworms
    if num_conformations > found_conformations:
        log.warning("Number of conformations is bigger than found solutions (%s > %s)" % (num_conformations,
                                                                                          found_conformations))
        log.warning("Clipping number of conformations to %s" % found_conformations)
        num_conformations = found_conformations

    # Destination path is the same as the lightdock output
    destination_path = os.path.dirname(args.lightdock_output)

    # If normal modes used, need to read them
    nmodes_rec = nmodes_lig = None
    if len(rec_extents):
        nm_path = os.path.abspath(os.path.dirname(args.receptor_structures))
        nmodes_rec = read_nmodes(os.path.join(nm_path, DEFAULT_REC_NM_FILE + '.npy'))
    if len(lig_extents):
        nm_path = os.path.abspath(os.path.dirname(args.ligand_structures))
        nmodes_lig = read_nmodes(os.path.join(nm_path, DEFAULT_LIG_NM_FILE + '.npy'))

    for i in range(num_conformations):
        receptor_pose = receptor.atom_coordinates[receptor_ids[i]].clone()
        ligand_pose = ligand.atom_coordinates[ligand_ids[i]].clone()

        # Use normal modes if provided:
        if len(rec_extents):
            try:
                for nm in range(num_anm_rec):
                    receptor_pose.coordinates += nmodes_rec[nm] * rec_extents[i][nm]
            except ValueError:
                log.error("Problem found on calculating ANM for receptor:")
                log.error("Number of atom coordinates is: %s" % str(receptor_pose.coordinates.shape))
                log.error("Number of ANM is: %s" % str(nmodes_rec.shape))
                raise SystemExit
            except IndexError:
                log.error("Problem found on calculating ANM for receptor:")
                log.error("If you have used anm_rec different than default, please use --setup")
                raise SystemExit
        if len(lig_extents):
            try:
                for nm in range(num_anm_lig):
                    ligand_pose.coordinates += nmodes_lig[nm] * lig_extents[i][nm]
            except ValueError:
                log.error("Problem found on calculating ANM for ligand:")
                log.error("Number of atom coordinates is: %s" % str(receptor_pose.coordinates.shape))
                log.error("Number of ANM is: %s" % str(nmodes_rec.shape))
                raise SystemExit
            except IndexError:
                log.error("Problem found on calculating ANM for ligand:")
                log.error("If you have used anm_lig different than default, please use --setup")
                raise SystemExit

        # We rotate first, ligand it's at initial position
        ligand_pose.rotate(rotations[i])
        ligand_pose.translate(translations[i])

        write_pdb_to_file(receptor, os.path.join(destination_path, 'lightdock_%s.pdb' % i), receptor_pose)
        write_pdb_to_file(ligand,  os.path.join(destination_path, 'lightdock_%s.pdb' % i), ligand_pose)
    log.info("Generated %d conformations" % num_conformations)
