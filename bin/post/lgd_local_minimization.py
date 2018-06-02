#!/usr/bin/env python

"""Runs a LightDock local minimization.

This is very experimental and may change without any warning."""

import argparse
import os
import numpy as np
import importlib
from lightdock.util.logger import LoggingManager
from lightdock.constants import DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG, \
    DEFAULT_REC_NM_FILE, DEFAULT_LIG_NM_FILE, DEFAULT_SCORING_FUNCTION
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.structure.nm import read_nmodes
from scipy.optimize import fmin_powell


log = LoggingManager.get_logger('minimizer')


def valid_file(file_name):
    """Checks if it is a valid file"""
    if not os.path.exists(file_name):
        raise argparse.ArgumentTypeError("The file does not exist")
    return file_name


def valid_integer_number(ivalue):
    """Checks for a valid integer"""
    try:
        ivalue = int(ivalue)
    except:
        raise argparse.ArgumentTypeError("%s is an invalid value" % ivalue)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid value" % ivalue)
    return ivalue


def parse_output_file(lightdock_output):
    """Parses a lightdock output file to get relevant simulation information for each glowworm"""
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
            translations.append([float(coord[0]),float(coord[1]),float(coord[2])])
            rotations.append(Quaternion(float(coord[3]),float(coord[4]),float(coord[5]),float(coord[6])))
            if len(coord) == (7 + DEFAULT_NMODES_REC + DEFAULT_NMODES_LIG):
                rec_extents.append(np.array([float(x) for x in coord[7:7+DEFAULT_NMODES_REC]]))
                lig_extents.append(np.array([float(x) for x in coord[-DEFAULT_NMODES_LIG:]]))
            raw_data = line[last+1:].split()
            receptor_id = int(raw_data[0])
            ligand_id = int(raw_data[1])
            receptor_ids.append(receptor_id)
            ligand_ids.append(ligand_id)
    log.info("Read %s coordinate lines" % counter)
    return translations, rotations, receptor_ids, ligand_ids, rec_extents, lig_extents


def set_scoring_function(function_name, receptor, ligand):
    """Initializes a scoring function given the receptor and ligand"""
    log.info("Loading scoring function...")
    if function_name:
        scoring_function_module = "lightdock.scoring.%s.driver" % function_name
    else:
        scoring_function_module = "lightdock.scoring.%s.driver" % DEFAULT_SCORING_FUNCTION
    module = importlib.import_module(scoring_function_module)

    log.info("Using %s scoring function" % module.DefinedScoringFunction.__name__)

    CurrentScoringFunction = getattr(module, "DefinedScoringFunction")
    CurrentModelAdapter = getattr(module, "DefinedModelAdapter")
    adapters = CurrentModelAdapter(receptor, ligand)
    scoring_function = CurrentScoringFunction()
    log.info("Done.")
    return scoring_function, adapters


def move_molecules(translation, rotation, nm_ext_rec, nm_ext_lig):
    receptor_pose = receptor.atom_coordinates[receptor_ids[0]].clone()
    ligand_pose = ligand.atom_coordinates[ligand_ids[0]].clone()
    if len(rec_extents):
        for nm in range(DEFAULT_NMODES_REC):
            receptor_pose.coordinates += nmodes_rec[nm] * nm_ext_rec[nm]
    if len(lig_extents):
        for nm in range(DEFAULT_NMODES_LIG):
            ligand_pose.coordinates += nmodes_lig[nm] * nm_ext_lig[nm]

    # We rotate first, ligand it's at initial position
    ligand_pose.rotate(rotation)
    ligand_pose.translate(translation)

    return receptor_pose, ligand_pose


def calculate_scoring(optimization_vector):
    #print optimization_vector
    translation = optimization_vector[:3]
    rotation = optimization_vector[3:7]
    q = Quaternion(rotation[0], rotation[1], rotation[2], rotation[3])
    nm_ext_rec = optimization_vector[7:7+DEFAULT_NMODES_REC]
    nm_ext_lig = optimization_vector[-DEFAULT_NMODES_LIG:]
    receptor_pose, ligand_pose = move_molecules(translation, q, nm_ext_rec, nm_ext_lig)
    energy = -1. * scoring_function(adapters.receptor_model, receptor_pose, adapters.ligand_model, ligand_pose)
    print energy
    return energy


if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser(prog="minimizer")
    # Receptor
    parser.add_argument("receptor_structure", help="receptor structure", type=valid_file, metavar="receptor_structure")
    # Ligand
    parser.add_argument("ligand_structure", help="ligand structure", type=valid_file, metavar="ligand_structure")
    # Lightdock output file
    parser.add_argument("lightdock_output", help="lightdock output file", type=valid_file, metavar="lightdock_output")
    # Solution to minimize
    parser.add_argument("glowworm", help="glowworm to minimize", type=valid_integer_number)
    # Scoring function
    parser.add_argument("-s", "--scoring_function", help="scoring function used", dest="scoring_function")
    args = parser.parse_args()

    # Read receptor structure
    structure = args.receptor_structure
    log.info("Reading %s receptor PDB file..." % structure)
    atoms, residues, chains = parse_complex_from_file(structure)
    receptor = Complex(chains, atoms, residues, structure_file_name=structure)
    log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    # Read ligand structure
    structure = args.ligand_structure
    log.info("Reading %s ligand PDB file..." % structure)
    atoms, residues, chains = parse_complex_from_file(structure)
    ligand = Complex(chains, atoms, residues, structure_file_name=structure)
    log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    # Output file
    translations, rotations, receptor_ids, ligand_ids, rec_extents, lig_extents = parse_output_file(args.lightdock_output)

    # Destination path is the same as the lightdock output
    destination_path = os.path.dirname(args.lightdock_output)

    # If normal modes used, need to read them
    if len(rec_extents):
        nm_path = os.path.abspath(os.path.dirname(args.receptor_structure))
        nmodes_rec = read_nmodes(os.path.join(nm_path, DEFAULT_REC_NM_FILE + '.npy'))
    if len(lig_extents):
        nm_path = os.path.abspath(os.path.dirname(args.ligand_structure))
        nmodes_lig = read_nmodes(os.path.join(nm_path, DEFAULT_LIG_NM_FILE + '.npy'))

    initial_rec_extents = list([float(x) for x in rec_extents[args.glowworm]])
    initial_lig_extents = list([float(x) for x in lig_extents[args.glowworm]])

    # Prepare scoring function
    scoring_function, adapters = set_scoring_function(args.scoring_function, receptor, ligand)

    optimization_vector = []
    optimization_vector.extend(translations[args.glowworm])
    q = rotations[args.glowworm]
    optimization_vector.extend([q.w, q.x, q.y, q.z])
    optimization_vector.extend(initial_rec_extents)
    optimization_vector.extend(initial_lig_extents)
    optimization_vector = np.array(optimization_vector)

    # Minimize using Powell algorythm
    result = fmin_powell(calculate_scoring, optimization_vector, maxiter=2, full_output=1, xtol=0.05)
    print 'Target: ', optimization_vector
    print 'Minimized: ', result[0]
    print 'Energy: ', result[1]

    with open('minimized_%d.out' % args.glowworm, 'w') as output:
        coordinates = ', '.join(['%8.5f' % x for x in result[0]])
        output.write('(%s)' % (coordinates) + os.linesep)
