#!/usr/bin/env python3

"""Generates the simulated trajectory of a given glowworm in a swarm"""

import argparse
import os
import numpy as np
from lightdock.constants import DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG, DEFAULT_REC_NM_FILE, DEFAULT_LIG_NM_FILE
from lightdock.structure.nm import read_nmodes
from lightdock.util.logger import LoggingManager
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.error.lightdock_errors import LightDockError
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.prep.simulation import get_setup_from_file


log = LoggingManager.get_logger('lgd_generate_trajectory')


def valid_file(file_name):
    """Checks if it is a valid file"""
    if not os.path.exists(file_name):
        raise argparse.ArgumentTypeError("The file does not exist")
    return file_name


def parse_command_line():
    parser = argparse.ArgumentParser(prog='lgd_generate_trajectory')
    parser.add_argument("glowworm_id", help="glowworm to consider", type=int, metavar="glowworm_id")
    parser.add_argument("steps", help="steps to consider", type=int, metavar="steps")
    parser.add_argument("receptor_pdb", help="Receptor LightDock parsed PDB structure", type=valid_file, 
                        metavar="receptor_pdb")
    parser.add_argument("ligand_pdb", help="Ligand LightDock parsed PDB structure", type=valid_file, 
                        metavar="ligand_pdb")
    parser.add_argument("setup_file", help="Simulation setup file", 
                        metavar="setup_file", type=valid_file)
    return parser.parse_args()


def parse_output_file(lightdock_output, glowworm_id, num_anm_rec, num_anm_lig):
    """Parses a LightDock simulation output file and returns only data for given glowworm_id"""
    with open(lightdock_output) as data_file:
        lines = data_file.readlines()
        num_glowworms = 0
        for line in lines:
            line = line.rstrip(os.linesep)
            if line.startswith('('):
                if num_glowworms == glowworm_id:
                    last = line.index(')')
                    coord = line[1:last].split(',')
                    translation = [float(coord[0]), float(coord[1]), float(coord[2])]
                    rotation = Quaternion(float(coord[3]), float(coord[4]), float(coord[5]), float(coord[6]))
                    rec_extent = None
                    lig_extent = None
                    if len(coord) > 7:
                        rec_extent = np.array([float(x) for x in coord[7:7+num_anm_rec]])
                        lig_extent = np.array([float(x) for x in coord[-num_anm_lig:]])
                    return translation, rotation, rec_extent, lig_extent
                num_glowworms += 1
    return None


if __name__ == "__main__":
    # Parse arguments
    args = parse_command_line()

    # Load setup configuration if provided
    setup = get_setup_from_file(args.setup_file) if args.setup_file else None

    num_anm_rec = DEFAULT_NMODES_REC
    num_anm_lig = DEFAULT_NMODES_LIG
    if setup and setup['use_anm']:
        num_anm_rec = setup['anm_rec']
        num_anm_lig = setup['anm_lig']
    
    # Read receptor
    log.info("Reading %s receptor PDB file..." % args.receptor_pdb)
    atoms, residues, chains = parse_complex_from_file(args.receptor_pdb)
    receptor = Complex(chains, atoms)
    log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))
    
    # Read ligand
    log.info("Reading %s ligand PDB file..." % args.ligand_pdb)
    atoms, residues, chains = parse_complex_from_file(args.ligand_pdb)
    ligand = Complex(chains, atoms)
    log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    try:
        nm_path = os.path.abspath(os.path.dirname(args.receptor_pdb))
        nmodes_rec = read_nmodes(os.path.join(nm_path, DEFAULT_REC_NM_FILE + '.npy'))
    except:
        nmodes_rec = None
    try:
        nm_path = os.path.abspath(os.path.dirname(args.ligand_pdb))
        nmodes_lig = read_nmodes(os.path.join(nm_path, DEFAULT_LIG_NM_FILE + '.npy'))
    except:
        nmodes_lig = None

    for step in range(0, args.steps+1):
        try:
            # Parse each stored step file
            file_name = 'gso_%d.out' % step
            if os.path.exists(file_name):
                translation, rotation, rec_extent, lig_extent = parse_output_file(file_name, args.glowworm_id, 
                                                                                  num_anm_rec, num_anm_lig)
                receptor_pose = receptor.atom_coordinates[0].clone()
                ligand_pose = ligand.atom_coordinates[0].clone()

                if nmodes_rec is not None:
                    try:
                        for nm in range(num_anm_rec):
                            receptor_pose.coordinates += nmodes_rec[nm] * rec_extent[nm]
                    except ValueError:
                        log.error("Problem found on calculating ANM for receptor:")
                        log.error("Number of atom coordinates is: %s" % str(receptor_pose.coordinates.shape))
                        log.error("Number of ANM is: %s" % str(nmodes_rec.shape))
                        raise SystemExit
                    except IndexError:
                        log.error("Problem found on calculating ANM for receptor:")
                        log.error("If you have used anm_rec different than default, please use --setup")
                        raise SystemExit
                if nmodes_lig is not None:
                    try:
                        for nm in range(num_anm_lig):
                            ligand_pose.coordinates += nmodes_lig[nm] * lig_extent[nm]
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
                ligand_pose.rotate(rotation)
                ligand_pose.translate(translation)

                output_file_name = 'trajectory_%s_step_%s.pdb' % (args.glowworm_id, step)
                write_pdb_to_file(receptor, output_file_name, receptor_pose)
                write_pdb_to_file(ligand, output_file_name, ligand_pose)
                log.info("Generated trajectory for step %s" % (step))
        except IOError:
            # Ignore not generated steps
            pass
