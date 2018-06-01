#!/usr/bin/env python

"""Generates the trajectory of a given glowworm in a swarm simulation"""

import sys
import os
import numpy as np
from lightdock.constants import DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG, DEFAULT_REC_NM_FILE, DEFAULT_LIG_NM_FILE
from lightdock.structure.nm import read_nmodes
from lightdock.util.logger import LoggingManager
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.error.lightdock_errors import LightDockError
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex

log = LoggingManager.get_logger('trajectory')


def usage():
    """Displays usage parameters and exists with error"""
    log.error("Wrong arguments")
    raise SystemExit("usage: %s glowworm_id steps receptor_pdb ligand_pdb" % (sys.argv[0]))


def parse_command_line():
    # Arguments parsing
    if len(sys.argv) != 5:
        usage()
    try:
        glowworm_id = int(sys.argv[1])
    except Exception, e:
        log.error(str(e))
    
    try:
        num_steps = int(sys.argv[2])
    except Exception, e:
        log.error(str(e))
        
    if os.path.isfile(sys.argv[3]):
        receptor_pdb = sys.argv[3]
    else:
        raise LightDockError("Wrong receptor PDB file")
    
    if os.path.isfile(sys.argv[4]):
        ligand_pdb = sys.argv[4]
    else:
        raise LightDockError("Wrong ligand PDB file")
        
    return glowworm_id, num_steps, receptor_pdb, ligand_pdb


def parse_file(file_name, glowworm_id):
    """Parses a given GSO step output file. Returns the translation, the quaternion
    and the normal modes (if any)
    """
    try:
        lines = open(file_name).readlines()[1:]
    except IOError:
        return None
    
    for line_id, line in enumerate(lines):
        if line_id == glowworm_id:
            if line[0] == '(':
                coord = line[1:].split(')')[0].split(',')
                translation = [float(coord[0]), float(coord[1]), float(coord[2])]
                q = Quaternion(float(coord[3]), float(coord[4]), float(coord[5]), float(coord[6]))
                nm_rec = None
                nm_lig = None
                if len(coord) > 7:
                    nm_rec = np.array([float(x) for x in coord[7:7+DEFAULT_NMODES_REC]])
                    nm_lig = np.array([float(x) for x in coord[7:7+DEFAULT_NMODES_LIG]])
                return translation, q, nm_rec, nm_lig
    return None


if __name__ == "__main__":
    # Parse arguments
    glowworm_id, num_steps, receptor_pdb, ligand_pdb = parse_command_line()
    
    # Read receptor
    log.info("Reading %s receptor PDB file..." % receptor_pdb)
    atoms, residues, chains = parse_complex_from_file(receptor_pdb)
    receptor = Complex(chains, atoms)
    log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))
    
    # Read ligand
    log.info("Reading %s ligand PDB file..." % ligand_pdb)
    atoms, residues, chains = parse_complex_from_file(ligand_pdb)
    ligand = Complex(chains, atoms)
    log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    try:
        nm_path = os.path.abspath(os.path.dirname(receptor_pdb))
        nmodes_rec = read_nmodes(os.path.join(nm_path, DEFAULT_REC_NM_FILE + '.npy'))
    except:
        nmodes_rec = None
    try:
        nm_path = os.path.abspath(os.path.dirname(ligand_pdb))
        nmodes_lig = read_nmodes(os.path.join(nm_path, DEFAULT_LIG_NM_FILE + '.npy'))
    except:
        nmodes_lig = None

    for step in xrange(0, num_steps+1):
        try:
            # Parse each stored step file
            file_name = 'gso_%d.out' % step
            translation, rotation, nm_rec, nm_lig = parse_file(file_name, glowworm_id)

            receptor_pose = receptor.atom_coordinates[0].clone()
            ligand_pose = ligand.atom_coordinates[0].clone()

            num_rec_nmodes = len(nm_rec)
            if num_rec_nmodes:
                for i in range(num_rec_nmodes):
                    receptor_pose.coordinates += nmodes_rec[i] * nm_rec[i]
            num_lig_nmodes = len(nm_lig)
            if num_lig_nmodes:
                for i in range(num_lig_nmodes):
                    ligand_pose.coordinates += nmodes_lig[i] * nm_lig[i]

            # We rotate first, ligand it's at initial position
            ligand_pose.rotate(rotation)
            ligand_pose.translate(translation)

            output_file_name = 'trajectory_%s_step_%s.pdb' % (glowworm_id, step)
            write_pdb_to_file(receptor, output_file_name, receptor_pose)
            write_pdb_to_file(ligand, output_file_name, ligand_pose)
            log.info("Generated trajectory for step %s" % (step))
        except TypeError:
            pass
