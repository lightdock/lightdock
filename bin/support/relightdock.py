#!/user/bin/env python

"""Prototype to re-start lightdock calculations from a gso output file.

This is high-experimental, this code can change without any warning."""

import os
import sys
from lightdock.util.logger import LoggingManager
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.gso.boundaries import Boundary, BoundingBox
from lightdock.mathutil.lrandom import MTGenerator
from lightdock.gso.parameters import GSOParameters
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.scoring.dfire.driver import DFIRE, DFIREAdapter
from lightdock.gso.algorithm import GSO
from lightdock.gso.searchspace.landscape import DockingLandscapePosition
from lightdock.gso.swarm import Swarm


log = LoggingManager.get_logger('relightdock')

MAX_TRANSLATION = 30                # In Angstroms
MAX_ROTATION = 1.0                  # Quaternion default value for its components
DEFAULT_TRANSLATION_STEP = 0.5      # Normalized step
DEFAULT_ROTATION_STEP = 0.5         # Normalized SLERP step. 1 means full jump, 0 means no movement
GSO_SEED = 324324                   # Seed for the random number generator in the GSO algorithm
STARTING_POINTS_SEED = 324324


def parse_output_file(lightdock_output):
    translations = []
    rotations = []
    luciferin = []
    neighbors = []
    vision_range = []
    scoring = []

    data_file = open(lightdock_output)
    lines = [line.rstrip(os.linesep) for line in data_file.readlines()]
    data_file.close()

    counter = 0
    for line in lines:
        if line[0] == '(':
            counter += 1
            last = line.index(')')
            coord = line[1:last].split(',')
            translations.append([float(coord[0]), float(coord[1]), float(coord[2])])
            rotations.append(Quaternion(float(coord[3]), float(coord[4]), float(coord[5]), float(coord[6])))
            values = line[last + 1:].split()
            luciferin.append(float(values[0]))
            neighbors.append(int(values[1]))
            vision_range.append(float(values[2]))
            scoring.append(float(values[3]))

    log.info("Read %s coordinate lines" % counter)
    return translations, rotations, luciferin, neighbors, vision_range, scoring


if __name__ == "__main__":
    starting_file = sys.argv[1]
    receptor_file = sys.argv[2]
    ligand_file = sys.argv[3]
    steps = int(sys.argv[4])
    configuration_file = sys.argv[5]

    log.info("Starting file: %s" % starting_file)
    log.info("Receptor: %s" % receptor_file)
    log.info("Ligand: %s" % ligand_file)
    log.info("Steps: %d" % steps)
    log.info("Configuration file: %s" % configuration_file)
    print

    # Read structures (already in the center)
    log.info("Reading %s receptor PDB file..." % receptor_file)
    atoms, residues, chains = parse_complex_from_file(receptor_file)
    receptor = Complex(chains, atoms)
    log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    log.info("Reading %s ligand PDB file..." % ligand_file)
    atoms, residues, chains = parse_complex_from_file(ligand_file)
    ligand = Complex(chains, atoms)
    log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    # Start from results positions
    log.info("Reading calculated data from %s" % starting_file)
    translations, rotations, luciferin, neighbors, vision_range, scoring = parse_output_file(starting_file)
    num_glowworms = len(translations)
    log.info("%d glowworms loaded" % num_glowworms)

    adapter = DFIREAdapter(receptor, ligand)
    scoring_function = DFIRE()
    log.info("Loaded DFIRE scoring function")

    bounding_box = BoundingBox([Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                                Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                                Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                                Boundary(-MAX_ROTATION, MAX_ROTATION),
                                Boundary(-MAX_ROTATION, MAX_ROTATION),
                                Boundary(-MAX_ROTATION, MAX_ROTATION),
                                Boundary(-MAX_ROTATION, MAX_ROTATION)])
    random_number_generator = MTGenerator(GSO_SEED)
    gso_parameters = GSOParameters(configuration_file)
    log.info("Parameters read")

    positions = []
    for i in range(num_glowworms):
        coordinates = [translations[i][0], translations[i][1], translations[i][2], rotations[i].w, rotations[i].x, rotations[i].y, rotations[i].z]
        positions.append(DockingLandscapePosition(scoring_function, coordinates, adapter.receptor_model, adapter.ligand_model))
    log.info("%d positions loaded" % len(positions))

    swarm = Swarm(positions, gso_parameters)
    for i, glowworm in enumerate(swarm.glowworms):
        glowworm.luciferin = luciferin[i]
        glowworm.vision_range = vision_range[i]
        glowworm.scoring = scoring[i]
    for i, glowworm in enumerate(swarm.glowworms):
        glowworm.search_neighbors(swarm.glowworms)
    log.info("Swarm created")

    gso = GSO(swarm, gso_parameters, random_number_generator)

    log.info("Starting simulation")
    gso.run(steps, verbose=True, save_intermediary=True)
    log.info("Finished.")
