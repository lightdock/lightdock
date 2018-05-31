#!/usr/bin/env python

"""Calculates the minimum ellipsoid containing a given structure in order to free memory for
the main simulation.

This step is recommended before starting the simulation.
"""

import argparse
import numpy as np
from lightdock.util.parser import SetupCommandLineParser
from lightdock.prep.simulation import read_input_structure, save_lightdock_structure, \
                                      calculate_starting_positions, prepare_results_environment
from lightdock.constants import DEFAULT_LIGHTDOCK_PREFIX, DEFAULT_ELLIPSOID_DATA_EXTENSION
from lightdock.mathutil.ellipsoid import MinimumVolumeEllipsoid
from lightdock.util.logger import LoggingManager
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.error.lightdock_errors import LightDockError


log = LoggingManager.get_logger('lightdock_setup')


if __name__ == "__main__":

    try:
        parser = SetupCommandLineParser()
        args = parser.args

        # Read input structures
        receptor = read_input_structure(args.receptor_pdb, args.noxt)
        ligand = read_input_structure(args.ligand_pdb, args.noxt)
        
        # Move structures to origin
        receptor.move_to_origin()
        ligand.move_to_origin()

        # Calculate reference points for receptor
        log.info("Calculating reference points for receptor %s..." % args.receptor_pdb)
        rec_ellipsoid = MinimumVolumeEllipsoid(receptor.representative().coordinates)
        ellipsoid_data_file = "%s%s" % (DEFAULT_LIGHTDOCK_PREFIX % receptor.structure_file_names[0],
                                        DEFAULT_ELLIPSOID_DATA_EXTENSION)
        np.save(ellipsoid_data_file, np.array([rec_ellipsoid.center.copy()]))
        log.info("Done.")

        # Calculate reference points for ligand
        log.info("Calculating reference points for ligand %s..." % args.ligand_pdb)
        lig_ellipsoid = MinimumVolumeEllipsoid(ligand.representative().coordinates)
        ellipsoid_data_file = "%s%s" % (DEFAULT_LIGHTDOCK_PREFIX % ligand.structure_file_names[0],
                                        DEFAULT_ELLIPSOID_DATA_EXTENSION)
        np.save(ellipsoid_data_file, np.array([lig_ellipsoid.center.copy()]))
        log.info("Done.")

        # Save to file parsed structures
        save_lightdock_structure(receptor)
        save_lightdock_structure(ligand)

        # Calculate surface points (swarm centers) over receptor structure
        starting_points_files = calculate_starting_positions(receptor, ligand, 
                                                             args.swarms, args.glowworms, 
                                                             args.starting_points_seed, 
                                                             args.ftdock_file, args.use_anm, args.anm_seed)

        # Create simulation folders
        prepare_results_environment(args.swarms)

        log.info("LightDock setup OK")

    except LightDockError, error:
        log.error("LightDock setup failed. Please see:")
        log.error(error)
