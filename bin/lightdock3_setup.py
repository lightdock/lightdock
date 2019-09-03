#!/usr/bin/env python3

"""Before launching the LightDock simulation, a setup step is required.

This step parses the input PDB structures, calculates the minimum ellipsoid 
containing each of them, calculates the swarms on the surface of the 
receptor and populates each swarm with random coordinates for each glowworm's
optimization vector.
"""

import argparse
import numpy as np
from lightdock.util.parser import SetupCommandLineParser
from lightdock.prep.simulation import read_input_structure, save_lightdock_structure, \
                                      calculate_starting_positions, prepare_results_environment, \
                                      create_setup_file, calculate_anm, parse_restraints_file, \
                                      get_restraints
from lightdock.constants import DEFAULT_LIGHTDOCK_PREFIX, DEFAULT_ELLIPSOID_DATA_EXTENSION, \
                                DEFAULT_NMODES_REC, DEFAULT_REC_NM_FILE, DEFAULT_NMODES_LIG, DEFAULT_LIG_NM_FILE
from lightdock.mathutil.ellipsoid import MinimumVolumeEllipsoid
from lightdock.util.logger import LoggingManager
from lightdock.error.lightdock_errors import LightDockError


log = LoggingManager.get_logger('lightdock_setup')


if __name__ == "__main__":

    try:
        parser = SetupCommandLineParser()
        args = parser.args

        # Read input structures
        receptor = read_input_structure(args.receptor_pdb, args.noxt, args.noh, args.verbose_parser)
        ligand = read_input_structure(args.ligand_pdb, args.noxt, args.noh, args.verbose_parser)
        
        # Move structures to origin
        rec_translation = receptor.move_to_origin()
        lig_translation = ligand.move_to_origin()

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

        # Calculate and save ANM if required
        if args.use_anm:
            if args.anm_rec > 0:
                log.info("Calculating ANM for receptor molecule...")
                calculate_anm(receptor, args.anm_rec, DEFAULT_REC_NM_FILE)
            if args.anm_lig > 0:
                log.info("Calculating ANM for ligand molecule...")
                calculate_anm(ligand, args.anm_lig, DEFAULT_LIG_NM_FILE)

        # Parse restraints if any:
        receptor_restraints = ligand_restraints = None
        if args.restraints:
            log.info("Reading restraints from %s" % args.restraints)
            restraints = parse_restraints_file(args.restraints)

            # Calculate number of restraints in order to check them
            num_rec_active = len(restraints['receptor']['active'])
            num_rec_passive = len(restraints['receptor']['passive'])
            num_lig_active = len(restraints['ligand']['active'])
            num_lig_passive = len(restraints['ligand']['passive'])

            # Complain if not a single restraint has been defined, but restraints are enabled
            if not num_rec_active and not num_rec_passive and not num_lig_active and not num_lig_passive:
                raise LightDockError("Restraints file specified, but not a single restraint found")

            # Check if restraints correspond with real residues
            receptor_restraints = get_restraints(receptor, restraints['receptor'])
            args.receptor_restraints = restraints['receptor']
            ligand_restraints = get_restraints(ligand, restraints['ligand'])
            args.ligand_restraints = restraints['ligand']

            log.info("Number of receptor restraints is: %d (active), %d (passive)" % (num_rec_active, num_rec_passive))
            log.info("Number of ligand restraints is: %d (active), %d (passive)" % (num_lig_active, num_lig_passive))

        rec_restraints = None
        try:
            rec_restraints = receptor_restraints['active'] + receptor_restraints['passive']
        except:
            pass

        lig_restraints = None
        try:
            lig_restraints = ligand_restraints['active'] + ligand_restraints['passive']
        except:
            pass

        # Calculate surface points (swarm centers) over receptor structure
        starting_points_files = calculate_starting_positions(receptor, ligand, 
                                                             args.swarms, args.glowworms, 
                                                             args.starting_points_seed,
                                                             rec_restraints, lig_restraints,
                                                             rec_translation, lig_translation,
                                                             args.ftdock_file, args.use_anm, args.anm_seed,
                                                             args.anm_rec, args.anm_lig, args.membrane)
        if len(starting_points_files) != args.swarms:
            args.swarms = len(starting_points_files)
            log.info('Number of swarms is %d after applying restraints' % args.swarms)

        # Create simulation folders
        prepare_results_environment(args.swarms)

        # Dump to a setup file the actual configuration
        create_setup_file(args)

        log.info("LightDock setup OK")

    except LightDockError as error:
        log.error("LightDock setup failed. Please see:")
        log.error(error)
