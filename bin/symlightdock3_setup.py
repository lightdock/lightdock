#!/usr/bin/env python3

"""SymLightDock setup step"""

from pathlib import Path
import numpy as np
from symlightdock.util.parser import SetupCommandLineParser
from symlightdock.prep.simulation import read_input_structure, save_lightdock_structure, \
                                         calculate_starting_positions, prepare_results_environment, \
                                         create_setup_file, calculate_anm, parse_restraints_file, get_restraints
from symlightdock.constants import DEFAULT_SYMLIGHTDOCK_PREFIX, DEFAULT_ELLIPSOID_DATA_EXTENSION, DEFAULT_NM_FILE
from lightdock.mathutil.ellipsoid import MinimumVolumeEllipsoid
from lightdock.util.logger import LoggingManager
from lightdock.error.lightdock_errors import LightDockError


log = LoggingManager.get_logger('symlightdock_setup')


if __name__ == "__main__":

    try:
        parser = SetupCommandLineParser()
        args = parser.args

        # Read input structures
        structure = read_input_structure(args.structure_pdb, args.noxt, args.noh, args.verbose_parser)

        # Move structure to origin
        translation = structure.move_to_origin()

        # Calculate reference points for receptor
        log.info(f"Calculating reference points for structure {args.structure_pdb}...")
        ellipsoid_data_file = "%s%s" % (DEFAULT_SYMLIGHTDOCK_PREFIX % structure.structure_file_names[0],
                                        DEFAULT_ELLIPSOID_DATA_EXTENSION)
        if not Path(ellipsoid_data_file).exists():
            log.info("Reference points for structure found, skipping")
        else:
            ellipsoid = MinimumVolumeEllipsoid(structure.representative().coordinates)
            np.save(ellipsoid_data_file, np.array([ellipsoid.center.copy()]))
        log.info("Done.")

        # Save to file parsed structure
        save_lightdock_structure(structure)

        # Calculate and save ANM if required
        if args.use_anm:
            if args.num_anm > 0:
                log.info("Calculating ANM for molecule...")
                calculate_anm(structure, args.num_anm, DEFAULT_NM_FILE)

        # Parse restraints if any:
        restraints = None
        if args.restraints:
            log.info(f"Reading restraints from {args.restraints}")
            restraints = parse_restraints_file(args.restraints)

            # Calculate number of restraints in order to check them
            num_active = len(restraints['receptor']['active'])
            num_passive = len(restraints['receptor']['passive'])
            num_blocked = len(restraints['receptor']['blocked'])

            # Complain if not a single restraint has been defined, but restraints are enabled
            if not num_active and not num_passive and not num_blocked:
                raise LightDockError("Restraints file specified, but not a single restraint found")

            # Check if restraints correspond with real residues
            restraints = get_restraints(structure, restraints['receptor'])
            args.restraints = restraints['receptor']

            log.info(f"Number of restraints is: {num_active} (active), {num_passive} (passive)")

        # Calculate starting conditions
        starting_points_files = calculate_starting_positions(structure,
                                                             args.trials, args.glowworms,
                                                             args.starting_points_seed,
                                                             restraints, translation,
                                                             args.use_anm, args.anm_seed, args.num_anm,
                                                             args.write_starting_positions)

        # Create simulation folders
        prepare_results_environment(args.trials)

        # Dump to a setup file the actual configuration
        create_setup_file(args)

        log.info("SymLightDock setup OK")

    except LightDockError as error:
        log.error("SymLightDock setup failed. Please see:")
        log.error(error)
