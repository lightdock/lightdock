#!/usr/bin/env python3

"""Before launching the LightDock simulation, a setup step is required.

This step parses the input PDB structures, calculates the minimum ellipsoid
containing each of them, calculates the swarms on the surface of the
receptor and populates each swarm with random coordinates for each glowworm's
optimization vector.
"""

from pathlib import Path
import numpy as np
from lightdock.util.parser import SetupCommandLineParser
from lightdock.prep.simulation import (
    read_input_structure,
    save_lightdock_structure,
    calculate_starting_positions,
    prepare_results_environment,
    create_setup_file,
    calculate_anm,
    parse_restraints_file,
    get_restraints,
)
from lightdock.constants import (
    DEFAULT_LIGHTDOCK_PREFIX,
    DEFAULT_ELLIPSOID_DATA_EXTENSION,
    DEFAULT_REC_NM_FILE,
    DEFAULT_LIG_NM_FILE,
)
from lightdock.mathutil.ellipsoid import MinimumVolumeEllipsoid
from lightdock.util.logger import LoggingManager
from lightdock.error.lightdock_errors import LightDockError
from lightdock.version import CURRENT_VERSION


log = LoggingManager.get_logger("lightdock3_setup")


if __name__ == "__main__":

    try:
        parser = SetupCommandLineParser()
        args = parser.args

        # Read input structures
        receptor = read_input_structure(
            args.receptor_pdb, args.noxt, args.noh, args.now, args.verbose_parser
        )
        ligand = read_input_structure(
            args.ligand_pdb, args.noxt, args.noh, args.now, args.verbose_parser
        )

        # Move structures to origin
        rec_translation = receptor.move_to_origin()
        lig_translation = ligand.move_to_origin()

        # Calculate reference points for receptor
        log.info(f"Calculating reference points for receptor {args.receptor_pdb}...")
        ellipsoid_data_file = "%s%s" % (
            DEFAULT_LIGHTDOCK_PREFIX % receptor.structure_file_names[0],
            DEFAULT_ELLIPSOID_DATA_EXTENSION,
        )
        if not Path(ellipsoid_data_file).exists():
            log.info("Reference points for receptor found, skipping")
        else:
            rec_ellipsoid = MinimumVolumeEllipsoid(
                receptor.representative().coordinates
            )
            np.save(ellipsoid_data_file, np.array([rec_ellipsoid.center.copy()]))
        log.info("Done.")

        # Calculate reference points for ligand
        log.info("Calculating reference points for ligand %s..." % args.ligand_pdb)
        ellipsoid_data_file = "%s%s" % (
            DEFAULT_LIGHTDOCK_PREFIX % ligand.structure_file_names[0],
            DEFAULT_ELLIPSOID_DATA_EXTENSION,
        )
        if not Path(ellipsoid_data_file).exists():
            log.info("Reference points for ligand found, skipping")
        else:
            lig_ellipsoid = MinimumVolumeEllipsoid(ligand.representative().coordinates)
            np.save(ellipsoid_data_file, np.array([lig_ellipsoid.center.copy()]))
        log.info("Done.")

        # Save to file parsed structures
        save_lightdock_structure(receptor)
        save_lightdock_structure(ligand)

        # Calculate and save ANM if required
        if args.use_anm:
            if args.anm_rec > 0:
                log.info("Calculating ANM for receptor molecule...")
                calculate_anm(receptor, args.anm_rec, args.anm_rec_rmsd, args.anm_seed, DEFAULT_REC_NM_FILE)
            if args.anm_lig > 0:
                log.info("Calculating ANM for ligand molecule...")
                calculate_anm(ligand, args.anm_lig, args.anm_lig_rmsd, args.anm_seed, DEFAULT_LIG_NM_FILE)

        # Parse restraints if any:
        receptor_restraints = ligand_restraints = None
        if args.restraints:
            log.info(f"Reading restraints from {args.restraints}")
            restraints = parse_restraints_file(args.restraints)

            # Calculate number of restraints in order to check them
            num_rec_active = len(restraints["receptor"]["active"])
            num_rec_passive = len(restraints["receptor"]["passive"])
            num_rec_blocked = len(restraints["receptor"]["blocked"])
            num_lig_active = len(restraints["ligand"]["active"])
            num_lig_passive = len(restraints["ligand"]["passive"])

            # Complain if not a single restraint has been defined, but restraints are enabled
            if (
                not num_rec_active
                and not num_rec_passive
                and not num_rec_blocked
                and not num_lig_active
                and not num_lig_passive
            ):
                raise LightDockError(
                    "Restraints file specified, but not a single restraint found"
                )

            # Check if restraints correspond with real residues
            receptor_restraints = get_restraints(receptor, restraints["receptor"])
            args.receptor_restraints = restraints["receptor"]
            ligand_restraints = get_restraints(ligand, restraints["ligand"])
            args.ligand_restraints = restraints["ligand"]

            log.info(
                f"Number of receptor restraints is: {num_rec_active} (active), {num_rec_passive} (passive)"
            )
            log.info(
                f"Number of ligand restraints is: {num_lig_active} (active), {num_lig_passive} (passive)"
            )

        try:
            lig_restraints = ligand_restraints["active"] + ligand_restraints["passive"]
        except (KeyError, TypeError):
            lig_restraints = None

        # Calculate surface points (swarm centers) over receptor structure
        starting_points_files = calculate_starting_positions(
            receptor,
            ligand,
            args.swarms,
            args.glowworms,
            args.starting_points_seed,
            receptor_restraints,
            lig_restraints,
            rec_translation,
            lig_translation,
            args.surface_density,
            args.use_anm,
            args.anm_seed,
            args.anm_rec,
            args.anm_lig,
            args.membrane,
            args.transmembrane,
            args.write_starting_positions,
            args.swarm_radius,
            args.flip,
            args.fixed_distance,
            args.swarms_per_restraint,
            args.dense_sampling,
        )
        if len(starting_points_files) != args.swarms:
            args.swarms = len(starting_points_files)
            log.info(f"Number of calculated swarms is {args.swarms}")

        # Create simulation folders
        prepare_results_environment(args.swarms)

        # Add manually setup version
        args.setup_version = CURRENT_VERSION

        # Dump to a setup file the actual configuration
        create_setup_file(args)

        log.info("LightDock setup OK")

    except LightDockError as error:
        log.error("LightDock setup failed. Please see:")
        log.error(error)
