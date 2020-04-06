"""LightDock simulation using the MPI4py library for parallelization"""

import os
import importlib
import glob
from mpi4py import MPI

from lightdock.util.logger import LoggingManager
from lightdock.util.parser import CommandLineParser
from lightdock.prep.simulation import get_setup_from_file, create_simulation_info_file, read_input_structure, \
                                      load_starting_positions, get_default_box
from lightdock.gso.algorithm import LightdockGSOBuilder
from lightdock.mathutil.lrandom import MTGenerator
from lightdock.gso.parameters import GSOParameters
from lightdock.constants import DEFAULT_SCORING_FUNCTION, DEFAULT_SWARM_FOLDER, \
                                DEFAULT_REC_NM_FILE, DEFAULT_LIG_NM_FILE, NUMPY_FILE_SAVE_EXTENSION, \
                                DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG
from lightdock.parallel.kraken import Kraken
from lightdock.parallel.util import GSOClusterTask
from lightdock.scoring.multiple import ScoringConfiguration
from lightdock.structure.nm import read_nmodes
from lightdock.error.lightdock_errors import NotSupportedInScoringError


log = LoggingManager.get_logger('lightdock')


def set_gso(number_of_glowworms, adapters, scoring_functions, initial_positions, seed,
            step_translation, step_rotation, configuration_file=None, 
            use_anm=False, nmodes_step=0.1, anm_rec=DEFAULT_NMODES_REC, anm_lig=DEFAULT_NMODES_LIG,
            local_minimization=False):
    """Creates a lightdock GSO simulation object"""

    bounding_box = get_default_box(use_anm, anm_rec, anm_lig)

    random_number_generator = MTGenerator(seed)
    if configuration_file:
        gso_parameters = GSOParameters(configuration_file)
    else:
        gso_parameters = GSOParameters()
    builder = LightdockGSOBuilder()
    if not use_anm:
        anm_rec = anm_lig = 0
    gso = builder.create_from_file(number_of_glowworms, random_number_generator, gso_parameters,
                                   adapters, scoring_functions, bounding_box, initial_positions,
                                   step_translation, step_rotation, nmodes_step, local_minimization)
    return gso


def set_scoring_function(parser, receptor, ligand, minion_id):
    """Set scoring function and docking models"""
    scoring_functions = []
    adapters = []
    if parser.args.scoring_function and os.path.exists(parser.args.scoring_function):
        # Multiple scoring functions found
        functions = ScoringConfiguration.parse_file(parser.args.scoring_function)
    else:
        if parser.args.scoring_function:
            functions = {parser.args.scoring_function: '1.0'}
        else:
            functions = {DEFAULT_SCORING_FUNCTION: '1.0'}

    for scoring_function, weight in functions.items():
        log.info("[Minion %d] Loading scoring function..." % minion_id)
        scoring_function_module = "lightdock.scoring.%s.driver" % scoring_function
        module = importlib.import_module(scoring_function_module)

        log.info("[Minion %d] Using %s scoring function" % (minion_id, module.DefinedScoringFunction.__name__))

        CurrentScoringFunction = getattr(module, "DefinedScoringFunction")
        CurrentModelAdapter = getattr(module, "DefinedModelAdapter")

        receptor_restraints = ligand_restraints = None
        try:
            receptor_restraints = parser.args.receptor_restraints['active']
        except:
            pass
        try:
            ligand_restraints = parser.args.ligand_restraints['active']
        except:
            pass
        adapter = CurrentModelAdapter(receptor, ligand, receptor_restraints, ligand_restraints)
        scoring_function = CurrentScoringFunction(weight)
        adapters.append(adapter)
        scoring_functions.append(scoring_function)
        log.info("[Minion %d] Done." % minion_id)
    return scoring_functions, adapters


def run_simulation(parser):
    """Main program, includes MPI directives"""
    try:
        comm = MPI.COMM_WORLD

        parser = CommandLineParser()
        args = parser.args

        # Read setup and add it to the actual args object
        setup = get_setup_from_file(args.setup_file)
        for k, v in setup.items():
            setattr(args, k, v)

        minion_id = comm.rank
        if minion_id == 0:
            info_file = create_simulation_info_file(args)
            log.info("simulation parameters saved to %s" % info_file)
        comm.Barrier()

        # Read input structures (use parsed ones)
        parsed_lightdock_receptor = os.path.join(os.path.dirname(args.receptor_pdb),
                                       DEFAULT_LIGHTDOCK_PREFIX % os.path.basename(args.receptor_pdb))
        receptor = read_input_structure(parsed_lightdock_receptor, args.noxt, args.noh, args.verbose_parser)
        parsed_lightdock_ligand = os.path.join(os.path.dirname(args.ligand_pdb),
                                       DEFAULT_LIGHTDOCK_PREFIX % os.path.basename(args.ligand_pdb))
        ligand = read_input_structure(parsed_lightdock_ligand, args.noxt, args.noh, args.verbose_parser)

        # CRITICAL to not break compatibility with previous results
        receptor.move_to_origin()
        ligand.move_to_origin()

        if args.use_anm:
            try:
                receptor.n_modes = read_nmodes("%s%s" % (DEFAULT_REC_NM_FILE, NUMPY_FILE_SAVE_EXTENSION) )
            except:
                log.warning("No ANM found for receptor molecule")
                receptor.n_modes = None
            try:
                ligand.n_modes = read_nmodes("%s%s" % (DEFAULT_LIG_NM_FILE, NUMPY_FILE_SAVE_EXTENSION) )
            except:
                log.warning("No ANM found for ligand molecule")
                ligand.n_modes = None

        starting_points_files = load_starting_positions(args.swarms, args.glowworms, args.use_anm,
                                                        args.anm_rec, args.anm_lig)

        comm.Barrier()

        num_workers = comm.size
        for worker_id in range(num_workers):
            if worker_id == minion_id:
                starting_points_files = glob.glob('init/initial_positions*.dat')
                scoring_functions, adapters = set_scoring_function(parser, receptor, ligand, minion_id)
                for id_swarm in range(parser.args.swarms):
                    if worker_id == (id_swarm % num_workers):
                        print('GSO cluster %d - Minion %d' % (id_swarm, minion_id))
                        gso = set_gso(parser.args.glowworms, adapters, scoring_functions,
                                      starting_points_files[id_swarm],
                                      parser.args.gso_seed, parser.args.translation_step,
                                      parser.args.rotation_step, parser.args.configuration_file,
                                      parser.args.use_anm, parser.args.nmodes_step,
                                      parser.args.anm_rec, parser.args.anm_lig,
                                      parser.args.local_minimization)
                        saving_path = "%s%d" % (DEFAULT_SWARM_FOLDER, id_swarm)
                        task = GSOClusterTask(id_swarm, gso, parser.args.steps, saving_path)
                        task.run()
        comm.Barrier()

    except NotSupportedInScoringError as score_error:
        log.error("Error found in selected scoring function:")
        log.error(score_error)

    except KeyboardInterrupt:
        log.info("Caught interrupt...")
        log.info("bye.")
