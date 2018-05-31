"""LightDock simulation using the multiprocessing library for parallelization"""

import os
import importlib
import glob

from lightdock.util.logger import LoggingManager
from lightdock.util.parser import CommandLineParser
from lightdock.prep.simulation import get_setup_from_file, create_simulation_info_file, read_input_structure, \
                                      load_starting_positions
from lightdock.gso.boundaries import Boundary, BoundingBox
from lightdock.gso.algorithm import LightdockGSOBuilder
from lightdock.mathutil.lrandom import MTGenerator
from lightdock.gso.parameters import GSOParameters
from lightdock.constants import MAX_TRANSLATION, MAX_ROTATION, DEFAULT_SCORING_FUNCTION, DEFAULT_SWARM_FOLDER
from lightdock.parallel.kraken import Kraken
from lightdock.parallel.util import GSOClusterTask
from lightdock.scoring.multiple import ScoringConfiguration


log = LoggingManager.get_logger('lightdock')


def set_gso(number_of_glowworms, adapters, scoring_functions, initial_positions, seed,
            step_translation, step_rotation, configuration_file=None, nm=False, nmodes_step=0.1,
            local_minimization=False):
    """Creates a lightdock GSO simulation object"""
    # Only dimension is relevant, initial positions are not randomized, but generated
    boundaries = [Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                  Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                  Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                  Boundary(-MAX_ROTATION, MAX_ROTATION),
                  Boundary(-MAX_ROTATION, MAX_ROTATION),
                  Boundary(-MAX_ROTATION, MAX_ROTATION),
                  Boundary(-MAX_ROTATION, MAX_ROTATION)]
    if nm:
        boundaries.extend([Boundary(MIN_EXTENT, MAX_EXTENT) for _ in xrange(DEFAULT_NMODES_REC)])
        boundaries.extend([Boundary(MIN_EXTENT, MAX_EXTENT) for _ in xrange(DEFAULT_NMODES_LIG)])

    bounding_box = BoundingBox(boundaries)

    random_number_generator = MTGenerator(seed)
    if configuration_file:
        gso_parameters = GSOParameters(configuration_file)
    else:
        gso_parameters = GSOParameters()
    builder = LightdockGSOBuilder()
    gso = builder.create_from_file(number_of_glowworms, random_number_generator, gso_parameters,
                                   adapters, scoring_functions, bounding_box, initial_positions,
                                   step_translation, step_rotation, nmodes_step, local_minimization)
    return gso



def set_scoring_function(parser, receptor, ligand):
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

    for scoring_function, weight in functions.iteritems():
        log.info("Loading scoring function...")
        scoring_function_module = "lightdock.scoring.%s.driver" % scoring_function
        module = importlib.import_module(scoring_function_module)

        log.info("Using %s scoring function" % module.DefinedScoringFunction.__name__)

        CurrentScoringFunction = getattr(module, "DefinedScoringFunction")
        CurrentModelAdapter = getattr(module, "DefinedModelAdapter")
        adapter = CurrentModelAdapter(receptor, ligand)
        scoring_function = CurrentScoringFunction(weight)
        adapters.append(adapter)
        scoring_functions.append(scoring_function)
        log.info("Done.")
    return scoring_functions, adapters


def prepare_gso_tasks(parser, adapters, scoring_functions, starting_points_files):
    """Creates the parallel GSOTasks objects to be executed by the scheduler"""
    tasks = []
    for id_swarm in range(parser.args.swarms):
        gso = set_gso(parser.args.glowworms, adapters, scoring_functions, starting_points_files[id_swarm],
                      parser.args.gso_seed, parser.args.translation_step,
                      parser.args.rotation_step, parser.args.configuration_file,
                      parser.args.use_anm, parser.args.nmodes_step,
                      parser.args.local_minimization)
        saving_path = "%s%d" % (DEFAULT_SWARM_FOLDER, id_swarm)
        task = GSOClusterTask(id_swarm, gso, parser.args.steps, saving_path)
        tasks.append(task)
    return tasks


def set_normal_modes(receptor, ligand):
    """Only calculates normal modes for representative structure"""
    from lightdock.structure.nm import calculate_nmodes, write_nmodes
    modes = calculate_nmodes(receptor.structure_file_names[receptor.representative_id],
                             DEFAULT_NMODES_REC, receptor)
    receptor.n_modes = modes
    write_nmodes(modes, DEFAULT_REC_NM_FILE)
    log.info("%d normal modes calculated for receptor" % DEFAULT_NMODES_REC)

    modes = calculate_nmodes(ligand.structure_file_names[ligand.representative_id],
                             DEFAULT_NMODES_LIG, ligand)
    ligand.n_modes = modes
    write_nmodes(modes, DEFAULT_LIG_NM_FILE)
    log.info("%d normal modes calculated for ligand" % DEFAULT_NMODES_LIG)


def run_simulation(parser):
    """Main program"""
    try:
        parser = CommandLineParser()
        args = parser.args

        # Read setup and add it to the actual args object
        setup = get_setup_from_file(args.setup_file)
        for k, v in setup.iteritems():
            setattr(args, k, v)

        info_file = create_simulation_info_file(args)
        log.info("simulation parameters saved to %s" % info_file)

        # Read input structures
        receptor = read_input_structure(args.receptor_pdb, args.noxt)
        ligand = read_input_structure(args.ligand_pdb, args.noxt)

        if args.use_anm:
            set_normal_modes(receptor, ligand)

        starting_points_files = load_starting_positions(args.swarms, args.glowworms, args.use_anm)

        scoring_functions, adapters = set_scoring_function(parser, receptor, ligand)

        tasks = prepare_gso_tasks(parser, adapters, scoring_functions, starting_points_files)

        # Preparing the parallel execution
        kraken = Kraken(tasks, parser.args.cores, parser.args.profiling)
        log.info("Monster spotted")
        reports_queue = kraken.release()
        log.info("Finished.")

    except KeyboardInterrupt:
        log.info("Caught interrupt...")
        try:
            kraken.sink()
        except:
            pass
        log.info("bye.")
