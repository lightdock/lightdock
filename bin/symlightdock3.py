#!/usr/bin/env python3

"""SymLightDock setup step"""

import traceback
from pathlib import Path
from lightdock.util.logger import LoggingManager
from lightdock.error.lightdock_errors import NotSupportedInScoringError
from lightdock.parallel.kraken import Kraken
from lightdock.parallel.util import GSOClusterTask
from lightdock.mathutil.lrandom import MTGenerator
from lightdock.gso.parameters import GSOParameters
from symlightdock.util.parser import CommandLineParser
from symlightdock.prep.simulation import get_setup_from_file, create_simulation_info_file, read_input_structure, \
                                         load_starting_positions
from symlightdock.constants import DEFAULT_SYMLIGHTDOCK_PREFIX, DEFAULT_NM_FILE, NUMPY_FILE_SAVE_EXTENSION
from symlightdock.scoring.fastdfire.driver import DFIRE, DFIREAdapter

log = LoggingManager.get_logger('symlightdock')


def set_gso(number_of_glowworms, adapters, scoring_functions, initial_positions, seed,
            step_translation, step_rotation, configuration_file=None,
            use_anm=False, nmodes_step=0.1, anm_rec=DEFAULT_NMODES_REC, anm_lig=DEFAULT_NMODES_LIG,
            local_minimization=False):
    """Creates a symlightdock GSO simulation object"""
    random_number_generator = MTGenerator(seed)
    if configuration_file:
        gso_parameters = GSOParameters(configuration_file)
    else:
        gso_parameters = GSOParameters()
    builder = SymLightdockGSOBuilder()
    if not use_anm:
        num_anm = 0
    gso = builder.create_from_file(number_of_glowworms, random_number_generator, gso_parameters,
                                   adapters, scoring_functions, initial_positions,
                                   step_translation, step_rotation, nmodes_step, num_anm)
    return gso


def prepare_gso_tasks(parser, adapters, scoring_functions, starting_points_files):
    """Creates the parallel GSOTasks objects to be executed by the scheduler"""
    tasks = []
    # Prepare tasks depending on swarms to simulate
    if parser.args.trials_list:
        trial_ids = parser.args.trials_list
        if min(trial_ids) < 0 or max(trial_ids) >= parser.args.trials:
            raise Exception("Wrong list of trials")
    else:
        trial_ids = list(range(parser.args.trials))

    for trial_id in trial_ids:
        gso = set_gso(parser.args.glowworms, adapters, scoring_functions, starting_points_files[trial_id],
                      parser.args.gso_seed, parser.args.translation_step,
                      parser.args.rotation_step, parser.args.configuration_file,
                      parser.args.use_anm, parser.args.nmodes_step, parser.args.num_anm)
        saving_path = "%s%d" % (DEFAULT_TRIAL_FOLDER, trial_id)
        task = GSOClusterTask(trial_id, gso, parser.args.steps, saving_path)
        tasks.append(task)
    return tasks


if __name__ == "__main__":

    try:
        parser = CommandLineParser()
        args = parser.args

        # Read setup and add it to the actual args object
        setup = get_setup_from_file(args.setup_file)
        for k, v in setup.items():
            setattr(args, k, v)

        info_file = create_simulation_info_file(args)
        log.info("simulation parameters saved to %s" % info_file)

        # Read input structures (use parsed ones)
        parsed_lightdock = Path(args.structure_pdb).absolute().parent / \
                                Path(DEFAULT_SYMLIGHTDOCK_PREFIX % Path(args.structure_pdb).name)
        structure = read_input_structure(parsed_lightdock, args.noxt, args.noh, args.verbose_parser)

        structure.move_to_origin()

        if args.use_anm:
            try:
                structure.n_modes = read_nmodes(f"{DEFAULT_NM_FILE}{NUMPY_FILE_SAVE_EXTENSION}")
            except:
                log.warning("No ANM found for structure")
                structure.n_modes = None

        starting_points_files = load_starting_positions(args.trials, args.glowworms, args.use_anm, args.num_anm)

        log.info("Loading scoring function...")
        restraints = None # TODO
        adapter = DFIREAdapter(structure, restraints)
        scoring_function = DFIRE()

        tasks = prepare_gso_tasks(parser, [adapter], [scoring_function], starting_points_files)

        # Preparing the parallel execution
        kraken = Kraken(tasks, parser.args.cores, parser.args.profiling)
        log.info("Monster spotted")
        reports_queue = kraken.release()
        log.info("Finished.")

    except NotSupportedInScoringError as score_error:
        log.error("Error found in selected scoring function:")
        log.error(score_error)

    except KeyboardInterrupt:
        log.info("Caught interrupt...")
        try:
            kraken.sink()
        except:
            pass
        log.info("bye.")

    except OSError as e:
        log.error("OS error found")
        try:
            kraken.sink()
        except:
            pass
        raise e
