"""LightDock simulation using the MPI4py library for parallelization"""

import os
import importlib
import glob
from mpi4py import MPI

from lightdock.util.logger import LoggingManager
from lightdock.prep.poses import calculate_initial_poses
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.gso.boundaries import Boundary, BoundingBox
from lightdock.gso.algorithm import LightdockGSOBuilder
from lightdock.mathutil.lrandom import MTGenerator
from lightdock.gso.parameters import GSOParameters
from lightdock.constants import MAX_TRANSLATION, MAX_ROTATION, DEFAULT_SCORING_FUNCTION,\
    DEFAULT_POSITIONS_FOLDER, DEFAULT_CLUSTER_FOLDER, DEFAULT_LIST_EXTENSION, DEFAULT_LIGHTDOCK_PREFIX, \
    DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG, DEFAULT_REC_NM_FILE, DEFAULT_LIG_NM_FILE, MIN_EXTENT, MAX_EXTENT
from lightdock.parallel.util import GSOClusterTask
from lightdock.error.lightdock_errors import LightDockError
from lightdock.util.simulation_info import show_parameters, create_simulation_info_file
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


def get_pdb_files(input_file):
    """Get a list of the PDB files in the input_file"""
    file_names = []
    with open(input_file) as input:
        for line in input:
            file_name = line.rstrip(os.linesep)
            if os.path.exists(file_name):
                file_names.append(file_name)
    return file_names


def read_input_structures(parser, minion_id):
    """Reads the input structures.

    It is able to read PDB structures or list of PDB structures in case conformers support is needed.
    """
    atoms_to_ignore = []
    if parser.args.noxt:
        atoms_to_ignore.append('OXT')

    receptor_structures = []
    file_names = []
    file_name, file_extension = os.path.splitext(parser.args.info_receptor_pdb)
    if file_extension == DEFAULT_LIST_EXTENSION:
        file_names.extend(get_pdb_files(parser.args.info_receptor_pdb))
    else:
        file_names.append(parser.args.info_receptor_pdb)
    for file_name in file_names:
        log.info("[Minion %d] Reading %s receptor PDB file..." % (minion_id, file_name))
        atoms, residues, chains = parse_complex_from_file(file_name, atoms_to_ignore)
        receptor_structures.append({'atoms': atoms, 'residues': residues, 'chains': chains, 'file_name': file_name})
        log.info("[Minion %d] %s atoms, %s residues read." % (minion_id, len(atoms), len(residues)))

    ligand_structures = []
    file_names = []
    file_name, file_extension = os.path.splitext(parser.args.info_ligand_pdb)
    if file_extension == DEFAULT_LIST_EXTENSION:
        file_names.extend(get_pdb_files(parser.args.info_ligand_pdb))
    else:
        file_names.append(parser.args.info_ligand_pdb)
    for file_name in file_names:
        log.info("[Minion %d] Reading %s ligand PDB file..." % (minion_id, file_name))
        atoms, residues, chains = parse_complex_from_file(file_name, atoms_to_ignore)
        ligand_structures.append({'atoms': atoms, 'residues': residues, 'chains': chains, 'file_name': file_name})
        log.info("[Minion %d] %s atoms, %s residues read." % (minion_id, len(atoms), len(residues)))

    # Representatives are now the first structure, but this could change in the future
    receptor = Complex.from_structures(receptor_structures)
    ligand = Complex.from_structures(ligand_structures)
    return receptor, ligand


def save_lightdock_structures(receptor, ligand):
    """Saves the parsed PDB structures"""
    log.info("Saving processed structures to PDB files...")
    for structure_index, file_name in enumerate(receptor.structure_file_names):
        moved_file_name = os.path.join(os.path.dirname(file_name),
                                       DEFAULT_LIGHTDOCK_PREFIX % os.path.basename(file_name))
        write_pdb_to_file(receptor, moved_file_name, receptor[structure_index])

    for structure_index, file_name in enumerate(ligand.structure_file_names):
        moved_file_name = os.path.join(os.path.dirname(file_name),
                                       DEFAULT_LIGHTDOCK_PREFIX % os.path.basename(file_name))
        write_pdb_to_file(ligand, moved_file_name, ligand[structure_index])
    log.info("Done.")


def calculate_starting_positions(parser, receptor, ligand):
    """Calculates the starting positions for any of the glowworm agents"""
    log.info("Calculating starting positions...")
    init_folder = DEFAULT_POSITIONS_FOLDER
    if not os.path.isdir(init_folder):
        os.mkdir(init_folder)
        starting_points_files = calculate_initial_poses(receptor, ligand,
                                                        parser.args.clusters, parser.args.glowworms,
                                                        parser.args.starting_points_seed, init_folder,
                                                        parser.args.ftdock_file, parser.args.nm,
                                                        parser.args.nm_seed, DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG)
        log.info("Generated %d positions files" % len(starting_points_files))
    else:
        log.warning("Folder %s already exists, skipping calculation" % init_folder)
        starting_points_files = glob.glob('init/initial_positions*.dat')
        if len(starting_points_files) != parser.args.clusters:
            raise LightDockError("The number of initial positions files does not correspond with the number of cluster")
    log.info("Done.")
    return starting_points_files


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

    for scoring_function, weight in functions.iteritems():
        log.info("[Minion %d] Loading scoring function..." % minion_id)
        scoring_function_module = "lightdock.scoring.%s.driver" % scoring_function
        module = importlib.import_module(scoring_function_module)

        log.info("[Minion %d] Using %s scoring function" % (minion_id, module.DefinedScoringFunction.__name__))

        CurrentScoringFunction = getattr(module, "DefinedScoringFunction")
        CurrentModelAdapter = getattr(module, "DefinedModelAdapter")
        adapter = CurrentModelAdapter(receptor, ligand)
        scoring_function = CurrentScoringFunction(weight)
        adapters.append(adapter)
        scoring_functions.append(scoring_function)
        log.info("[Minion %d] Done." % minion_id)
    return scoring_functions, adapters


def prepare_results_environment(parser):
    """Prepares the folder structure required by the simulation"""
    log.info("Preparing environment")
    for id_cluster in range(parser.args.clusters):
        saving_path = "%s%d" % (DEFAULT_CLUSTER_FOLDER, id_cluster)
        if os.path.isdir(saving_path):
            raise LightDockError("Results folder %s already exists" % saving_path)
        else:
            os.mkdir(saving_path)
    log.info("Done.")


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
    """Main program, includes MPI directives"""
    try:
        comm = MPI.COMM_WORLD

        minion_id = comm.rank
        if minion_id == 0:
            show_parameters(log, parser)
            info_file = create_simulation_info_file(parser)
            log.info("lightdock parameters saved to %s" % info_file)
        comm.Barrier()

        receptor, ligand = read_input_structures(parser, minion_id)
        # Move structures to origin
        receptor.move_to_origin()
        ligand.move_to_origin()

        if parser.args.nm:
            set_normal_modes(receptor, ligand)

        if minion_id == 0:
            save_lightdock_structures(receptor, ligand)
            calculate_starting_positions(parser, receptor, ligand)
            prepare_results_environment(parser)
        comm.Barrier()

        num_workers = comm.size
        for worker_id in xrange(num_workers):
            if worker_id == minion_id:
                starting_points_files = glob.glob('init/initial_positions*.dat')
                scoring_functions, adapters = set_scoring_function(parser, receptor, ligand, minion_id)
                for id_cluster in xrange(parser.args.clusters):
                    if worker_id == (id_cluster % num_workers):
                        print 'GSO cluster %d - Minion %d' % (id_cluster, minion_id)
                        gso = set_gso(parser.args.glowworms, adapters, scoring_functions,
                                      starting_points_files[id_cluster],
                                      parser.args.gso_seed, parser.args.translation_step,
                                      parser.args.rotation_step, parser.args.configuration_file,
                                      parser.args.nm, parser.args.nmodes_step,
                                      parser.args.local_minimization)
                        saving_path = "%s%d" % (DEFAULT_CLUSTER_FOLDER, id_cluster)
                        task = GSOClusterTask(id_cluster, gso, parser.args.steps, saving_path)
                        task.run()
        comm.Barrier()

    except KeyboardInterrupt:
        log.info("Caught interrupt...")
        log.info("bye.")
