import os
import glob
import json
import time
from lightdock.prep.poses import calculate_initial_poses
from lightdock.constants import DEFAULT_POSITIONS_FOLDER, DEFAULT_SWARM_FOLDER, DEFAULT_LIST_EXTENSION, \
    DEFAULT_LIGHTDOCK_PREFIX, DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG, DEFAULT_REC_NM_FILE, DEFAULT_LIG_NM_FILE, \
    MIN_EXTENT, MAX_EXTENT, DEFAULT_SETUP_FILE, DEFAULT_LIGHTDOCK_INFO, DEFAULT_POSITIONS_FOLDER, DEFAULT_STARTING_PREFIX
from lightdock.util.logger import LoggingManager
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.structure.nm import calculate_nmodes, write_nmodes
from lightdock.error.lightdock_errors import LightDockError


log = LoggingManager.get_logger('lightdock_setup')


def get_pdb_files(input_file):
    """Get a list of the PDB files in the input_file"""
    file_names = []
    with open(input_file) as input:
        for line in input:
            file_name = line.rstrip(os.linesep)
            if os.path.exists(file_name):
                file_names.append(file_name)
    return file_names


def read_input_structure(pdb_file_name, ignore_oxt=True):
    """Reads the input structure.

    The arguments pdb_file_name can be a PDB file or a file 
    containing a list of PDB files.

    ignore_oxt flag avoids saving OXT atoms.
    """
    atoms_to_ignore = []
    if ignore_oxt:
        atoms_to_ignore.append('OXT')

    structures = []
    file_names = []
    file_name, file_extension = os.path.splitext(pdb_file_name)
    if file_extension == DEFAULT_LIST_EXTENSION:
        file_names.extend(get_pdb_files(pdb_file_name))
    else:
        file_names.append(pdb_file_name)
    for file_name in file_names:
        log.info("Reading structure from %s PDB file..." % file_name)
        atoms, residues, chains = parse_complex_from_file(file_name, atoms_to_ignore)
        structures.append({'atoms': atoms, 'residues': residues, 'chains': chains, 'file_name': file_name})
        log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    # Representatives are now the first structure, but this could change in the future
    structure = Complex.from_structures(structures)
    return structure


def save_lightdock_structure(structure):
    """Saves the structure parsed by LightDock"""
    log.info("Saving processed structure to PDB file...")
    for structure_index, file_name in enumerate(structure.structure_file_names):
        moved_file_name = os.path.join(os.path.dirname(file_name),
                                       DEFAULT_LIGHTDOCK_PREFIX % os.path.basename(file_name))
        write_pdb_to_file(structure, moved_file_name, structure[structure_index])
    log.info("Done.")


def calculate_anm(structure, num_nmodes, file_name):
    """Calculates ANM for representative structure"""
    modes = calculate_nmodes(structure.structure_file_names[structure.representative_id],
                             num_nmodes, structure)
    structure.n_modes = modes
    write_nmodes(modes, file_name)
    log.info("%d normal modes calculated" % num_nmodes)


def check_starting_file(file_name, glowworms, use_anm, anm_rec, anm_lig):
    """Check if the file_name file contains the required starting coordinates"""
    with open(file_name) as input_file:
        lines = [line.rstrip(os.linesep) for line in input_file.readlines()]
        num_glowworms = 0
        expected_size = 7 + anm_rec + anm_lig if use_anm else 7
        for line in lines:
            fields = line.split()
            if len(fields) != expected_size:
                return False
            num_glowworms += 1
        return num_glowworms == glowworms


def calculate_starting_positions(receptor, ligand, swarms, glowworms, starting_points_seed, 
    ftdock_file=None, use_anm=False, anm_seed=0, anm_rec=DEFAULT_NMODES_REC, anm_lig=DEFAULT_NMODES_LIG):
    """Defines the starting positions of each glowworm in the simulation.

    If the init_folder already exists, uses the starting positions from this folder.
    """
    log.info("Calculating starting positions...")
    init_folder = DEFAULT_POSITIONS_FOLDER
    if not os.path.isdir(init_folder):
        os.mkdir(init_folder)
        starting_points_files = calculate_initial_poses(receptor, ligand,
                                                        swarms, glowworms,
                                                        starting_points_seed, init_folder,
                                                        ftdock_file, use_anm,
                                                        anm_seed, anm_rec, anm_lig)
        log.info("Generated %d positions files" % len(starting_points_files))
    else:
        log.warning("Folder %s already exists, skipping calculation" % init_folder)
        pattern = os.path.join(DEFAULT_POSITIONS_FOLDER, "%s*.dat" % DEFAULT_STARTING_PREFIX)
        starting_points_files = glob.glob(pattern)
        if len(starting_points_files) != swarms:
            raise LightDockError("The number of initial positions files does not correspond with the number of swarms")
        for starting_point_file in starting_points_files:
            if not check_starting_file(starting_point_file, glowworms, use_anm, anm_rec, anm_lig):
                raise LightDockError("Error reading starting coordinates from file %s" % starting_point_file)
    log.info("Done.")
    return starting_points_files


def load_starting_positions(swarms, glowworms, use_anm, anm_rec=DEFAULT_NMODES_REC, anm_lig=DEFAULT_NMODES_LIG):
    """Gets the list of starting positions of this simulation"""
    pattern = os.path.join(DEFAULT_POSITIONS_FOLDER, "%s*.dat" % DEFAULT_STARTING_PREFIX)
    starting_points_files = glob.glob(pattern)
    if len(starting_points_files) != swarms:
        raise LightDockError("The number of initial positions files does not correspond with the number of swarms")
    for starting_point_file in starting_points_files:
        if not check_starting_file(starting_point_file, glowworms, use_anm, anm_rec, anm_lig):
            raise LightDockError("Error reading starting coordinates from file %s" % starting_point_file)
    return starting_points_files


def prepare_results_environment(swarms=10):
    """Prepares the folder structure required by the simulation"""
    log.info("Preparing environment")
    for id_swarm in range(swarms):
        saving_path = "%s%d" % (DEFAULT_SWARM_FOLDER, id_swarm)
        if os.path.isdir(saving_path):
            raise LightDockError("Simulation folder %s already exists" % saving_path)
        else:
            os.mkdir(saving_path)
    log.info("Done.")


def create_setup_file(args):
    """Dumps the args object into a setup file in JSON format"""
    with open(DEFAULT_SETUP_FILE, 'w') as fp:
        json.dump(vars(args), fp, indent=4)


def get_setup_from_file(file_name):
    """Reads the simulation setup from the file_name"""
    with open(file_name) as input_file:
        return json.load(input_file)


def create_simulation_info_file(args, path='.', file_name=DEFAULT_LIGHTDOCK_INFO):
    """Creates a simulation file from which recover from in a new simulation"""
    # Create the simulation info file. If it exists, includes a number
    # in the extension to avoid collision
    output_file_name = os.path.join(path, file_name)
    if os.path.isfile(output_file_name):
        original_file_name = output_file_name
        i = 1
        while os.path.isfile(output_file_name) and i < 255:
            output_file_name = "%s.%d" % (original_file_name, i)
            i += 1
        if i == 255:
            raise LightDockError('Too many simulation files')

    # Data to store
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    data = {'start_time': now}
    data.update(vars(args))

    # Store the data in the file sorted alphabetically
    with open(output_file_name, 'w') as fp:
        json.dump(vars(args), fp, indent=4, sort_keys=True)

    return output_file_name
