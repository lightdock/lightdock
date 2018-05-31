import os
import glob
from lightdock.prep.poses import calculate_initial_poses
from lightdock.constants import DEFAULT_POSITIONS_FOLDER, DEFAULT_SWARM_FOLDER, DEFAULT_LIST_EXTENSION, \
    DEFAULT_LIGHTDOCK_PREFIX, DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG, DEFAULT_REC_NM_FILE, DEFAULT_LIG_NM_FILE, MIN_EXTENT, MAX_EXTENT
from lightdock.util.logger import LoggingManager
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
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


def check_starting_file(file_name, glowworms, use_anm):
    """Check if the file_name file contains the required starting coordinates"""
    with open(file_name) as input_file:
        lines = [line.rstrip(os.linesep) for line in input_file.readlines()]
        num_glowworms = 0
        expected_size = 7 + DEFAULT_NMODES_REC + DEFAULT_NMODES_LIG if use_anm else 7
        for line in lines:
            fields = line.split()
            if len(fields) != expected_size:
                return False
            num_glowworms += 1
        return num_glowworms == glowworms


def calculate_starting_positions(receptor, ligand, swarms, glowworms, starting_points_seed, 
    ftdock_file=None, use_anm=False, anm_seed=0):
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
                                                        anm_seed, DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG)
        log.info("Generated %d positions files" % len(starting_points_files))
    else:
        log.warning("Folder %s already exists, skipping calculation" % init_folder)
        starting_points_files = glob.glob('init/initial_positions*.dat')
        if len(starting_points_files) != swarms:
            raise LightDockError("The number of initial positions files does not correspond with the number of swarms")
        for starting_point_file in starting_points_files:
            if not check_starting_file(starting_point_file, glowworms, use_anm):
                raise LightDockError("Error reading starting coordinates from file %s" % starting_point_file)
    log.info("Done.")
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