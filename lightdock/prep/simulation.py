import os
from lightdock.prep.poses import calculate_initial_poses
from lightdock.constants import DEFAULT_POSITIONS_FOLDER, DEFAULT_CLUSTER_FOLDER, DEFAULT_LIST_EXTENSION, \
    DEFAULT_LIGHTDOCK_PREFIX, DEFAULT_NMODES_REC, DEFAULT_NMODES_LIG, DEFAULT_REC_NM_FILE, DEFAULT_LIG_NM_FILE, MIN_EXTENT, MAX_EXTENT
from lightdock.util.logger import LoggingManager
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.error.lightdock_errors import LightDockError


log = LoggingManager.get_logger('lightdock')


def get_pdb_files(input_file):
    """Get a list of the PDB files in the input_file"""
    file_names = []
    with open(input_file) as input:
        for line in input:
            file_name = line.rstrip(os.linesep)
            if os.path.exists(file_name):
                file_names.append(file_name)
    return file_names


def read_input_structures(parser):
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
        log.info("Reading %s receptor PDB file..." % file_name)
        atoms, residues, chains = parse_complex_from_file(file_name, atoms_to_ignore)
        receptor_structures.append({'atoms': atoms, 'residues': residues, 'chains': chains, 'file_name': file_name})
        log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    ligand_structures = []
    file_names = []
    file_name, file_extension = os.path.splitext(parser.args.info_ligand_pdb)
    if file_extension == DEFAULT_LIST_EXTENSION:
        file_names.extend(get_pdb_files(parser.args.info_ligand_pdb))
    else:
        file_names.append(parser.args.info_ligand_pdb)
    for file_name in file_names:
        log.info("Reading %s ligand PDB file..." % file_name)
        atoms, residues, chains = parse_complex_from_file(file_name, atoms_to_ignore)
        ligand_structures.append({'atoms': atoms, 'residues': residues, 'chains': chains, 'file_name': file_name})
        log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

    # Representatives are now the first structure, but this could change in the future
    receptor = Complex.from_structures(receptor_structures)
    ligand = Complex.from_structures(ligand_structures)
    return receptor, ligand


def save_lightdock_structures(receptor, ligand):
    """Saves the structures parsed by LightDock"""
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
    """Defines the starting positions of each glowworm in the simulation.

    If the init_folder already exists, uses the starting positions from this folder.
    """
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