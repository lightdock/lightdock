"""Simulation common functions"""

import os
import glob
import json
import time
from pathlib import Path
from lightdock.util.logger import LoggingManager
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file, create_pdb_from_points
from lightdock.prep.geometry import create_bild_file
from lightdock.structure.complex import Complex
from lightdock.structure.nm import calculate_nmodes, write_nmodes
from lightdock.gso.boundaries import Boundary, BoundingBox
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.mathutil.lrandom import MTGenerator, NormalGenerator
from lightdock.error.lightdock_errors import LightDockError
from lightdock.prep.poses import get_random_point_within_sphere, create_file_from_poses
from symlightdock.constants import DEFAULT_POSITIONS_FOLDER, DEFAULT_TRIAL_FOLDER, DEFAULT_LIST_EXTENSION, \
                                   DEFAULT_SYMLIGHTDOCK_PREFIX, DEFAULT_NMODES, \
                                   MIN_EXTENT, MAX_EXTENT, DEFAULT_SETUP_FILE, DEFAULT_SYMLIGHTDOCK_INFO, \
                                   DEFAULT_STARTING_PREFIX, MAX_TRANSLATION, MAX_ROTATION, \
                                   DEFAULT_EXTENT_MU, DEFAULT_EXTENT_SIGMA, DEFAULT_PDB_STARTING_PREFIX, \
                                   DEFAULT_BILD_STARTING_PREFIX

log = LoggingManager.get_logger('symlightdock_setup')


def get_pdb_files(input_file):
    """Get a list of the PDB files in the input_file"""
    file_names = []
    with open(input_file) as handle:
        for line in handle:
            file_name = line.rstrip(os.linesep)
            if Path(file_name).exists():
                file_names.append(file_name)
    return file_names


def read_input_structure(pdb_file_name, ignore_oxt=True, ignore_hydrogens=False, verbose_parser=False):
    """Reads the input structure.

    The arguments pdb_file_name can be a PDB file or a file
    containing a list of PDB files.

    ignore_oxt flag avoids saving OXT atoms.
    """
    atoms_to_ignore = []
    if ignore_oxt:
        atoms_to_ignore.append("OXT")
        log.info("Ignoring OXT atoms")
    if ignore_hydrogens:
        atoms_to_ignore.append("H")
        log.info("Ignoring Hydrogen atoms")

    structures = []
    file_names = []
    file_name, file_extension = os.path.splitext(pdb_file_name)
    if file_extension == DEFAULT_LIST_EXTENSION:
        file_names.extend(get_pdb_files(pdb_file_name))
    else:
        file_names.append(pdb_file_name)
    for file_name in file_names:
        log.info(f"Reading structure from {file_name} PDB file...")
        atoms, residues, chains = parse_complex_from_file(file_name, atoms_to_ignore, verbose_parser)
        structures.append({'atoms': atoms, 'residues': residues, 'chains': chains, 'file_name': file_name})
        log.info(f"{len(atoms)} atoms, {len(residues)} residues read.")

    # Representatives are now the first structure, but this could change in the future
    structure = Complex.from_structures(structures)
    return structure


def save_lightdock_structure(structure):
    """Saves the structure parsed by LightDock"""
    log.info("Saving processed structure to PDB file...")
    for structure_index, file_name in enumerate(structure.structure_file_names):
        moved_file_name = Path(file_name).parent / Path(DEFAULT_SYMLIGHTDOCK_PREFIX % Path(file_name).name)
        if moved_file_name.exists():
            raise LightDockError(f"{moved_file_name} already exists, please delete previous setup generated files")
        write_pdb_to_file(structure, moved_file_name, structure[structure_index])
    log.info("Done.")


def calculate_anm(structure, num_nmodes, file_name):
    """Calculates ANM for representative structure"""
    original_file_name = structure.structure_file_names[structure.representative_id]
    # We have to use the parsed structure by LightDock
    parsed_lightdock_structure = Path(original_file_name).parent / \
                                 Path(DEFAULT_SYMLIGHTDOCK_PREFIX % Path(original_file_name).name)
    modes = calculate_nmodes(parsed_lightdock_structure, num_nmodes, structure)
    structure.n_modes = modes
    write_nmodes(modes, file_name)
    log.info(f"{num_nmodes} normal modes calculated")


def check_starting_file(file_name, glowworms, use_anm, num_anm):
    """Check if the file_name file contains the required starting coordinates"""
    with open(file_name) as input_file:
        lines = [line.rstrip(os.linesep) for line in input_file.readlines()]
        num_glowworms = 0
        expected_size = 5 + num_anm if use_anm else 5
        for line in lines:
            fields = line.split()
            if len(fields) != expected_size:
                return False
            num_glowworms += 1
        return num_glowworms == glowworms


def calculate_initial_poses(structure, num_trials, num_glowworms,
                            seed, restraints, translation,
                            dest_folder, nm_mode=False, nm_seed=0, num_anm=0,
                            writing_starting_positions=False):
    """Calculates the starting points for each of the glowworms"""

    # Random number generator for poses
    rng = MTGenerator(seed)

    # Random number generator for NM
    if nm_mode:
        rng_nm = NormalGenerator(nm_seed, mu=DEFAULT_EXTENT_MU, sigma=DEFAULT_EXTENT_SIGMA)
    else:
        rng_nm = None

    positions_files = []
    for trial_id in range(num_trials):
        poses = []
        for _ in range(num_glowworms):
            x = rng(0., 10.)

            q = Quaternion.random(rng)

            # Glowworm's optimization vector
            op_vector = [x, q.w, q.x, q.y, q.z]

            # If ANM is enabled, we need to create random components for the extents
            if nm_mode:
                if num_anm > 0:
                    op_vector.extend([rng_nm() for _ in range(num_anm)])

            poses.append(op_vector)

        if writing_starting_positions:
            # Save poses as pdb file
            pdb_file_name = os.path.join(dest_folder, '%s_%s.pdb' % (DEFAULT_PDB_STARTING_PREFIX, trial_id))
            create_pdb_from_points(pdb_file_name, [[pose[0], pose[1], pose[2]] for pose in poses[:num_glowworms]])
            # Generate bild files for glowworm orientations
            bild_file_name = os.path.join(dest_folder, '%s_%s.bild' % (DEFAULT_BILD_STARTING_PREFIX, trial_id))
            create_bild_file(bild_file_name, poses)

        # Save poses as initial_positions file
        pos_file_name = os.path.join(dest_folder, '%s_%s.dat' % (DEFAULT_STARTING_PREFIX, trial_id))
        create_file_from_poses(pos_file_name, poses[:num_glowworms])
        positions_files.append(pos_file_name)

    return positions_files



def calculate_starting_positions(structure, trials, glowworms, starting_points_seed,
    restraints, translation, use_anm=False, anm_seed=0, num_anm=DEFAULT_NMODES,
    write_starting_positions=False):
    """Defines the starting positions of each glowworm in the simulation.

    If the init folder already exists, uses the starting positions from this folder.
    """
    log.info("Calculating starting positions...")
    init_folder = DEFAULT_POSITIONS_FOLDER
    if not os.path.isdir(init_folder):
        os.mkdir(init_folder)
        starting_points_files = calculate_initial_poses(structure, trials, glowworms,
                                                        starting_points_seed,
                                                        restraints, translation,
                                                        init_folder, use_anm, anm_seed, num_anm,
                                                        write_starting_positions)
        log.info(f"Generated {len(starting_points_files)} positions files")
    else:
        if restraints:
            log.warning(f"Folder {init_folder} already exists and restraints apply. Please check for consistency.")
        else:
            log.warning(f"Folder {init_folder} already exists, skipping calculation")

        pattern = str(Path(DEFAULT_POSITIONS_FOLDER) / Path(f"{DEFAULT_STARTING_PREFIX}*.dat"))
        starting_points_files = glob.glob(pattern)
        for starting_point_file in starting_points_files:
            if not check_starting_file(starting_point_file, glowworms, use_anm, num_anm):
                raise LightDockError(f"Error reading starting coordinates from file {starting_point_file}")
    log.info("Done.")
    return starting_points_files


def load_starting_positions(trials, glowworms, use_anm, num_anm=DEFAULT_NMODES):
    """Gets the list of starting positions of this simulation"""
    pattern = str(Path(DEFAULT_POSITIONS_FOLDER) / Path(f"{DEFAULT_STARTING_PREFIX}*.dat"))
    starting_points_files = sorted(glob.glob(pattern))
    if len(starting_points_files) != trials:
        raise LightDockError("The number of initial positions files does not correspond with the number of trials")
    for trial_id in range(len(starting_points_files)):
        starting_point_file = Path(DEFAULT_POSITIONS_FOLDER) / Path(f"{DEFAULT_STARTING_PREFIX}_{trial_id}.dat")
        if not check_starting_file(starting_point_file, glowworms, use_anm, num_anm):
            raise LightDockError(f"Error reading starting coordinates from file {starting_point_file}")
    return starting_points_files


def prepare_results_environment(trials=10):
    """Prepares the folder structure required by the simulation"""
    log.info("Preparing environment")
    for trial_id in range(trials):
        saving_path = f"{DEFAULT_TRIAL_FOLDER}{trial_id}"
        if Path(saving_path).is_dir():
            raise LightDockError(f"Simulation folder {saving_path} already exists")
        os.mkdir(saving_path)
    log.info("Done.")


def create_setup_file(args):
    """Dumps the args object into a setup file in JSON format"""
    with open(DEFAULT_SETUP_FILE, 'w') as fp:
        json.dump(vars(args), fp, sort_keys=True, indent=4)


def get_setup_from_file(file_name):
    """Reads the simulation setup from the file_name"""
    with open(file_name) as input_file:
        return json.load(input_file)


def create_simulation_info_file(args, path='.', file_name=DEFAULT_SYMLIGHTDOCK_INFO):
    """Creates a simulation file from which recover from in a new simulation"""
    # Create the simulation info file. If it exists, includes a number
    # in the extension to avoid collision
    output_file_name = Path(path) / Path(file_name)
    if output_file_name.is_file():
        original_file_name = output_file_name
        i = 1
        while os.path.isfile(output_file_name):
            output_file_name = f"{original_file_name}.{i}"
            i += 1

    # Data to store
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    data = {'start_time': now}
    data.update(vars(args))

    # Store the data in the file sorted alphabetically
    with open(output_file_name, 'w') as fp:
        json.dump(vars(args), fp, indent=4, sort_keys=True)

    return output_file_name


def get_default_box(use_anm, anm_rec, anm_lig):
    """Get the default bounding box"""
    boundaries = [Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                  Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                  Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                  Boundary(-MAX_ROTATION, MAX_ROTATION),
                  Boundary(-MAX_ROTATION, MAX_ROTATION),
                  Boundary(-MAX_ROTATION, MAX_ROTATION),
                  Boundary(-MAX_ROTATION, MAX_ROTATION)]
    if use_anm:
        boundaries.extend([Boundary(MIN_EXTENT, MAX_EXTENT) for _ in range(anm_rec)])
        boundaries.extend([Boundary(MIN_EXTENT, MAX_EXTENT) for _ in range(anm_lig)])

    return BoundingBox(boundaries)


def parse_restraints_file(restraints_file_name):
    """Parse a restraints file, returns a dictionary for receptor and ligand"""
    with open(restraints_file_name) as input_restraints:
        raw_restraints = [line.rstrip(os.linesep) for line in input_restraints.readlines()]
        restraints = {'receptor': {'active':[], 'passive':[], 'blocked':[]},
                      'ligand': {'active':[], 'passive':[], 'blocked':[]}}
        for restraint in raw_restraints:
            if restraint and restraint[0] in ['R', 'L']:
                try:
                    fields = restraint.split()
                    residue_identifier = fields[1].split('.')
                    # Only consider first character if many
                    chain_id = residue_identifier[0][0].upper()
                    # Only considering 3 chars if many
                    residue = residue_identifier[1][0:3].upper()
                    # Check for integer
                    residue_number = int(residue_identifier[2])
                    parsed_restraint = f"{chain_id}.{residue}.{residue_number}"
                    # Check type of restraint:
                    active = passive = blocked = False
                    try:
                        restraint_type = fields[2][0].upper()
                        passive = (restraint_type == 'P')
                        blocked = (restraint_type == 'B')
                        active = (restraint_type == 'A')
                    except (IndexError, AttributeError):
                        active = True

                    if fields[0] == 'R':
                        if parsed_restraint not in restraints['receptor']['active'] and \
                            parsed_restraint not in restraints['receptor']['passive']:
                            if active:
                                restraints['receptor']['active'].append(parsed_restraint)
                            elif passive:
                                restraints['receptor']['passive'].append(parsed_restraint)
                            elif blocked:
                                restraints['receptor']['blocked'].append(parsed_restraint)
                            else:
                                pass
                    else:
                        if parsed_restraint not in restraints['ligand']['active'] and \
                            parsed_restraint not in restraints['ligand']['passive']:
                            if active:
                                restraints['ligand']['active'].append(parsed_restraint)
                            elif passive:
                                restraints['ligand']['passive'].append(parsed_restraint)
                            elif blocked:
                                restraints['ligand']['blocked'].append(parsed_restraint)
                            else:
                                pass
                except (AttributeError, IndexError):
                    log.warning(f"Ignoring malformed restraint {restraint}")

        return restraints


def get_restraints(structure, restraints):
    """Check for each restraint in the format Chain.ResidueName.ResidueNumber in
    restraints if they exist in the given structure.
    """
    residues = {'active':[], 'passive':[], 'blocked':[]}
    for restraint_type in ['active', 'passive', 'blocked']:
        for restraint in restraints[restraint_type]:
            chain_id, residue_name, residue_number = restraint.split('.')
            residue = structure.get_residue(chain_id, residue_name, residue_number)
            if not residue:
                raise LightDockError(f"Restraint {restraint} not found in structure")
            residues[restraint_type].append(residue)
    return residues
