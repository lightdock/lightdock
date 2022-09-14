"""Simulation common functions"""

import os
import glob
import json
import time
import numpy as np
from pathlib import Path
from lightdock.prep.poses import calculate_initial_poses
from lightdock.constants import (
    DEFAULT_POSITIONS_FOLDER,
    DEFAULT_SWARM_FOLDER,
    DEFAULT_LIST_EXTENSION,
    DEFAULT_LIGHTDOCK_PREFIX,
    DEFAULT_NMODES_REC,
    DEFAULT_NMODES_LIG,
    MIN_EXTENT,
    MAX_EXTENT,
    DEFAULT_SETUP_FILE,
    DEFAULT_LIGHTDOCK_INFO,
    DEFAULT_STARTING_PREFIX,
    MAX_TRANSLATION,
    MAX_ROTATION,
    DEFAULT_SWARM_RADIUS,
    DEFAULT_MASK_FILE,
    DEFAULT_SWARM_DISTANCE,
    DEFAULT_SWARMS_PER_RESTRAINT,
)
from lightdock.util.logger import LoggingManager
from lightdock.pdbutil.PDBIO import parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.structure.nm import calculate_nmodes, write_nmodes
from lightdock.gso.boundaries import Boundary, BoundingBox
from lightdock.error.lightdock_errors import LightDockError
from lightdock.version import CURRENT_VERSION


log = LoggingManager.get_logger("lightdock3_setup")


def get_pdb_files(input_file):
    """Get a list of the PDB files in the input_file"""
    file_names = []
    with open(input_file) as handle:
        for line in handle:
            file_name = line.rstrip(os.linesep)
            if Path(file_name).exists():
                file_names.append(file_name)
    return file_names


def read_input_structure(
    pdb_file_name,
    ignore_oxt=True,
    ignore_hydrogens=False,
    ignore_water=False,
    verbose_parser=False,
):
    """Reads the input structure.

    The arguments pdb_file_name can be a PDB file or a file
    containing a list of PDB files.

    ignore_oxt flag avoids saving OXT atoms.
    """
    atoms_to_ignore = []
    residues_to_ignore = []
    if ignore_oxt:
        atoms_to_ignore.append("OXT")
        log.info("Ignoring OXT atoms")
    if ignore_hydrogens:
        atoms_to_ignore.append("H")
        log.info("Ignoring Hydrogen atoms")
    if ignore_water:
        residues_to_ignore.append("HOH")
        log.info("Ignoring water")

    structures = []
    file_names = []
    file_name, file_extension = os.path.splitext(pdb_file_name)
    if file_extension == DEFAULT_LIST_EXTENSION:
        file_names.extend(get_pdb_files(pdb_file_name))
    else:
        file_names.append(pdb_file_name)
    for file_name in file_names:
        log.info(f"Reading structure from {file_name} PDB file...")
        atoms, residues, chains = parse_complex_from_file(
            file_name, atoms_to_ignore, residues_to_ignore, verbose_parser
        )
        structures.append(
            {
                "atoms": atoms,
                "residues": residues,
                "chains": chains,
                "file_name": file_name,
            }
        )
        log.info(f"{len(atoms)} atoms, {len(residues)} residues read.")

    # Representatives are now the first structure, but this could change in the future
    structure = Complex.from_structures(structures)
    return structure


def write_mask_to_file(nm_mask, mask_file_name):
    """Saves the indexes of atoms involved in ANM"""
    np.save(mask_file_name, nm_mask)


def save_lightdock_structure(structure):
    """Saves the structure parsed by LightDock"""
    log.info("Saving processed structure to PDB file...")
    for structure_index, file_name in enumerate(structure.structure_file_names):
        moved_file_name = Path(file_name).parent / Path(
            DEFAULT_LIGHTDOCK_PREFIX % Path(file_name).name
        )
        if moved_file_name.exists():
            raise LightDockError(
                f"{moved_file_name} already exists, please delete previous setup generated files"
            )
        write_pdb_to_file(structure, moved_file_name, structure[structure_index])
        mask_file_name = Path(file_name).parent / Path(
            DEFAULT_MASK_FILE % Path(file_name).stem
        )
        write_mask_to_file(structure.nm_mask, mask_file_name)
    log.info("Done.")


def calculate_anm(structure, num_nmodes, rmsd, seed, file_name):
    """Calculates ANM for representative structure"""
    original_file_name = structure.structure_file_names[structure.representative_id]
    # We have to use the parsed structure by LightDock
    parsed_lightdock_structure = Path(original_file_name).parent / Path(
        DEFAULT_LIGHTDOCK_PREFIX % Path(original_file_name).name
    )
    modes = calculate_nmodes(parsed_lightdock_structure, num_nmodes, rmsd, seed, structure)
    structure.n_modes = modes
    write_nmodes(modes, file_name)
    log.info(f"{num_nmodes} normal modes calculated")


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


def calculate_starting_positions(
    receptor,
    ligand,
    swarms,
    glowworms,
    starting_points_seed,
    receptor_restraints,
    ligand_restraints,
    rec_translation,
    lig_translation,
    surface_density,
    use_anm=False,
    anm_seed=0,
    anm_rec=DEFAULT_NMODES_REC,
    anm_lig=DEFAULT_NMODES_LIG,
    is_membrane=False,
    is_transmembrane=False,
    write_starting_positions=False,
    swarm_radius=DEFAULT_SWARM_RADIUS,
    flip=False,
    swarms_at_fixed_distance=DEFAULT_SWARM_DISTANCE,
    swarms_per_restraint=DEFAULT_SWARMS_PER_RESTRAINT,
    dense_sampling=False
):
    """Defines the starting positions of each glowworm in the simulation.

    If the init folder already exists, uses the starting positions from this folder.
    """
    log.info("Calculating starting positions...")
    log.info(f"  * Surface density: TotalSASA/{surface_density:.2f}")
    log.info(f"  * Swarm radius: {swarm_radius:.2f} Å")
    log.info(f"  * 180° flip of 50% of starting poses: {flip}")
    init_folder = DEFAULT_POSITIONS_FOLDER
    if not os.path.isdir(init_folder):
        os.mkdir(init_folder)
        starting_points_files = calculate_initial_poses(
            receptor,
            ligand,
            swarms,
            glowworms,
            starting_points_seed,
            receptor_restraints,
            ligand_restraints,
            rec_translation,
            lig_translation,
            surface_density,
            init_folder,
            use_anm,
            anm_seed,
            anm_rec,
            anm_lig,
            is_membrane,
            is_transmembrane,
            write_starting_positions,
            swarm_radius,
            flip,
            swarms_at_fixed_distance,
            swarms_per_restraint,
            dense_sampling,
        )
        log.info(f"Generated {len(starting_points_files)} positions files")
    else:
        if receptor_restraints:
            log.warning(
                f"Folder {init_folder} already exists and restraints apply. Please check for consistency."
            )
        else:
            log.warning(f"Folder {init_folder} already exists, skipping calculation")

        pattern = str(
            Path(DEFAULT_POSITIONS_FOLDER) / Path(f"{DEFAULT_STARTING_PREFIX}*.dat")
        )
        starting_points_files = glob.glob(pattern)
        for starting_point_file in starting_points_files:
            if not check_starting_file(
                starting_point_file, glowworms, use_anm, anm_rec, anm_lig
            ):
                raise LightDockError(
                    f"Error reading starting coordinates from file {starting_point_file}"
                )
    log.info("Done.")
    return starting_points_files


def load_starting_positions(
    swarms, glowworms, use_anm, anm_rec=DEFAULT_NMODES_REC, anm_lig=DEFAULT_NMODES_LIG
):
    """Gets the list of starting positions of this simulation"""
    pattern = str(
        Path(DEFAULT_POSITIONS_FOLDER) / Path(f"{DEFAULT_STARTING_PREFIX}*.dat")
    )
    starting_points_files = sorted(glob.glob(pattern))
    if len(starting_points_files) != swarms:
        raise LightDockError(
            "The number of initial positions files does not correspond with the number of swarms"
        )
    for swarm_id in range(len(starting_points_files)):
        starting_point_file = Path(DEFAULT_POSITIONS_FOLDER) / Path(
            f"{DEFAULT_STARTING_PREFIX}_{swarm_id}.dat"
        )
        if not check_starting_file(
            starting_point_file, glowworms, use_anm, anm_rec, anm_lig
        ):
            raise LightDockError(
                f"Error reading starting coordinates from file {starting_point_file}"
            )
    return starting_points_files


def prepare_results_environment(swarms=10):
    """Prepares the folder structure required by the simulation"""
    log.info("Preparing environment")
    for id_swarm in range(swarms):
        saving_path = "%s%d" % (DEFAULT_SWARM_FOLDER, id_swarm)
        if Path(saving_path).is_dir():
            raise LightDockError(f"Simulation folder {saving_path} already exists")
        os.mkdir(saving_path)
    log.info("Done.")


def create_setup_file(args):
    """Dumps the args object into a setup file in JSON format"""
    with open(DEFAULT_SETUP_FILE, "w") as fp:
        json.dump(vars(args), fp, sort_keys=True, indent=4)


def get_setup_from_file(file_name):
    """Reads the simulation setup from the file_name"""
    with open(file_name) as input_file:
        return json.load(input_file)


def create_simulation_info_file(args, path=".", file_name=DEFAULT_LIGHTDOCK_INFO):
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
    data = {"start_time": now, "simulation_version": CURRENT_VERSION}
    data.update(vars(args))

    # Store the data in the file sorted alphabetically
    with open(output_file_name, "w") as fp:
        json.dump(data, fp, indent=4, sort_keys=True)

    return output_file_name


def get_default_box(use_anm, anm_rec, anm_lig):
    """Get the default bounding box"""
    boundaries = [
        Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
        Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
        Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
        Boundary(-MAX_ROTATION, MAX_ROTATION),
        Boundary(-MAX_ROTATION, MAX_ROTATION),
        Boundary(-MAX_ROTATION, MAX_ROTATION),
        Boundary(-MAX_ROTATION, MAX_ROTATION),
    ]
    if use_anm:
        boundaries.extend([Boundary(MIN_EXTENT, MAX_EXTENT) for _ in range(anm_rec)])
        boundaries.extend([Boundary(MIN_EXTENT, MAX_EXTENT) for _ in range(anm_lig)])

    return BoundingBox(boundaries)


def parse_restraints_file(restraints_file_name):
    """Parse a restraints file, returns a dictionary for receptor and ligand"""
    with open(restraints_file_name) as input_restraints:
        raw_restraints = [
            line.rstrip(os.linesep) for line in input_restraints.readlines()
        ]
        restraints = {
            "receptor": {"active": [], "passive": [], "blocked": []},
            "ligand": {"active": [], "passive": [], "blocked": []},
        }
        for restraint in raw_restraints:
            if restraint and restraint[0] in ["R", "L"]:
                try:
                    fields = restraint.split()
                    residue_identifier = fields[1].split(".")
                    # Only consider first character if many
                    chain_id = residue_identifier[0][0].upper()
                    # Only considering 3 chars if many
                    residue = residue_identifier[1][0:3].upper()
                    # Check for integer
                    try:
                        residue_number = int(residue_identifier[2])
                        residue_insertion = ""
                    except ValueError:
                        # Possible residue insertion
                        residue_number = int(residue_identifier[2][:-1])
                        residue_insertion = residue_identifier[2][-1].upper()
                    parsed_restraint = (
                        f"{chain_id}.{residue}.{residue_number}{residue_insertion}"
                    )
                    # Check type of restraint:
                    active = passive = blocked = False
                    try:
                        restraint_type = fields[2][0].upper()
                        passive = restraint_type == "P"
                        blocked = restraint_type == "B"
                        active = restraint_type == "A"
                    except (IndexError, AttributeError):
                        active = True

                    if fields[0] == "R":
                        if (
                            parsed_restraint not in restraints["receptor"]["active"]
                            and parsed_restraint
                            not in restraints["receptor"]["passive"]
                        ):
                            if active:
                                restraints["receptor"]["active"].append(
                                    parsed_restraint
                                )
                            elif passive:
                                restraints["receptor"]["passive"].append(
                                    parsed_restraint
                                )
                            elif blocked:
                                restraints["receptor"]["blocked"].append(
                                    parsed_restraint
                                )
                            else:
                                pass
                    else:
                        if (
                            parsed_restraint not in restraints["ligand"]["active"]
                            and parsed_restraint not in restraints["ligand"]["passive"]
                        ):
                            if active:
                                restraints["ligand"]["active"].append(parsed_restraint)
                            elif passive:
                                restraints["ligand"]["passive"].append(parsed_restraint)
                            elif blocked:
                                restraints["ligand"]["blocked"].append(parsed_restraint)
                            else:
                                pass
                except (AttributeError, IndexError):
                    log.warning(f"Ignoring malformed restraint {restraint}")

        return restraints


def get_restraints(structure, restraints):
    """Check for each restraint in the format Chain.ResidueName.ResidueNumber in
    restraints if they exist in the given structure.
    """
    residues = {"active": [], "passive": [], "blocked": []}
    for restraint_type in ["active", "passive", "blocked"]:
        for restraint in restraints[restraint_type]:
            chain_id, residue_name, residue_number = restraint.split(".")
            if residue_number[-1].isalpha():
                residue_insertion = residue_number[-1]
                residue_number = residue_number[:-1]
            else:
                residue_insertion = ""
            residue = structure.get_residue(
                chain_id, residue_name, residue_number, residue_insertion
            )
            if not residue:
                raise LightDockError(f"Restraint {restraint} not found in structure")
            residues[restraint_type].append(residue)
    return residues
