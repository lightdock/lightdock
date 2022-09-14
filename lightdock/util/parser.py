"""Module to parse the LightDock program command line options"""

import argparse
import os
from lightdock.constants import (
    DEFAULT_NUM_GLOWWORMS,
    GSO_SEED,
    STARTING_POINTS_SEED,
    DEFAULT_TRANSLATION_STEP,
    DEFAULT_ROTATION_STEP,
    STARTING_NM_SEED,
    DEFAULT_NMODES_STEP,
    DEFAULT_LIST_EXTENSION,
    DEFAULT_LIGHTDOCK_PREFIX,
    DEFAULT_NMODES_REC,
    DEFAULT_NMODES_LIG,
    DEFAULT_SURFACE_DENSITY,
    DEFAULT_SWARM_RADIUS,
    DEFAULT_ANM_RMSD,
    DEFAULT_SWARM_DISTANCE,
    DEFAULT_SWARMS_PER_RESTRAINT,
)
from lightdock.error.lightdock_errors import LightDockError
from lightdock.version import CURRENT_VERSION


def get_lightdock_structures(input_file):
    """Get a list of the PDB files in the input_file"""
    input_file_name, input_file_extension = os.path.splitext(input_file)
    file_names = []
    if input_file_extension == DEFAULT_LIST_EXTENSION:
        with open(input_file) as input_lines:
            for line in input_lines:
                file_name = line.rstrip(os.linesep)
                lightdock_structure = os.path.join(
                    os.path.dirname(file_name),
                    DEFAULT_LIGHTDOCK_PREFIX % os.path.basename(file_name),
                )
                if os.path.exists(lightdock_structure):
                    file_names.append(lightdock_structure)
    else:
        file_name = input_file
        lightdock_structure = os.path.join(
            os.path.dirname(file_name),
            DEFAULT_LIGHTDOCK_PREFIX % os.path.basename(file_name),
        )
        if os.path.exists(lightdock_structure):
            file_names.append(lightdock_structure)
        else:
            raise LightDockError(f"Structure file {lightdock_structure} not found")
    return file_names


# Validators
def valid_file(file_name):
    if not os.path.exists(file_name):
        raise argparse.ArgumentTypeError("The file does not exist")
    return file_name


def valid_integer_number(int_value):
    try:
        int_value = int(int_value)
    except:
        raise argparse.ArgumentTypeError(f"{int_value} is an invalid value")
    if int_value <= 0:
        raise argparse.ArgumentTypeError(f"{int_value} is an invalid value")
    return int_value


def valid_natural_number(int_value):
    try:
        int_value = int(int_value)
    except:
        raise argparse.ArgumentTypeError(f"{int_value} is an invalid value")
    if int_value < 0:
        raise argparse.ArgumentTypeError(f"{int_value} is an invalid value")
    return int_value


def valid_float_number(float_value):
    try:
        float_value = float(float_value)
    except:
        raise argparse.ArgumentTypeError(f"{float_value} is an invalid value")
    if float_value <= 0.0:
        raise argparse.ArgumentTypeError(f"{float_value} is an invalid value")
    return float_value


class SetupCommandLineParser(object):
    """Parses the command line of lightdock_setup"""

    def __init__(self, input_args=None):
        parser = argparse.ArgumentParser(prog="lightdock3_setup")

        # Receptor
        parser.add_argument(
            "receptor_pdb",
            help="Receptor structure PDB file",
            type=valid_file,
            metavar="receptor_pdb_file",
        )
        # Ligand
        parser.add_argument(
            "ligand_pdb",
            help="Ligand structure PDB file",
            type=valid_file,
            metavar="ligand_pdb_file",
        )
        # Clusters
        parser.add_argument(
            "-s",
            "--s",
            "--swarms",
            help="Fixed number of swarms of the simulation",
            dest="swarms",
            type=valid_integer_number,
            default=0,
        )
        # Glowworms
        parser.add_argument(
            "-g",
            "--g",
            "--glowworms",
            help="Number of glowworms per swarm",
            dest="glowworms",
            type=valid_integer_number,
            default=DEFAULT_NUM_GLOWWORMS,
        )
        # Starting points seed
        parser.add_argument(
            "--seed_points",
            help="Random seed used in starting positions calculation",
            dest="starting_points_seed",
            type=int,
            default=STARTING_POINTS_SEED,
        )
        # Dealing with OXT atoms
        parser.add_argument(
            "--noxt",
            "-noxt",
            help="Remove OXT atoms",
            dest="noxt",
            action="store_true",
            default=False,
        )
        # Dealing with hydrogen atoms
        parser.add_argument(
            "--noh",
            "-noh",
            help="Remove Hydrogen atoms",
            dest="noh",
            action="store_true",
            default=False,
        )
        # Dealing with hydrogen atoms
        parser.add_argument(
            "--now",
            "-now",
            help="Remove H2O",
            dest="now",
            action="store_true",
            default=False,
        )
        # Verbose PDB parser
        parser.add_argument(
            "--verbose_parser",
            "-verbose_parser",
            help="PDB parsing verbose mode",
            dest="verbose_parser",
            action="store_true",
            default=False,
        )
        # Normal modes
        parser.add_argument(
            "-anm",
            "--anm",
            help="Activates the use of ANM backbone flexibility",
            dest="use_anm",
            action="store_true",
            default=False,
        )
        # Normal modes extent seed
        parser.add_argument(
            "--seed_anm",
            "-seed_anm",
            help="Random seed used in ANM intial extent",
            dest="anm_seed",
            type=int,
            default=STARTING_NM_SEED,
        )
        # ANM RMSD interval for receptor
        parser.add_argument(
            "--anm_rec_rmsd",
            "-anm_rec_rmsd",
            help="RMSD interval to generate random ANM conformations",
            dest="anm_rec_rmsd",
            type=valid_float_number,
            default=DEFAULT_ANM_RMSD,
        )
        # ANM RMSD interval for ligand
        parser.add_argument(
            "--anm_lig_rmsd",
            "-anm_lig_rmsd",
            help="RMSD interval to generate random ANM conformations",
            dest="anm_lig_rmsd",
            type=valid_float_number,
            default=DEFAULT_ANM_RMSD,
        )
        # Default number of non-trivial modes for receptor
        parser.add_argument(
            "-ar",
            "--ar",
            "-anm_rec",
            "--anm_rec",
            help="Number of ANM modes for receptor",
            type=valid_natural_number,
            dest="anm_rec",
            default=DEFAULT_NMODES_REC,
        )
        # Default number of non-trivial modes for ligand
        parser.add_argument(
            "-al",
            "--al",
            "-anm_lig",
            "--anm_lig",
            help="Number of ANM modes for ligand",
            type=valid_natural_number,
            dest="anm_lig",
            default=DEFAULT_NMODES_LIG,
        )
        # Restraints file
        parser.add_argument(
            "-r",
            "--r",
            "-rst",
            "--rst",
            help="Restraints file",
            dest="restraints",
            type=valid_file,
            metavar="restraints",
            default=None,
        )
        # Membrane setup
        parser.add_argument(
            "-membrane",
            "--membrane",
            help="Enables the extra filter for membrane restraints",
            dest="membrane",
            action="store_true",
            default=False,
        )
        # Transmembrane setup
        parser.add_argument(
            "-transmembrane",
            "--transmembrane",
            help="Enables the extra filter for transmembrane restraints",
            dest="transmembrane",
            action="store_true",
            default=False,
        )
        # Create bild files
        parser.add_argument(
            "-sp",
            "--sp",
            "-starting_positions",
            "--starting_positions",
            help="Enables the generation of support files in init directory",
            dest="write_starting_positions",
            action="store_true",
            default=False,
        )
        # Surface density
        parser.add_argument(
            "-sd",
            "--sd",
            "-surface_density",
            "--surface_density",
            help="Surface density used for calculating the number of swarms",
            type=valid_float_number,
            dest="surface_density",
            default=DEFAULT_SURFACE_DENSITY,
        )
        # Swarm radius
        parser.add_argument(
            "-sr",
            "--sr",
            "-swarm_radius",
            "--swarm_radius",
            help="Swarm radius at generating glowworm poses (translation)",
            type=valid_float_number,
            dest="swarm_radius",
            default=DEFAULT_SWARM_RADIUS,
        )
        # Autoflip mode
        parser.add_argument(
            "-flip",
            "--flip",
            help="Activates the 180° flip of 50%% of starting poses when using restraints",
            dest="flip",
            action="store_true",
            default=False,
        )
        # Swarm fixed distance
        parser.add_argument(
            "-fd",
            "--fd",
            "-fixed_distance",
            "--fixed_distance",
            help="Use a fixed distance (Å) instead of using ligand to calculate distance of swarms to receptor surface",
            type=valid_float_number,
            dest="fixed_distance",
            default=DEFAULT_SWARM_DISTANCE,
        )
        # Number of swarms to keep per restraint
        parser.add_argument(
            "-spr",
            "--spr",
            "-swarms_per_restraint",
            "--swarms_per_restraint",
            help="Number of swarms to keep per restraint",
            type=valid_natural_number,
            dest="swarms_per_restraint",
            default=DEFAULT_SWARMS_PER_RESTRAINT,
        )
        # Enable dense sampling
        parser.add_argument(
            "--ds",
            "-ds",
            help="Enable dense sampling for restraints-swarm filtering",
            dest="dense_sampling",
            action="store_true",
            default=False,
        )
        # Version
        parser.add_argument(
            "-V",
            "-v",
            "--version",
            help="show version",
            action="version",
            version="%s %s" % (parser.prog, CURRENT_VERSION),
        )
        if input_args:
            self.args = parser.parse_args(input_args)
        else:
            self.args = parser.parse_args()


class ListScoringAction(argparse.Action):
    """Support class for listing the available scoring functions"""

    def __call__(self, parser, namespace, values, option_string=None):
        from lightdock import scoring

        scoring_path = os.path.dirname(scoring.__file__)
        scoring_functions = [
            s
            for s in os.listdir(scoring_path)
            if os.path.isdir(os.path.join(scoring_path, s)) and not s.startswith("__")
        ]
        print("Available scoring functions are: ", ", ".join(scoring_functions))
        raise SystemExit


class CommandLineParser(object):
    """Parses the command line"""

    def __init__(self, input_args=None):
        parser = argparse.ArgumentParser(prog="lightdock3")

        # Receptor
        parser.add_argument(
            "setup_file", help="Setup file name", type=valid_file, metavar="setup_file"
        )
        # Steps
        parser.add_argument(
            "steps", help="number of steps of the simulation", type=valid_integer_number
        )
        # Configuration file
        parser.add_argument(
            "-f",
            "--file",
            help="algorithm configuration file",
            dest="configuration_file",
            type=valid_file,
            metavar="configuration_file",
        )
        # Scoring function
        parser.add_argument(
            "-s",
            "--scoring_function",
            help="scoring function used",
            dest="scoring_function",
        )
        # GSO seed
        parser.add_argument(
            "-sg",
            "--seed",
            help="seed used in the algorithm",
            dest="gso_seed",
            type=int,
            default=GSO_SEED,
        )
        # Translation step
        parser.add_argument(
            "-t",
            "--translation_step",
            help="translation step",
            type=valid_float_number,
            default=DEFAULT_TRANSLATION_STEP,
        )
        # Rotation step
        parser.add_argument(
            "-r",
            "--rotation_step",
            help="normalized rotation step",
            type=valid_float_number,
            default=DEFAULT_ROTATION_STEP,
        )
        # Version
        parser.add_argument(
            "-V",
            "-v",
            "--version",
            help="show version",
            action="version",
            version="%s %s" % (parser.prog, CURRENT_VERSION),
        )
        # Number of cpu cores to use
        parser.add_argument(
            "-c",
            "--cores",
            help="number of cpu cores to use",
            dest="cores",
            type=valid_integer_number,
        )
        # Profiling
        parser.add_argument(
            "--profile",
            help="Output profiling data",
            dest="profiling",
            action="store_true",
            default=False,
        )
        # MPI parallel execution
        parser.add_argument(
            "-mpi",
            "--mpi",
            help="activates the MPI parallel execution",
            dest="mpi",
            action="store_true",
            default=False,
        )
        # Normal modes step
        parser.add_argument(
            "-ns",
            "--nmodes_step",
            help="normalized normal modes step",
            type=valid_float_number,
            default=DEFAULT_NMODES_STEP,
        )
        # Local minimization
        parser.add_argument(
            "-min",
            "--min",
            help="activates the local minimization",
            dest="local_minimization",
            action="store_true",
            default=False,
        )
        # List of available scoring functions
        parser.add_argument(
            "--listscoring",
            help="list all available scoring functions",
            action=ListScoringAction,
            nargs=0,
        )
        # List of swarms to simulate
        parser.add_argument(
            "-l",
            "--list",
            help="List of swarms to simulate",
            dest="swarm_list",
            nargs="+",
            type=int,
            required=False,
        )
        if input_args:
            self.args = parser.parse_args(input_args)
        else:
            self.args = parser.parse_args()
