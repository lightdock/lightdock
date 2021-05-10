"""Module to parse the SymLightDock program command line options"""

import argparse
from lightdock.util.parser import valid_file, valid_integer_number, valid_natural_number, valid_float_number
from symlightdock.constants import DEFAULT_NUM_TRIALS, DEFAULT_NUM_GLOWWORMS, STARTING_POINTS_SEED, \
                                   STARTING_NM_SEED, DEFAULT_NMODES, GSO_SEED, DEFAULT_TRANSLATION_STEP, \
                                   DEFAULT_ROTATION_STEP, DEFAULT_NMODES_STEP
from lightdock.version import CURRENT_VERSION


class SetupCommandLineParser(object):
    """Parses the command line of symlightdock_setup"""
    def __init__(self, input_args=None):
        parser = argparse.ArgumentParser(prog="symlightdock_setup")

        # Structure
        parser.add_argument("structure_pdb", help="Monomer structure PDB file",
                            type=valid_file, metavar="stucture_pdb")
        # Swarms
        parser.add_argument("-t", "--t", "--trials", help="Number of trials",
                            dest="trials", type=valid_natural_number, default=DEFAULT_NUM_TRIALS)
        # Glowworms
        parser.add_argument("-g", "--g", "--glowworms", help="Number of glowworms per swarm",
                            dest="glowworms", type=valid_integer_number, default=DEFAULT_NUM_GLOWWORMS)
        # Starting points seed
        parser.add_argument("--seed_points", help="Random seed used in starting positions calculation",
                            dest="starting_points_seed", type=int, default=STARTING_POINTS_SEED)
        # Dealing with OXT atoms
        parser.add_argument("--noxt", "-noxt", help="Remove OXT atoms",
                            dest="noxt", action='store_true', default=False)
        # Dealing with hydrogen atoms
        parser.add_argument("--noh", "-noh", help="Remove Hydrogen atoms",
                            dest="noh", action='store_true', default=False)
        # Verbose PDB parser
        parser.add_argument("--verbose_parser", "-verbose_parser", help="PDB parsing verbose mode",
                            dest="verbose_parser", action='store_true', default=False)
        # Normal modes
        parser.add_argument("-anm", "--anm", help="Activates the use of ANM backbone flexibility",
                            dest="use_anm", action='store_true', default=False)
        # Normal modes extent seed
        parser.add_argument("--seed_anm", "-seed_anm", help="Random seed used in ANM intial extent",
                            dest="anm_seed", type=int, default=STARTING_NM_SEED)
        parser.add_argument("-num_anm", "--num_anm", help="Number of ANM modes for the structure",
                            type=valid_natural_number, dest="num_anm", default=DEFAULT_NMODES)

        # Restraints file
        parser.add_argument("-r", "--r", "-rst", "--rst", help="Restraints file",
                            dest="restraints", type=valid_file,
                            metavar="restraints", default=None)
        # Create bild files
        parser.add_argument("-sp", "--sp", "-starting_positions", "--starting_positions",
                            help="Enables the generation of support files in init directory",
                            dest="write_starting_positions", action='store_true', default=False)
        if input_args:
            self.args = parser.parse_args(input_args)
        else:
            self.args = parser.parse_args()


class CommandLineParser(object):
    """Parses the command line"""
    def __init__(self):
        parser = argparse.ArgumentParser(prog="symlightdock")

        # Receptor
        parser.add_argument("setup_file", help="Setup file name",
                            type=valid_file, metavar="setup_file")
        # Steps
        parser.add_argument("steps", help="number of steps of the simulation",
                            type=valid_integer_number)
        # Configuration file
        parser.add_argument("-f", "--file", help="algorithm configuration file",
                            dest="configuration_file", type=valid_file,
                            metavar="configuration_file")
        # Scoring function
        parser.add_argument("-s", "--scoring_function", help="scoring function used",
                            dest="scoring_function")
        # GSO seed
        parser.add_argument("-sg", "--seed", help="seed used in the algorithm",
                            dest="gso_seed", type=int, default=GSO_SEED)
        # Translation step
        parser.add_argument("-t", "--translation_step", help="translation step",
                            type=valid_float_number, default=DEFAULT_TRANSLATION_STEP)
        # Rotation step
        parser.add_argument("-r", "--rotation_step", help="normalized rotation step",
                            type=valid_float_number, default=DEFAULT_ROTATION_STEP)
        # Version
        parser.add_argument("-V", "-v", "--version", help="show version",
                            action='version', version="%s %s" % (parser.prog, CURRENT_VERSION))
        # Number of cpu cores to use
        parser.add_argument("-c", "--cores", help="number of cpu cores to use",
                            dest="cores", type=valid_integer_number)
        # Profiling
        parser.add_argument("--profile", help="Output profiling data",
                            dest="profiling", action='store_true', default=False)
        # MPI parallel execution
        parser.add_argument("-mpi", "--mpi", help="activates the MPI parallel execution",
                            dest="mpi", action='store_true', default=False)
        # Normal modes step
        parser.add_argument("-ns", "--nmodes_step", help="normalized normal modes step",
                            type=valid_float_number, default=DEFAULT_NMODES_STEP)
        # List of trials to simulate
        parser.add_argument("-l", "--list", help="List of trials to simulate",
                            dest="trial_list", nargs="+", type=int, required=False)

        self.args = parser.parse_args()
