"""Tests for the simulation module"""

import os
import shutil
from pathlib import Path
from glob import glob
from nose.tools import raises
from lightdock.error.lightdock_errors import LightDockError
from lightdock.prep.simulation import (
    parse_restraints_file,
    get_restraints,
    get_default_box,
    get_setup_from_file,
    create_setup_file,
    prepare_results_environment,
    get_pdb_files,
    read_input_structure,
    load_starting_positions,
    create_simulation_info_file,
    check_starting_file,
)
from lightdock.structure.complex import Complex
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.util.parser import SetupCommandLineParser, CommandLineParser
from lightdock.version import CURRENT_VERSION
from lightdock.test.support import compare_two_files


class TestParsingRestraintsFile:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_parsing_restraints"
        self.golden_data_path = self.path / "golden_data"

    def setup(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.mkdir(self.test_path)

    def teardown(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass

    def test_simple_restraints_file(self):
        input_file = self.golden_data_path / "rst_1.lst"

        restraints = parse_restraints_file(input_file)

        expected = {
            "receptor": {"active": ["A.ALA.1"], "passive": [], "blocked": []},
            "ligand": {"active": [], "passive": ["B.TYR.1"], "blocked": []},
        }

        assert restraints == expected

    def test_restraints_file_with_duplicate(self):
        input_file = self.golden_data_path / "rst_2.lst"

        restraints = parse_restraints_file(input_file)

        expected = {
            "receptor": {"active": ["A.ALA.1"], "passive": [], "blocked": []},
            "ligand": {"active": [], "passive": ["B.TYR.1"], "blocked": []},
        }

        assert restraints == expected

    def test_empty_restraints_file(self):
        input_file = self.golden_data_path / "rst_3.lst"

        restraints = parse_restraints_file(input_file)

        expected = {
            "receptor": {"active": [], "passive": [], "blocked": []},
            "ligand": {"active": [], "passive": [], "blocked": []},
        }

        assert restraints == expected

    def test_blocked_restraints_file(self):
        input_file = self.golden_data_path / "rst_4.lst"

        restraints = parse_restraints_file(input_file)

        expected = {
            "receptor": {"active": ["A.ALA.1"], "passive": [], "blocked": ["A.LYS.2"]},
            "ligand": {"active": [], "passive": ["B.TYR.1"], "blocked": ["B.TRP.2"]},
        }

        assert restraints == expected


class TestRestraints:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_restraints"
        self.golden_data_path = self.path / "golden_data"

    def setup(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.mkdir(self.test_path)

    def teardown(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass

    def test_get_restraints(self):
        input_file = self.golden_data_path / "2UUY_lig.pdb"
        _, _, chains = parse_complex_from_file(input_file)
        structure = Complex(chains)
        restraints = {"active": ["B.ALA.21"], "passive": ["B.GLY.75"], "blocked": []}

        residues = get_restraints(structure, restraints)

        assert len(residues["active"]) == 1 and len(residues["passive"]) == 1
        assert (
            residues["active"][0].name == "ALA" and residues["active"][0].number == 21
        )
        assert (
            residues["passive"][0].name == "GLY" and residues["passive"][0].number == 75
        )

    @raises(LightDockError)
    def test_get_restraints_with_error(self):
        input_file = self.golden_data_path / "2UUY_lig.pdb"
        _, _, chains = parse_complex_from_file(input_file)
        structure = Complex(chains)
        restraints = {"active": ["B.VAL.21"], "passive": ["B.GLY.75"], "blocked": []}

        residues = get_restraints(structure, restraints)

        assert residues is not None


class TestSimulation:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_simulation"
        self.golden_data_path = self.path / "golden_data"

    def setup(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.mkdir(self.test_path)

    def teardown(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass

    def test_bounding_box(self):
        box = get_default_box(use_anm=False, anm_rec=0, anm_lig=0)

        assert box.dimension == 7

    def test_bounding_box_anm(self):
        box = get_default_box(use_anm=True, anm_rec=5, anm_lig=3)

        assert box.dimension == 7 + 5 + 3

    def test_get_setup_from_file(self):
        read = get_setup_from_file(self.golden_data_path / "setup.json")

        expected = {
            "anm_lig": 10,
            "anm_lig_rmsd": 0.5,
            "anm_rec": 10,
            "anm_rec_rmsd": 0.5,
            "anm_seed": 324324,
            "dense_sampling": False,
            "fixed_distance": 0.0,
            "flip": False,
            "glowworms": 10,
            "ligand_pdb": "2UUY_lig.pdb",
            "membrane": False,
            "noh": False,
            "now": False,
            "noxt": True,
            "receptor_pdb": "2UUY_rec.pdb",
            "restraints": None,
            "setup_version": CURRENT_VERSION,
            "starting_points_seed": 324324,
            "surface_density": 50.0,
            "swarm_radius": 10.0,
            "swarms": 5,
            "swarms_per_restraint": 20,
            "transmembrane": False,
            "use_anm": True,
            "verbose_parser": False,
            "write_starting_positions": False,
        }

        assert read == expected

    def test_create_setup_file(self):
        shutil.copyfile(
            self.golden_data_path / "2UUY_rec.pdb", self.test_path / "2UUY_rec.pdb"
        )
        shutil.copyfile(
            self.golden_data_path / "2UUY_lig.pdb", self.test_path / "2UUY_lig.pdb"
        )
        os.chdir(self.test_path)
        parser = SetupCommandLineParser(
            ["2UUY_rec.pdb", "2UUY_lig.pdb", "-s 5", "-g 10", "-anm", "--noxt"]
        )

        create_setup_file(parser.args)

        assert compare_two_files(
            self.test_path / "setup.json", self.golden_data_path / "setup.json",
            ignore=["setup_version", "start_time"]
        )

    def test_prepare_results_environment(self):
        os.chdir(self.test_path)

        prepare_results_environment(20)

        directories = glob(str(self.test_path / "swarm_*/"))

        assert len(directories) == 20

    @raises(LightDockError)
    def test_prepare_results_environment_with_existing_data(self):
        os.chdir(self.test_path)

        prepare_results_environment(20)

        prepare_results_environment()

        assert False

    def test_get_pdb_files(self):
        os.chdir(self.golden_data_path)
        file_names = get_pdb_files(self.golden_data_path / "pdb_files.list")

        expected = ["2UUY_rec.pdb", "2UUY_lig.pdb"]

        assert file_names == expected

    def test_read_input_structure(self):
        os.chdir(self.golden_data_path)

        structure = read_input_structure(
            "2UUY_lig.pdb",
            ignore_oxt=True,
            ignore_hydrogens=False,
            verbose_parser=False,
        )

        assert len(structure.atoms) == 415

    def test_read_multiple_input_structure(self):
        os.chdir(self.golden_data_path)

        structure = read_input_structure(
            "pdb_files.list",
            ignore_oxt=True,
            ignore_hydrogens=False,
            verbose_parser=False,
        )

        assert structure.num_structures == 2

    def test_load_starting_positions(self):
        working_path = self.golden_data_path / "load_starting_positions" / "ok"
        os.chdir(working_path)

        swarms = 2
        glowworms = 10
        use_anm = False
        positions = load_starting_positions(swarms, glowworms, use_anm)

        assert positions == [
            "init/initial_positions_0.dat",
            "init/initial_positions_1.dat",
        ]

    @raises(LightDockError)
    def test_load_starting_positions_wrong_num_swarms(self):
        working_path = self.golden_data_path / "load_starting_positions" / "ok"
        os.chdir(working_path)

        swarms = 5
        glowworms = 10
        use_anm = False
        _ = load_starting_positions(swarms, glowworms, use_anm)

        assert False

    @raises(LightDockError)
    def test_load_starting_positions_wrong_dat_files(self):
        working_path = self.golden_data_path / "load_starting_positions" / "wrong_dat"
        os.chdir(working_path)

        swarms = 2
        glowworms = 10
        use_anm = False
        _ = load_starting_positions(swarms, glowworms, use_anm)

        assert False

    def test_create_simulation_info_file(self):
        setup = get_setup_from_file(
            self.golden_data_path / "create_simulation_info_file" / "setup.json"
        )

        class Args:
            pass

        args = Args()
        for k, v in setup.items():
            setattr(args, k, v)

        file_name = create_simulation_info_file(args, path=self.test_path)

        assert Path(file_name).name == "lightdock.info"

        file_name = create_simulation_info_file(args, path=self.test_path)

        assert Path(file_name).name == "lightdock.info.1"

    def test_check_starting_file(self):
        file_name = (
            self.golden_data_path
            / "load_starting_positions"
            / "ok"
            / "init"
            / "initial_positions_0.dat"
        )
        glowworms = 10
        use_anm = False
        anm_rec = 0
        anm_lig = 0
        assert check_starting_file(file_name, glowworms, use_anm, anm_rec, anm_lig)

        glowworms = 9
        assert not check_starting_file(file_name, glowworms, use_anm, anm_rec, anm_lig)


    def test_simulation_parser(self):
        shutil.copyfile(
            self.golden_data_path / "setup.json", self.test_path / "setup.json"
        )

        os.chdir(self.test_path)
        parser = CommandLineParser(
            ["setup.json", "10", "-c 1", "-l 0", "-s fastdfire"]
        )

        file_name = create_simulation_info_file(parser.args, path=self.test_path)

        assert Path(file_name).name == "lightdock.info"
