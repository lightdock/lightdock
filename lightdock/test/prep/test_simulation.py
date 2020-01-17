"""Tests for the simulation module"""

import os
import shutil
import filecmp
from glob import glob
from nose.tools import assert_almost_equal, raises
from lightdock.error.lightdock_errors import LightDockError
from lightdock.prep.simulation import parse_restraints_file, get_restraints, get_default_box,\
    get_setup_from_file, create_setup_file, prepare_results_environment, get_pdb_files, \
    read_input_structure
from lightdock.structure.complex import Complex
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.util.parser import SetupCommandLineParser


class TestParsingRestraintsFile:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = os.path.join(self.path, 'scratch_simulation')
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.join(os.path.normpath(os.path.dirname(os.path.realpath(__file__))),
                                             'golden_data')

    def tearDown(self):
        try:
            shutil.rmtree(self.test_path)
        except:
            pass

    def test_simple_restraints_file(self):
        input_file = os.path.join(self.golden_data_path, 'rst_1.lst')

        restraints = parse_restraints_file(input_file)

        expected = {'receptor': {'active': ['A.ALA.1'], 'passive': []}, 'ligand': {'active': [], 'passive': ['B.TYR.1']}}
        
        assert restraints == expected

    def test_restraints_file_with_duplicate(self):
        input_file = os.path.join(self.golden_data_path, 'rst_2.lst')

        restraints = parse_restraints_file(input_file)

        expected = {'receptor': {'active': ['A.ALA.1'], 'passive': []}, 'ligand': {'active': [], 'passive': ['B.TYR.1']}}
        
        assert restraints == expected

    def test_empty_restraints_file(self):
        input_file = os.path.join(self.golden_data_path, 'rst_3.lst')

        restraints = parse_restraints_file(input_file)

        expected = {'receptor': {'active': [], 'passive': []}, 'ligand': {'active': [], 'passive': []}}
        
        assert restraints == expected


class TestRestraints:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = os.path.join(self.path, 'scratch')
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.join(os.path.normpath(os.path.dirname(os.path.realpath(__file__))),
                                             'golden_data')

    def tearDown(self):
        try:
            shutil.rmtree(self.test_path)
        except:
            pass

    def test_get_restraints(self):
        input_file = os.path.join(self.golden_data_path, '2UUY_lig.pdb')
        atoms, residues, chains = parse_complex_from_file(input_file)
        structure = Complex(chains)
        restraints = {'active':['B.ALA.21'], 'passive':['B.GLY.75']}
        
        residues = get_restraints(structure, restraints)

        assert len(residues['active']) == 1 and len(residues['passive']) == 1
        assert residues['active'][0].name == 'ALA' and residues['active'][0].number == 21
        assert residues['passive'][0].name == 'GLY' and residues['passive'][0].number == 75

    @raises(LightDockError)
    def test_get_restraints_with_error(self):
        input_file = os.path.join(self.golden_data_path, '2UUY_lig.pdb')
        atoms, residues, chains = parse_complex_from_file(input_file)
        structure = Complex(chains)
        restraints = {'active':['B.VAL.21'], 'passive':['B.GLY.75']}
        
        residues = get_restraints(structure, restraints)

        assert False


class TestSimulation:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = os.path.join(self.path, 'scratch')
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.join(os.path.normpath(os.path.dirname(os.path.realpath(__file__))),
                                             'golden_data')

    def tearDown(self):
        try:
            shutil.rmtree(self.test_path)
        except:
            pass

    def test_bounding_box(self):
        box = get_default_box(use_anm=False, anm_rec=0, anm_lig=0)

        assert box.dimension == 7

    def test_bounding_box_anm(self):
        box = get_default_box(use_anm=True, anm_rec=5, anm_lig=3)

        assert box.dimension == 7 + 5 + 3

    def test_get_setup_from_file(self):
        read = get_setup_from_file(os.path.join(self.golden_data_path, 'setup.json'))

        expected = {
            "anm_lig": 10,
            "anm_rec": 10,
            "anm_seed": 324324,
            "ftdock_file": None,
            "glowworms": 10,
            "ligand_pdb": "2UUY_lig.pdb",
            "membrane": False,
            "noh": False,
            "noxt": True,
            "receptor_pdb": "2UUY_rec.pdb",
            "restraints": None,
            "starting_points_seed": 324324,
            "swarms": 5,
            "use_anm": True,
            "verbose_parser": False
        }

        assert read == expected

    def test_create_setup_file(self):
        shutil.copyfile(os.path.join(self.golden_data_path, '2UUY_rec.pdb'), 
                        os.path.join(self.test_path, '2UUY_rec.pdb'))
        shutil.copyfile(os.path.join(self.golden_data_path, '2UUY_lig.pdb'), 
                        os.path.join(self.test_path, '2UUY_lig.pdb'))
        os.chdir(self.test_path)
        parser = SetupCommandLineParser(['2UUY_rec.pdb', '2UUY_lig.pdb', '5', '10', '-anm', '--noxt'])

        create_setup_file(parser.args)

        assert filecmp.cmp(os.path.join(self.test_path, 'setup.json'), 
                           os.path.join(self.golden_data_path, 'setup.json'))

    def test_prepare_results_environment(self):
        os.chdir(self.test_path)

        prepare_results_environment(20)

        directories = glob(os.path.join(self.test_path, 'swarm_*/'))

        assert len(directories) == 20

    @raises(LightDockError)
    def test_prepare_results_environment_with_existing_data(self):
        os.chdir(self.test_path)

        prepare_results_environment(20)

        prepare_results_environment()

        assert False

    def test_get_pdb_files(self):
        os.chdir(self.golden_data_path)
        file_names = get_pdb_files(os.path.join(self.golden_data_path, 'pdb_files.list'))

        expected = ['2UUY_rec.pdb', '2UUY_lig.pdb']

        assert file_names == expected

    def test_read_input_structure(self):
        os.chdir(self.golden_data_path)

        structure = read_input_structure('2UUY_lig.pdb', ignore_oxt=True, 
                                         ignore_hydrogens=False, verbose_parser=False)

        assert len(structure.atoms) == 415

    def test_read_multiple_input_structure(self):
        os.chdir(self.golden_data_path)

        structure = read_input_structure('pdb_files.list', ignore_oxt=True, 
                                         ignore_hydrogens=False, verbose_parser=False)

        assert structure.num_structures == 2
