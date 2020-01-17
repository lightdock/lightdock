"""Tests for the simulation module"""

import os
import shutil
import filecmp
from nose.tools import assert_almost_equal, raises
from lightdock.error.lightdock_errors import LightDockError
from lightdock.prep.simulation import parse_restraints_file, get_restraints, get_default_box,\
    get_setup_from_file, create_setup_file
from lightdock.structure.complex import Complex
from lightdock.pdbutil.PDBIO import parse_complex_from_file


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
