"""Tests for the simulation module"""

import os
import shutil
import filecmp
from nose.tools import assert_almost_equal
from lightdock.prep.simulation import parse_restraints_file


class TestParsingRestraintsFile:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_simulation/'
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'

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
