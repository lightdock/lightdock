"""Tests for parser module"""
import shutil
import os
import filecmp
import sys
import argparse

from lightdock.util.parser import get_lightdock_structures, valid_file, valid_integer_number, \
    valid_natural_number, valid_float_number
from nose.tools import raises



class TestParserUtils:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_parser_utils/'
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
    
    def test_get_lightdock_structures(self):
        os.chdir(self.golden_data_path)

        # First test list of structures
        list_of_pdbs = os.path.join(self.golden_data_path, 'pdbs.list')
        structures = get_lightdock_structures(list_of_pdbs)

        assert structures == ['lightdock_1czy_protein.pdb', 'lightdock_1czy_peptide.pdb']

        # Test single structure
        single_pdb = os.path.join(self.golden_data_path, '1czy_protein.pdb')
        structures = get_lightdock_structures(single_pdb)

        assert os.path.basename(structures[0]) == 'lightdock_1czy_protein.pdb'

    def test_valid_file_ok(self):
        os.chdir(self.golden_data_path)

        filename = valid_file('lightdock_1czy_protein.pdb')

        assert filename == 'lightdock_1czy_protein.pdb'

    @raises(argparse.ArgumentTypeError)
    def test_valid_file_ko(self):
        os.chdir(self.golden_data_path)

        filename = valid_file('1czy_protein.pdb')

        assert False

    def test_valid_integer_number_ok(self):
        assert 1 == valid_integer_number('1')

    @raises(argparse.ArgumentTypeError)
    def test_valid_integer_number_ko_1(self):
        assert 0 == valid_integer_number('aa')

    @raises(argparse.ArgumentTypeError)
    def test_valid_integer_number_ko_2(self):
        assert 0 == valid_integer_number('0')

    def test_valid_natural_number_ok(self):
        assert 1 == valid_natural_number('1')

    def test_valid_natural_number_ok_2(self):
        assert 0 == valid_natural_number('0')

    @raises(argparse.ArgumentTypeError)
    def test_valid_integer_number_ko(self):
        assert 0 == valid_natural_number('aa')

    def test_valid_float_number_ok(self):
        assert 1.0 == valid_float_number('1.0')

    @raises(argparse.ArgumentTypeError)
    def test_valid_float_number_ko_1(self):
        assert 0.0 == valid_float_number('aa')

    @raises(argparse.ArgumentTypeError)
    def test_valid_float_number_ko_2(self):
        assert -1.0 == valid_float_number('-1.0')
