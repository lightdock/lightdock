"""Tests for parser module"""
import shutil
import os
import argparse
from pathlib import Path
from lightdock.util.parser import (
    get_lightdock_structures,
    valid_file,
    valid_integer_number,
    valid_natural_number,
    valid_float_number,
)
from nose.tools import raises


class TestParserUtils:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_parser_utils"
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

    def test_get_lightdock_structures(self):
        os.chdir(self.golden_data_path)
        # First test list of structures
        list_of_pdbs = self.golden_data_path / "pdbs.list"
        structures = get_lightdock_structures(list_of_pdbs)

        assert structures == [
            "lightdock_1czy_protein.pdb",
            "lightdock_1czy_peptide.pdb",
        ]

        # Test single structure
        single_pdb = self.golden_data_path / "1czy_protein.pdb"
        structures = get_lightdock_structures(single_pdb)

        assert Path(structures[0]).name == "lightdock_1czy_protein.pdb"

    def test_valid_file_ok(self):
        os.chdir(self.golden_data_path)

        filename = valid_file("lightdock_1czy_protein.pdb")

        assert filename == "lightdock_1czy_protein.pdb"

    @raises(argparse.ArgumentTypeError)
    def test_valid_file_ko(self):
        os.chdir(self.golden_data_path)

        valid_file("1czy_protein.pdb")

        assert False

    def test_valid_integer_number_ok(self):
        assert valid_integer_number("1")

    @raises(argparse.ArgumentTypeError)
    def test_valid_integer_number_ko(self):
        assert valid_natural_number("aa") == 0

    @raises(argparse.ArgumentTypeError)
    def test_valid_integer_number_ko_1(self):
        assert not valid_integer_number("aa")

    @raises(argparse.ArgumentTypeError)
    def test_valid_integer_number_ko_2(self):
        assert valid_integer_number("0") == 0

    def test_valid_natural_number_ok(self):
        assert valid_natural_number("1") == 1

    def test_valid_natural_number_ok_2(self):
        assert valid_natural_number("0") == 0

    @raises(argparse.ArgumentTypeError)
    def test_valid_natural_number_ko(self):
        assert valid_natural_number("-1") == 0

    def test_valid_float_number_ok(self):
        assert valid_float_number("1.0") == 1.0

    @raises(argparse.ArgumentTypeError)
    def test_valid_float_number_ko_1(self):
        assert valid_float_number("aa") == 0.0

    @raises(argparse.ArgumentTypeError)
    def test_valid_float_number_ko_2(self):
        assert valid_float_number("-1.0") == -1.0
