"""Tests for PDBReader module"""

import shutil
import os
import filecmp
from pathlib import Path
from nose.tools import assert_almost_equals
from nose.tools import raises
from lightdock.pdbutil.PDBIO import (
    read_atom_line,
    parse_complex_from_file,
    write_pdb_to_file,
)
from lightdock.structure.complex import Complex
from lightdock.error.lightdock_errors import PDBParsingError


class TestPDBReader:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_pdbio"
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

    @raises(PDBParsingError)
    def test_read_wrong_atom_line(self):
        line = "ATOM     11  NH2 ARG A   1       2.559  16.752     AA   1.00 14.90           N"
        read_atom_line(line, "ATOM")
        assert False

    @raises(PDBParsingError)
    def test_read_nan_in_atom_line(self):
        line = "ATOM     11  NH2 ARG A   1       2.559  16.752     NaN  1.00 14.90           N"
        read_atom_line(line, "ATOM")
        assert False

    @raises(PDBParsingError)
    def test_wrong_atom_number(self):
        line = "ATOM     OO  NH2 ARG A   1       2.559  16.752    1.00  1.00 14.90           N"
        read_atom_line(line, "ATOM")
        assert False

    @raises(PDBParsingError)
    def test_wrong_residue_number(self):
        line = "ATOM     12  NH2 ARG A   A       2.559  16.752    1.00  1.00 14.90           N"
        read_atom_line(line, "ATOM")
        assert False

    def test_default_occupancy(self):
        line = "ATOM     12  NH2 ARG A   1       2.559  16.752    1.00       14.90           N"
        atom = read_atom_line(line, "ATOM")

        assert_almost_equals(1.0, atom.occupancy)

    def test_default_b_factor(self):
        line = "ATOM     12  NH2 ARG A   1       2.559  16.752    1.00  1.00                 N"
        atom = read_atom_line(line, "ATOM")

        assert_almost_equals(0.0, atom.b_factor)

    def test_default_no_line_type(self):
        line = "ATOM     12  NH2 ARG A   1       2.559  16.752    1.00  1.00                 N"
        atom = read_atom_line(line)

        assert atom.__class__.__name__ == "Atom"

    def test_parse_complex_from_file(self):
        atoms, residues, chains = parse_complex_from_file(
            self.golden_data_path / "1PPE_l_u.pdb"
        )
        assert len(atoms) == 224
        assert len(residues) == 31
        assert len(chains) == 2
        assert len(chains[0].residues) == 29
        assert len(chains[1].residues) == 2

    def test_parse_multi_model_from_file(self):
        atoms, residues, chains = parse_complex_from_file(
            self.golden_data_path / "multi_model.pdb"
        )
        assert len(atoms) == 11
        assert len(residues) == 1
        assert len(chains) == 1

    def test_write_pdb_to_file(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPE_l_u.pdb"
        )
        protein = Complex(chains)
        assert atoms == protein.atoms

        write_pdb_to_file(
            protein, self.test_path / "1PPE_l_u.pdb.parsed", protein.atom_coordinates[0]
        )

        assert filecmp.cmp(
            self.golden_data_path / "1PPE_l_u.pdb.parsed",
            self.test_path / "1PPE_l_u.pdb.parsed",
        )

    def test_parse_pdb_noh(self):
        atoms_to_ignore = ["H"]
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPE_lig.pdb.H", atoms_to_ignore
        )
        protein = Complex(chains)
        assert atoms == protein.atoms

        write_pdb_to_file(
            protein,
            self.test_path / "parsed_1PPE_lig.pdb.H",
            protein.atom_coordinates[0],
        )

        assert filecmp.cmp(
            self.golden_data_path / "parsed_1PPE_lig.pdb.H",
            self.test_path / "parsed_1PPE_lig.pdb.H",
        )
