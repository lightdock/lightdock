"""Tests for PDBReader module"""

import shutil
import os
import filecmp
from nose.tools import assert_almost_equals
from nose.tools import raises
from lightdock.pdbutil.PDBIO import read_atom_line, parse_complex_from_file, write_pdb_to_file
from lightdock.structure.complex import Complex
from lightdock.error.lightdock_errors import PDBParsingError


class TestPDBReader:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch/'
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
    
    @raises(PDBParsingError)
    def test_read_wrong_atom_line(self):
        line = 'ATOM     11  NH2 ARG A   1       2.559  16.752     AA   1.00 14.90           N'
        read_atom_line(line, "ATOM")
        assert False
        
    @raises(PDBParsingError)
    def test_read_nan_in_atom_line(self):
        line = 'ATOM     11  NH2 ARG A   1       2.559  16.752     NaN  1.00 14.90           N'
        read_atom_line(line, "ATOM")
        assert False
 
    @raises(PDBParsingError)
    def test_wrong_atom_number(self):
        line = 'ATOM     OO  NH2 ARG A   1       2.559  16.752    1.00  1.00 14.90           N'
        read_atom_line(line, "ATOM")
        assert False
        
    @raises(PDBParsingError)
    def test_wrong_residue_number(self):
        line = 'ATOM     12  NH2 ARG A   A       2.559  16.752    1.00  1.00 14.90           N'
        read_atom_line(line, "ATOM")
        assert False
        
    def test_default_occupancy(self):
        line = 'ATOM     12  NH2 ARG A   1       2.559  16.752    1.00       14.90           N'
        atom = read_atom_line(line, "ATOM")
        
        assert_almost_equals(1.0, atom.occupancy)
        
    def test_default_b_factor(self):
        line = 'ATOM     12  NH2 ARG A   1       2.559  16.752    1.00  1.00                 N'
        atom = read_atom_line(line, "ATOM")
        
        assert_almost_equals(0.0, atom.b_factor)
 
    def test_default_no_line_type(self):
        line = 'ATOM     12  NH2 ARG A   1       2.559  16.752    1.00  1.00                 N'
        atom = read_atom_line(line)
        
        assert "Atom" == atom.__class__.__name__
    
    def test_parse_complex_from_file(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPE_l_u.pdb')
        assert 224 == len(atoms)
        assert 31 == len(residues)
        assert 2 == len(chains)
        assert 29 == len(chains[0].residues)
        assert 2 == len(chains[1].residues)
        
    def test_parse_multi_model_from_file(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + 'multi_model.pdb')
        assert 11 == len(atoms)
        assert 1 == len(residues)
        assert 1 == len(chains)

    def test_write_pdb_to_file(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPE_l_u.pdb')
        protein = Complex(chains)
        assert atoms == protein.atoms
        
        write_pdb_to_file(protein, self.test_path + '1PPE_l_u.pdb.parsed', protein.atom_coordinates[0])
        
        assert filecmp.cmp(self.golden_data_path + '1PPE_l_u.pdb.parsed', self.test_path + '1PPE_l_u.pdb.parsed')

    def test_parse_pdb_noh(self):
        atoms_to_ignore = ['H']
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPE_lig.pdb.H', atoms_to_ignore)
        protein = Complex(chains)
        assert atoms == protein.atoms
        
        write_pdb_to_file(protein, self.test_path + 'parsed_1PPE_lig.pdb.H', protein.atom_coordinates[0])
        
        assert filecmp.cmp(self.golden_data_path + 'parsed_1PPE_lig.pdb.H', self.test_path + 'parsed_1PPE_lig.pdb.H')
