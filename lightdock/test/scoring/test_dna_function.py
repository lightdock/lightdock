"""Tests for CPyDockDNA scoring function module"""

from nose.tools import assert_almost_equal
import os
from lightdock.scoring.dna.driver import DNA, DNAAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestPyDockDNA:

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        self.dna = DNA()

    def teardown(self):
        pass

    def test_calculate_DNA_3MFK(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '3mfk_homodimer.pdb')
        receptor = Complex(chains, atoms, structure_file_name=(self.golden_data_path + '3mfk_homodimer.pdb'))
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '3mfk_dna.pdb')
        ligand = Complex(chains, atoms, structure_file_name=(self.golden_data_path + '3mfk_dna.pdb'))
        adapter = DNAAdapter(receptor, ligand)
        assert_almost_equal(-2716.68018700585, self.dna(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                                        adapter.ligand_model, adapter.ligand_model.coordinates[0]))

    def test_calculate_DNA_3MFK_with_hydrogens(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '3mfk_homodimer.pdb.H')
        receptor = Complex(chains, atoms, structure_file_name=(self.golden_data_path + '3mfk_homodimer.pdb.H'))
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '3mfk_dna.pdb')
        ligand = Complex(chains, atoms, structure_file_name=(self.golden_data_path + '3mfk_dna.pdb'))
        adapter = DNAAdapter(receptor, ligand)
        assert_almost_equal(688.1703668834168, self.dna(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                                        adapter.ligand_model, adapter.ligand_model.coordinates[0]))
