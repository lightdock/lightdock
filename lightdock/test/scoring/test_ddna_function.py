"""Tests for DDNA scoring function module"""

from nose.tools import assert_almost_equal
import os
from lightdock.scoring.ddna.driver import DDNA, DDNAAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestDDNA:

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        self.dna = DDNA()

    def teardown(self):
        pass

    def test_calculate_DNA_1AZP(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1azp_prot.pdb')
        receptor = Complex(chains, atoms, structure_file_name=(self.golden_data_path + '1azp_prot.pdb'))
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1azp_dna.pdb')
        ligand = Complex(chains, atoms, structure_file_name=(self.golden_data_path + '1azp_dna.pdb'))
        adapter = DDNAAdapter(receptor, ligand)
        assert_almost_equal(6.915295143021656, self.dna(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                                        adapter.ligand_model, adapter.ligand_model.coordinates[0]))
