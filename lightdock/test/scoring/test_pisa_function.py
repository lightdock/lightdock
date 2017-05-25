"""Tests for PISA scoring function module"""

from nose.tools import assert_almost_equal
import os
from lightdock.scoring.pisa.driver import PISAPotential, PISA, PISAAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestPISAPotential:

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_create_PISAPotential_interface(self):
        potential = PISAPotential()
        assert True
        

class TestPISA:
    """Original PISA scoring energy goes from negative to positive"""
    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        self.pisa = PISA()

    def teardown(self):
        pass
    
    def test_calculate_PISA_1PPE(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPErec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = PISAAdapter(receptor, ligand)
        assert_almost_equal(-0.4346, round(self.pisa(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                                     adapter.ligand_model, adapter.ligand_model.coordinates[0]), 4))
    
    def test_calculate_PISA_1EAW(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1EAWrec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1EAWlig.pdb')
        ligand = Complex(chains, atoms)
        adapter = PISAAdapter(receptor, ligand)
        assert_almost_equal(-0.2097, round(self.pisa(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                                     adapter.ligand_model, adapter.ligand_model.coordinates[0]), 4))
        
    def test_calculate_PISA_1AY7(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1AY7rec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1AY7lig.pdb')
        ligand = Complex(chains, atoms)
        adapter = PISAAdapter(receptor, ligand)
        assert_almost_equal(-0.2141, round(self.pisa(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                                     adapter.ligand_model, adapter.ligand_model.coordinates[0]), 4))
