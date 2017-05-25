"""Tests for C implementation of DFIRE scoring function module"""

from nose.tools import assert_almost_equal
import os
from lightdock.scoring.fastdfire.driver import DFIRE, DFIREAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
        

class TestFastDFIRE:
    
    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        self.dfire = DFIRE()

    def teardown(self):
        pass
    
    def test_calculate_FastDFIRE_1PPE(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPErec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = DFIREAdapter(receptor, ligand)
        assert_almost_equal(-17.3745706065, self.dfire(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                                       adapter.ligand_model, adapter.ligand_model.coordinates[0]))
    
    def test_calculate_FastDFIRE_1EAW(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1EAWrec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1EAWlig.pdb')
        ligand = Complex(chains, atoms)
        adapter = DFIREAdapter(receptor, ligand)
        assert_almost_equal(-16.2239702546, self.dfire(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                                       adapter.ligand_model, adapter.ligand_model.coordinates[0]))
        
    def test_calculate_FastDFIRE_1AY7(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1AY7rec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1AY7lig.pdb')
        ligand = Complex(chains, atoms)
        adapter = DFIREAdapter(receptor, ligand)
        assert_almost_equal(-20.7459619159, self.dfire(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                                       adapter.ligand_model, adapter.ligand_model.coordinates[0]))
