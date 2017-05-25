"""Tests for TOBI scoring function module"""

from nose.tools import assert_almost_equal #@UnresolvedImport
import os
from lightdock.scoring.tobi.driver import TOBIPotential, TOBI, TOBIAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestTOBIPotential:

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_create_TOBIPotential_interface(self):
        potential = TOBIPotential()
        
        assert 22 == len(potential.tobi_sc_1)
        assert 22 == len(potential.tobi_sc_2)
        
        assert_almost_equal(-0.59, potential.tobi_sc_1[0][0])
        assert_almost_equal(-0.09, potential.tobi_sc_1[-1][-1])
        assert_almost_equal(1.37, potential.tobi_sc_1[1][20])
        
        assert_almost_equal(-0.58, potential.tobi_sc_2[0][0])
        assert_almost_equal(-0.24, potential.tobi_sc_2[-1][-1])
        assert_almost_equal(0.39, potential.tobi_sc_2[3][20])
        

class TestTOBI:
    
    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        self.tobi = TOBI()

    def tearDown(self):
        pass
    
    def test_calculate_TOBI_1PPE(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPErec.pdb') #@UnusedVariable
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb') #@UnusedVariable
        ligand = Complex(chains, atoms)
        adapter = TOBIAdapter(receptor, ligand)
        assert_almost_equal(17.58, self.tobi(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                             adapter.ligand_model, adapter.ligand_model.coordinates[0]))
        # Fixed to include vdw radii in search
        # assert_almost_equal(-150.42, self.tobi(adapter.receptor_model, adapter.ligand_model))
    
    def test_calculate_TOBI_1EAW(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1EAWrec.pdb') #@UnusedVariable
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1EAWlig.pdb') #@UnusedVariable
        ligand = Complex(chains, atoms)
        adapter = TOBIAdapter(receptor, ligand)
        assert_almost_equal(-11.22, self.tobi(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                              adapter.ligand_model, adapter.ligand_model.coordinates[0]))
        # Fixed to include vdw radii in search
        # assert_almost_equal(-133.22, self.tobi(adapter.receptor_model, adapter.ligand_model))
        
    def test_calculate_TOBI_1AY7(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1AY7rec.pdb') #@UnusedVariable
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1AY7lig.pdb') #@UnusedVariable
        ligand = Complex(chains, atoms)
        adapter = TOBIAdapter(receptor, ligand)
        assert_almost_equal(2.34, self.tobi(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                            adapter.ligand_model, adapter.ligand_model.coordinates[0]))
        # Fixed to include vdw radii in search
        # assert_almost_equal(-136.66, self.tobi(adapter.receptor_model, adapter.ligand_model))
