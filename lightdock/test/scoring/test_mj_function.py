"""Tests for MJ3h potentials module"""

import os
from nose.tools import assert_almost_equal
from lightdock.scoring.mj3h.driver import MJPotential, MJ3h, MJ3hAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestMJPotential:

    def setup(self):
        pass

    def teardown(self):
        pass
    
    def test_create_MJPotential_interface(self):
        mj = MJPotential()
        
        assert 4 == len(mj.potentials)
    
    def test_create_MJ3h(self):
        mj3h = MJ3h()
        
        assert 20 == len(mj3h.potentials)
        
        assert_almost_equal(-0.84, mj3h.potentials[0][0])
        assert_almost_equal(0.05, mj3h.potentials[15][3])
        assert_almost_equal(0.76, mj3h.potentials[19][19])


class TestMJ3hAdapter:
    
    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        
    def teardown(self):
        pass
    
    def test_create_adapter(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPErec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(receptor, ligand)
        
        assert (223, 3) == adapter.receptor_model.coordinates[0].coordinates.shape
        assert (29, 3) == adapter.ligand_model.coordinates[0].coordinates.shape
        assert 223 == len(adapter.receptor_model.objects)
        assert 29 == len(adapter.ligand_model.objects)

        
class TestMJ3h:
    
    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        # Test for creating object in the interface tests for speeding up purposes
        self.mj3h = MJ3h()

    def teardown(self):
        pass
    
    def test_calculate_MJ3h_1PPE(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPErec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(receptor, ligand)
        assert_almost_equal(2.02, self.mj3h(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                              adapter.ligand_model, adapter.ligand_model.coordinates[0]))
        # Original potential without min cutoff
        # assert_almost_equal(-17.94/2, self.mj3h(receptor, receptor.residue_coordinates, ligand, ligand.residue_coordinates))
    
    def test_calculate_MJ3h_1EAW(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1EAWrec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1EAWlig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(receptor, ligand)
        assert_almost_equal(-4.22, self.mj3h(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                              adapter.ligand_model, adapter.ligand_model.coordinates[0]))
        # Original potential without min cutoff
        # assert_almost_equal(-3.14/2, self.mj3h(receptor, receptor.residue_coordinates, ligand, ligand.residue_coordinates))
 
    def test_calculate_MJ3h_1AY7(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1AY7rec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1AY7lig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(receptor, ligand)
        assert_almost_equal(-12.92, self.mj3h(adapter.receptor_model, adapter.receptor_model.coordinates[0],
                                              adapter.ligand_model, adapter.ligand_model.coordinates[0]))
        # Original potential without min cutoff
        # assert_almost_equal(9.06/2, self.mj3h(receptor, receptor.residue_coordinates, ligand, ligand.residue_coordinates))
