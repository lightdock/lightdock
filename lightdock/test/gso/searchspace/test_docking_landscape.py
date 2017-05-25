"""Tests for DockingLandscapePosition class"""

import os
from nose.tools import assert_almost_equal
import numpy as np
from lightdock.gso.searchspace.landscape import DockingLandscapePosition
from lightdock.gso.coordinates import Coordinates
from lightdock.scoring.mj3h.driver import MJ3h, MJ3hAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.scoring.tobi.driver import TOBIAdapter, TOBI


class TestDockingLandscapePosition:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPErec.pdb')
        self.receptor = Complex(chains, atoms)

    def tearDown(self):
        pass
    
    def test_evaluate_objective_function_no_movement(self):
        """The result of this test must be the same value as testing the MJ3h function.
        Translation is 0 and Quaternion [1,0,0,0] means no rotation.
        """
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(self.receptor, ligand)
        scoring_function = MJ3h()
        coordinates = Coordinates([0.,0.,0.,1.,0.,0.,0.])
        landscape_position = DockingLandscapePosition(scoring_function, coordinates, 
                                                      adapter.receptor_model, 
                                                      adapter.ligand_model)
        
        assert_almost_equal(2.02, landscape_position.evaluate_objective_function())

    def test_evaluate_objective_function_rotation_y_axis_180(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(self.receptor, ligand)
        scoring_function = MJ3h()
        coordinates = Coordinates([0.,0.,0.,0.,0.,1.,0.])
        landscape_position = DockingLandscapePosition(scoring_function, coordinates, 
                                                      adapter.receptor_model, 
                                                      adapter.ligand_model)
        
        assert_almost_equal(-1.4, landscape_position.evaluate_objective_function())
        
    def test_evaluate_objective_function_rotation_y_axis_180_translation_10(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(self.receptor, ligand)
        scoring_function = MJ3h()
        coordinates = Coordinates([10.,0.,0.,0.,0.,1.,0.])
        landscape_position = DockingLandscapePosition(scoring_function, coordinates, 
                                                      adapter.receptor_model, 
                                                      adapter.ligand_model)
        
        assert_almost_equal(6.39, landscape_position.evaluate_objective_function())

    def test_distance2_same_landscape_position(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(self.receptor, ligand)
        scoring_function = MJ3h()
        coordinates = Coordinates([0., 0., 0., 1., 0., 0., 0.])
        landscape_position1 = DockingLandscapePosition(scoring_function, coordinates, 
                                                      adapter.receptor_model, 
                                                      adapter.ligand_model)
        landscape_position2 = DockingLandscapePosition(scoring_function, coordinates, 
                                                      adapter.receptor_model, 
                                                      adapter.ligand_model)
        assert_almost_equal(0.0, landscape_position1.distance2(landscape_position2))
        
    def test_distance2_10A_translation_x(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(self.receptor, ligand)
        scoring_function = MJ3h()
        coordinates = Coordinates([0., 0., 0., 1., 0., 0., 0.])
        landscape_position1 = DockingLandscapePosition(scoring_function, coordinates,
                                                       adapter.receptor_model, adapter.ligand_model)
        adapter2 = MJ3hAdapter(self.receptor, ligand)
        adapter2.ligand_model.translate([10.0, 0.0, 0.0])
        landscape_position2 = DockingLandscapePosition(scoring_function, coordinates,
                                                       adapter2.receptor_model, adapter2.ligand_model)
        assert_almost_equal(100.0, landscape_position1.distance2(landscape_position2))
        assert_almost_equal(10.0, landscape_position1.distance(landscape_position2))

    def test_distance2_minus_10A_translation_y(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(self.receptor, ligand)
        scoring_function = MJ3h()
        coordinates = Coordinates([0.,0.,0.,1.,0.,0.,0.])
        landscape_position1 = DockingLandscapePosition(scoring_function, coordinates, 
                                                      adapter.receptor_model, 
                                                      adapter.ligand_model)
        adapter2 = MJ3hAdapter(self.receptor, ligand)
        adapter2.ligand_model.translate([0.0,-10.0,0.0])
        landscape_position2 = DockingLandscapePosition(scoring_function, coordinates, 
                                                      adapter2.receptor_model, 
                                                      adapter2.ligand_model)
        assert_almost_equal(100.0, landscape_position1.distance2(landscape_position2))
        assert_almost_equal(10.0, landscape_position1.distance(landscape_position2))
        
    def test_move_step_rot_full_step_trans_half(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = TOBIAdapter(self.receptor, ligand)
        scoring_function = TOBI()
        coordinates1 = Coordinates([0.,0.,0.,1.,0.,0.,0.])
        landscape_position1 = DockingLandscapePosition(scoring_function, coordinates1, 
                                                      adapter.receptor_model, 
                                                      adapter.ligand_model, step_translation=5.0)
        adapter2 = TOBIAdapter(self.receptor, ligand)
        coordinates2 = Coordinates([10.,0.,0.,1.,0.,0.,0.])
        landscape_position2 = DockingLandscapePosition(scoring_function, coordinates2, 
                                                      adapter2.receptor_model, 
                                                      adapter2.ligand_model)
        
        landscape_position1.move(landscape_position2)

        expected_translation = np.array([5.,0.,0.])
        expected_rotation = Quaternion()
        
        assert (expected_translation == landscape_position1.translation).all()
        assert expected_rotation == landscape_position1.rotation

    def test_move_step_rot_full_step_trans_full(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(self.receptor, ligand)
        scoring_function = MJ3h()
        coordinates1 = Coordinates([0.,0.,0.,1.,0.,0.,0.])
        landscape_position1 = DockingLandscapePosition(scoring_function, coordinates1, 
                                                      adapter.receptor_model, 
                                                      adapter.ligand_model, step_translation=10.0)
        adapter2 = MJ3hAdapter(self.receptor, ligand)
        coordinates2 = Coordinates([10.,0.,0.,1.,0.,0.,0.])
        landscape_position2 = DockingLandscapePosition(scoring_function, coordinates2, 
                                                      adapter2.receptor_model, 
                                                      adapter2.ligand_model)
        
        landscape_position1.move(landscape_position2)
        
        assert landscape_position1 == landscape_position2
    
    def test_move_step_rot_half_step_trans_half(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(self.receptor, ligand)
        scoring_function = MJ3h()
        coordinates1 = Coordinates([0.,0.,0.,1.,0.,0.,0.])
        landscape_position1 = DockingLandscapePosition(scoring_function, coordinates1, 
                                                      adapter.receptor_model, 
                                                      adapter.ligand_model, 
                                                      step_translation=5.0,
                                                      step_rotation=0.5)
        adapter2 = MJ3hAdapter(self.receptor, ligand)
        coordinates2 = Coordinates([10.,0.,0.,0.,0.,1.,0.])
        landscape_position2 = DockingLandscapePosition(scoring_function, coordinates2, 
                                                      adapter2.receptor_model, 
                                                      adapter2.ligand_model)
        
        landscape_position1.move(landscape_position2)
        
        expected_translation = np.array([5.,0.,0.])
        expected_rotation = Quaternion(0.707106781, 0.0, 0.707106781, 0.0)
        
        assert (expected_translation == landscape_position1.translation).all()
        assert expected_rotation == landscape_position1.rotation
