"""Tests for ftdock module"""

import os
import shutil
from nose.tools import assert_almost_equal
from lightdock.structure.complex import Complex
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.prep.ftdock import FTDockCoordinatesParser
from lightdock.mathutil.cython.quaternion import Quaternion


class TestFTDockCoordinatesParser:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch/'
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        
        self.quaternions = [Quaternion(-0.122207449, -0.084565307,  0.804585130,  0.574940708),
                            Quaternion( 0.334565328,  0.706772732,  0.229644422, -0.579484069),
                            Quaternion(-0.122207350,  0.541338054,  0.601216794, -0.574940729),
                            Quaternion(-0.239073775,  0.739073825, -0.329056839, -0.536968536),
                            Quaternion(-0.189936901, -0.272159950, -0.302264234,  0.893582267)]
        self.translations = [[-15.734048025344688, -16.449850861746139, 0.64470115806562056],
                             [-28.314891508567385, 12.036561774174739, 10.178799979062202],
                             [-13.844755498611772, -12.788839194052038, 23.000473566520469],
                             [21.7408596820244, -3.0846843852631736, 22.997219312770568],
                             [0.089556216428448465, -0.67458762985527621, -21.185751778087777]]

    def tearDown(self):
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
    
    def test_calculate_starting_points_no_need_to_minimize(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPE_rec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPE_lig.pdb')
        ligand = Complex(chains, atoms)
        
        # Move structures to origin
        receptor_center = receptor.center_of_coordinates()
        ligand_center = ligand.center_of_coordinates()
        receptor_translation = [-1.*c for c in receptor_center]
        ligand_translation = [-1.*c for c in ligand_center]
        receptor.translate(receptor_translation)
        ligand.translate(ligand_translation)
        
        poses = FTDockCoordinatesParser.get_list_of_poses(self.golden_data_path + '1PPE.ftdock',
                                                          ligand_center)
        
        assert 10000 == len(poses)
        
        for i in range(5):
            assert self.quaternions[i] == poses[i].q
            for j in range(3):
                assert_almost_equal(self.translations[i][j], poses[i].translation[j])

