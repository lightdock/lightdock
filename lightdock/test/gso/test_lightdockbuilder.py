"""Tests for GSOBuilder class using J1 function"""

import os
import shutil
import filecmp
from lightdock.gso.parameters import GSOParameters
from lightdock.gso.algorithm import LightdockGSOBuilder
from lightdock.gso.boundaries import Boundary, BoundingBox
from lightdock.mathutil.lrandom import MTGenerator
from lightdock.constants import MAX_TRANSLATION, MAX_ROTATION
from lightdock.structure.complex import Complex
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.scoring.mj3h.driver import MJ3h, MJ3hAdapter


class TestLightDockGSOBuilder:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch/'
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        self.gso_parameters = GSOParameters()
        
        self.bounding_box = BoundingBox([Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                                         Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                                         Boundary(-MAX_TRANSLATION, MAX_TRANSLATION),
                                         Boundary(-MAX_ROTATION, MAX_ROTATION),
                                         Boundary(-MAX_ROTATION, MAX_ROTATION),
                                         Boundary(-MAX_ROTATION, MAX_ROTATION),
                                         Boundary(-MAX_ROTATION, MAX_ROTATION)])
        self.random_number_generator = MTGenerator(324324)
        
    def tearDown(self):
        try:
            shutil.rmtree(self.test_path)
        except:
            pass

    def test_LightDockGSOBuilder_using_FromFileInitializer(self):
        builder = LightdockGSOBuilder()
        number_of_glowworms = 5
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPErec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(receptor, ligand)
        scoring_function = MJ3h()
        gso = builder.create_from_file(number_of_glowworms, self.random_number_generator, 
                                       self.gso_parameters, [adapter], [scoring_function], self.bounding_box,
                                       self.golden_data_path+'initial_positions_1PPE.txt',
                                       0.5, 0.5, 0.5, False, 10, 10)
        
        assert 5 == gso.swarm.get_size()
        
        gso.report(self.test_path + 'report.out')
        assert filecmp.cmp(self.golden_data_path + 'report_lightdockbuilder.out',
                           self.test_path + 'report.out')
