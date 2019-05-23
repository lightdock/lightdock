"""Test for lgd_prepare_new_simulation post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestPrepareNewSimulation(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_prepare_new_simulation/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/prepare_new_simulation/'

    def teardown(self):
        self.clean_test_path()

    def test_prepare_new_simulation(self):

        # Prepare folder structure for this test
        os.chdir(self.test_path)
        swarm_dir = 'swarm_0'
        os.mkdir(swarm_dir)
        shutil.copyfile(os.path.join(self.golden_data_path, swarm_dir, 'gso_10.out'),
                        os.path.join(self.test_path, swarm_dir, 'gso_10.out'))
        shutil.copyfile(os.path.join(self.golden_data_path, swarm_dir, 'cluster.repr'),
                        os.path.join(self.test_path, swarm_dir, 'cluster.repr'))

        
        command = "lgd_prepare_new_simulation.py 0 10 new -nm > test.out"
        os.system(command)

        assert filecmp.cmp(os.path.join(self.golden_data_path, 'new', 'init', 'cluster_centers.pdb'), 
                           os.path.join(self.test_path, 'new', 'init', 'cluster_centers.pdb'))
 