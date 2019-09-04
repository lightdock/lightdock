"""Test for lgd_generate_trajectory post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestGenerateTrajectory(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_generate_trajectory/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/generate_trajectory/'

    def teardown(self):
        self.clean_test_path()

    def test_generate_trajectory(self):

        # Prepare folder structure for this test
        os.chdir(self.test_path)
        swarm_dir = 'swarm_0'
        os.mkdir(swarm_dir)
        shutil.copyfile(os.path.join(self.golden_data_path, swarm_dir, 'gso_0.out'),
                        os.path.join(self.test_path, swarm_dir, 'gso_0.out'))
        shutil.copyfile(os.path.join(self.golden_data_path, swarm_dir, 'gso_10.out'),
                        os.path.join(self.test_path, swarm_dir, 'gso_10.out'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'lightdock_4IZ7_A_noh.pdb'),
                        os.path.join(self.test_path, 'lightdock_4IZ7_A_noh.pdb'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'lightdock_4IZ7_B_noh.pdb'),
                        os.path.join(self.test_path, 'lightdock_4IZ7_B_noh.pdb'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'setup.json'),
                        os.path.join(self.test_path, 'setup.json'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'lightdock_rec.nm.npy'),
                        os.path.join(self.test_path, 'lightdock_rec.nm.npy'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'lightdock_lig.nm.npy'),
                        os.path.join(self.test_path, 'lightdock_lig.nm.npy'))

        os.chdir(os.path.join(self.test_path, swarm_dir))
        command = "lgd_generate_trajectory.py 4 10 ../lightdock_4IZ7_A_noh.pdb ../lightdock_4IZ7_B_noh.pdb ../setup.json > test.out"
        os.system(command)

        assert filecmp.cmp(os.path.join(self.golden_data_path, 'swarm_0', 'trajectory_4_step_0.pdb'), 
                           os.path.join(self.test_path, 'swarm_0', 'trajectory_4_step_0.pdb'))
        assert filecmp.cmp(os.path.join(self.golden_data_path, 'swarm_0', 'trajectory_4_step_10.pdb'), 
                           os.path.join(self.test_path, 'swarm_0', 'trajectory_4_step_10.pdb'))
 