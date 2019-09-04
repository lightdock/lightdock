"""Test for lgd_top post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestGenerateTop(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_generate_top/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/generate_top/'

    def teardown(self):
        self.clean_test_path()

    def test_generate_trajectory(self):

        # Prepare folder structure for this test
        os.chdir(self.test_path)
        shutil.copyfile(os.path.join(self.golden_data_path, 'rank_by_scoring.list'),
                        os.path.join(self.test_path, 'rank_by_scoring.list'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'lightdock_4IZ7_A_noh.pdb'),
                        os.path.join(self.test_path, 'lightdock_4IZ7_A_noh.pdb'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'lightdock_4IZ7_B_noh.pdb'),
                        os.path.join(self.test_path, 'lightdock_4IZ7_B_noh.pdb'))
        shutil.copyfile(os.path.join(self.golden_data_path, '4IZ7_A_noh.pdb'),
                        os.path.join(self.test_path, '4IZ7_A_noh.pdb'))
        shutil.copyfile(os.path.join(self.golden_data_path, '4IZ7_B_noh.pdb'),
                        os.path.join(self.test_path, '4IZ7_B_noh.pdb'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'setup.json'),
                        os.path.join(self.test_path, 'setup.json'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'lightdock_rec.nm.npy'),
                        os.path.join(self.test_path, 'lightdock_rec.nm.npy'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'lightdock_lig.nm.npy'),
                        os.path.join(self.test_path, 'lightdock_lig.nm.npy'))

        command = "lgd_top.py 4IZ7_A_noh.pdb 4IZ7_B_noh.pdb rank_by_scoring.list 10 --setup setup.json > test.out"
        os.system(command)

        assert filecmp.cmp(os.path.join(self.golden_data_path, 'top_1.pdb'), 
                           os.path.join(self.test_path, 'top_1.pdb'))
        assert filecmp.cmp(os.path.join(self.golden_data_path, 'top_10.pdb'), 
                           os.path.join(self.test_path, 'top_10.pdb'))
 