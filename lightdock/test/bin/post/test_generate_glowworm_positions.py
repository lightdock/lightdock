"""Test for generate_glowworm_positions post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestGenerateGlowwormPositions(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_generate_glowworm_positions/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/generate_glowworm_positions/'

    def teardown(self):
        self.clean_test_path()

    def test_generate_conformations(self):
        os.chdir(self.test_path)
        shutil.copyfile(os.path.join(self.golden_data_path, 'gso_10.out'),
                        os.path.join(self.test_path, 'gso_10.out'))
        command = "lgd_generate_glowworm_positions.py %s > test.out" % (os.path.join(self.test_path, 'gso_10.out'))
        os.system(command)

        assert filecmp.cmp(os.path.join(self.golden_data_path, 'gso_10.pdb'),
                           os.path.join(self.test_path, 'gso_10.pdb'))
