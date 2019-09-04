"""Test for lgd_calculate_reference_points post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestCalculateReferencePoints(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_reference_points/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) \
                                + '/golden_data/reference_points/'

    def teardown(self):
        self.clean_test_path()

    def test_calculate_4IZ7_B(self):
        os.chdir(self.test_path)
        shutil.copyfile(os.path.join(self.golden_data_path, '4IZ7_B_noh.pdb'),
                        os.path.join(self.test_path, '4IZ7_B_noh.pdb'))

        command = "lgd_calculate_reference_points.py %s > test.out" % (os.path.join(self.test_path, '4IZ7_B_noh.pdb'))
        os.system(command)

        assert filecmp.cmp(os.path.join(self.golden_data_path, '4IZ7_B_noh.pdb.vol'), 
                           os.path.join(self.test_path, '4IZ7_B_noh.pdb.vol'))
        assert filecmp.cmp(os.path.join(self.golden_data_path, '4IZ7_B_noh.pdb.vol.pdb'), 
                           os.path.join(self.test_path, '4IZ7_B_noh.pdb.vol.pdb'))
