"""Test for lgd_create_membrane post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestCreateMembrane(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_lgd_create_membrane/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/create_membrane/'

    def teardown(self):
        self.clean_test_path()

    def test_rank_with_clusters(self):
        anchor_residue = "C.TRP.3087"

        # Prepare folder structure for this test
        os.chdir(self.test_path)

        shutil.copyfile(os.path.join(self.golden_data_path, 'receptor.pdb'),
                        os.path.join(self.test_path, 'receptor.pdb'))
        
        command = f'lgd_create_membrane.py receptor.pdb {anchor_residue} > test.out'
        os.system(command)

        assert filecmp.cmp(os.path.join(self.golden_data_path, 'membrane.pdb'), 
                           os.path.join(self.test_path, 'membrane.pdb'))
