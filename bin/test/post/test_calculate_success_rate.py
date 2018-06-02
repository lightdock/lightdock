"""Test for calculate_success_rate post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestGenerateConformations(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_calculate_success_rate/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) \
                                + '/golden_data/1PPE/results/clustered/'

    def teardown(self):
        self.clean_test_path()

    def test_calculate_success_rate(self):
        os.chdir(self.test_path)
        shutil.copyfile(self.golden_data_path + 'rank_by_scoring.list', self.test_path + 'rank_by_scoring.list')

        command = "lgd_success_rate.py %s > test.out" % (self.test_path + 'rank_by_scoring.list')
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'success_rate.out', self.test_path + 'test.out')
