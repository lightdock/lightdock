"""Test for calculate_success_rate post script"""

import os
import filecmp
import shutil
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestGenerateConformations(RegressionTest):

    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / 'scratch_calculate_success_rate'
        self.golden_data_path = self.path / 'golden_data' / '1PPE' / 'results' / 'clustered'

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_calculate_success_rate(self):
        os.chdir(self.test_path)
        shutil.copyfile(self.golden_data_path / 'rank_by_scoring.list', self.test_path / 'rank_by_scoring.list')

        command = f"lgd_success_rate.py {self.test_path / 'rank_by_scoring.list'} > test.out"
        os.system(command)

        assert filecmp.cmp(self.golden_data_path / 'success_rate.out', self.test_path / 'test.out')
