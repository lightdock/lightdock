"""Test for lgd_calculate_reference_points post script"""

import os
import filecmp
import shutil
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestCalculateReferencePoints(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_reference_points"
        self.golden_data_path = self.path / "golden_data" / "reference_points"

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_calculate_4IZ7_B(self):
        os.chdir(self.test_path)
        shutil.copyfile(
            self.golden_data_path / "4IZ7_B_noh.pdb", self.test_path / "4IZ7_B_noh.pdb"
        )

        command = f"lgd_calculate_reference_points.py {self.test_path / '4IZ7_B_noh.pdb'} > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "4IZ7_B_noh.pdb.vol",
            self.test_path / "4IZ7_B_noh.pdb.vol",
        )
        assert filecmp.cmp(
            self.golden_data_path / "4IZ7_B_noh.pdb.vol.pdb",
            self.test_path / "4IZ7_B_noh.pdb.vol.pdb",
        )
