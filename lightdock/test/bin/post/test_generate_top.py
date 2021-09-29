"""Test for lgd_top post script"""

import os
import filecmp
import shutil
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestGenerateTop(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_generate_top"
        self.golden_data_path = self.path / "golden_data" / "generate_top"

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_generate_trajectory(self):
        # Prepare folder structure for this test
        os.chdir(self.test_path)
        shutil.copyfile(
            self.golden_data_path / "rank_by_scoring.list",
            self.test_path / "rank_by_scoring.list",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_4IZ7_A_noh.pdb",
            self.test_path / "lightdock_4IZ7_A_noh.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_4IZ7_B_noh.pdb",
            self.test_path / "lightdock_4IZ7_B_noh.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "4IZ7_A_noh.pdb", self.test_path / "4IZ7_A_noh.pdb"
        )
        shutil.copyfile(
            self.golden_data_path / "4IZ7_B_noh.pdb", self.test_path / "4IZ7_B_noh.pdb"
        )
        shutil.copyfile(
            self.golden_data_path / "setup.json", self.test_path / "setup.json"
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_rec.nm.npy",
            self.test_path / "lightdock_rec.nm.npy",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_lig.nm.npy",
            self.test_path / "lightdock_lig.nm.npy",
        )

        command = "lgd_top.py 4IZ7_A_noh.pdb 4IZ7_B_noh.pdb rank_by_scoring.list 10 --setup setup.json > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "top_1.pdb", self.test_path / "top_1.pdb"
        )
        assert filecmp.cmp(
            self.golden_data_path / "top_10.pdb", self.test_path / "top_10.pdb"
        )
