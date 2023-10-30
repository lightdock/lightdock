"""Test for lgd_top post script"""

import os
import filecmp
import shutil
from pathlib import Path


class TestGenerateTop:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "generate_top"

    def test_generate_trajectory(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copyfile(
            self.golden_data_path / "rank_by_scoring.list",
            tmp_path / "rank_by_scoring.list",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_4IZ7_A_noh.pdb",
            tmp_path / "lightdock_4IZ7_A_noh.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_4IZ7_B_noh.pdb",
            tmp_path / "lightdock_4IZ7_B_noh.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "4IZ7_A_noh.pdb", tmp_path / "4IZ7_A_noh.pdb"
        )
        shutil.copyfile(
            self.golden_data_path / "4IZ7_B_noh.pdb", tmp_path / "4IZ7_B_noh.pdb"
        )
        shutil.copyfile(
            self.golden_data_path / "setup.json", tmp_path / "setup.json"
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_rec.nm.npy",
            tmp_path / "lightdock_rec.nm.npy",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_lig.nm.npy",
            tmp_path / "lightdock_lig.nm.npy",
        )

        command = "lgd_top.py 4IZ7_A_noh.pdb 4IZ7_B_noh.pdb rank_by_scoring.list 10 --setup setup.json > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "top_1.pdb", tmp_path / "top_1.pdb"
        )
        assert filecmp.cmp(
            self.golden_data_path / "top_10.pdb", tmp_path / "top_10.pdb"
        )
