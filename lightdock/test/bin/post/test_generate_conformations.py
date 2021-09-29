"""Test for generate_conformations post script"""

import os
import filecmp
import shutil
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestGenerateConformations(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_generate_conformations"
        self.golden_data_path = self.path / "golden_data" / "generate_conformations"

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_generate_conformations(self):
        os.chdir(self.test_path)
        num_conformations = 2
        shutil.copyfile(
            self.golden_data_path / "1PPE_rec.pdb", self.test_path / "1PPE_rec.pdb"
        )
        shutil.copyfile(
            self.golden_data_path / "1PPE_lig.pdb", self.test_path / "1PPE_lig.pdb"
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_1PPE_rec.pdb",
            self.test_path / "lightdock_1PPE_rec.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_1PPE_lig.pdb",
            self.test_path / "lightdock_1PPE_lig.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "gso_1.out", self.test_path / "gso_1.out"
        )
        command = "lgd_generate_conformations.py %s %s %s %d > test.out" % (
            self.test_path / "1PPE_rec.pdb",
            self.test_path / "1PPE_lig.pdb",
            self.test_path / "gso_1.out",
            num_conformations,
        )
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "lightdock_0.pdb",
            self.test_path / "lightdock_0.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_1.pdb",
            self.test_path / "lightdock_1.pdb",
        )
