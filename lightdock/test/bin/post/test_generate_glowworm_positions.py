"""Test for generate_glowworm_positions post script"""

import os
import filecmp
import shutil
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestGenerateGlowwormPositions(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_generate_glowworm_positions"
        self.golden_data_path = (
            self.path / "golden_data" / "generate_glowworm_positions"
        )

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_generate_conformations(self):
        os.chdir(self.test_path)
        shutil.copyfile(
            self.golden_data_path / "gso_10.out", self.test_path / "gso_10.out"
        )

        command = f"lgd_generate_glowworm_positions.py {self.test_path / 'gso_10.out'} > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "gso_10.pdb", self.test_path / "gso_10.pdb"
        )
