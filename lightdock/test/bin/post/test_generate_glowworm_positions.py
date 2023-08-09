"""Test for generate_glowworm_positions post script"""

import os
import filecmp
import shutil
from pathlib import Path


class TestGenerateGlowwormPositions:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "generate_glowworm_positions"

    def test_generate_conformations(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copyfile(
            self.golden_data_path / "gso_10.out", tmp_path / "gso_10.out"
        )

        command = f"lgd_generate_glowworm_positions.py {tmp_path / 'gso_10.out'} > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "gso_10.pdb", tmp_path / "gso_10.pdb"
        )
