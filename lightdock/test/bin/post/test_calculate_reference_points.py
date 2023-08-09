"""Test for lgd_calculate_reference_points post script"""

import os
import filecmp
import shutil
from pathlib import Path


class TestCalculateReferencePoints:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "reference_points"

    def test_calculate_4IZ7_B(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copyfile(
            self.golden_data_path / "4IZ7_B_noh.pdb", tmp_path / "4IZ7_B_noh.pdb"
        )

        command = f"lgd_calculate_reference_points.py {tmp_path / '4IZ7_B_noh.pdb'} > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "4IZ7_B_noh.pdb.vol",
            tmp_path / "4IZ7_B_noh.pdb.vol",
        )
        assert filecmp.cmp(
            self.golden_data_path / "4IZ7_B_noh.pdb.vol.pdb",
            tmp_path / "4IZ7_B_noh.pdb.vol.pdb",
        )
