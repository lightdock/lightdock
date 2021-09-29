"""Test for lgd_create_membrane post script"""

import os
import filecmp
import shutil
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestCreateMembrane(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_lgd_create_membrane"
        self.golden_data_path = self.path / "golden_data" / "create_membrane"

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_rank_with_clusters(self):
        anchor_residue = "C.TRP.3087"

        # Prepare folder structure for this test
        os.chdir(self.test_path)

        shutil.copyfile(
            self.golden_data_path / "receptor.pdb", self.test_path / "receptor.pdb"
        )

        command = f"lgd_create_membrane.py receptor.pdb {anchor_residue} > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "membrane.pdb", self.test_path / "membrane.pdb"
        )
