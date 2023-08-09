"""Test for lgd_create_membrane post script"""

import os
import filecmp
import shutil
from pathlib import Path


class TestCreateMembrane:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "create_membrane"

    def test_rank_with_clusters(self, tmp_path):
        anchor_residue = "C.TRP.3087"

        # Prepare folder structure for this test
        os.chdir(tmp_path)

        shutil.copyfile(
            self.golden_data_path / "receptor.pdb", tmp_path / "receptor.pdb"
        )

        command = f"lgd_create_membrane.py receptor.pdb {anchor_residue} > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "membrane.pdb", tmp_path / "membrane.pdb"
        )
