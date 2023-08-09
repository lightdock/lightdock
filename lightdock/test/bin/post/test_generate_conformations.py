"""Test for generate_conformations post script"""

import os
import filecmp
import shutil
from pathlib import Path


class TestGenerateConformations:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "generate_conformations"

    def test_generate_conformations(self, tmp_path):
        os.chdir(tmp_path)
        num_conformations = 2
        shutil.copyfile(
            self.golden_data_path / "1PPE_rec.pdb", tmp_path / "1PPE_rec.pdb"
        )
        shutil.copyfile(
            self.golden_data_path / "1PPE_lig.pdb", tmp_path / "1PPE_lig.pdb"
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_1PPE_rec.pdb",
            tmp_path / "lightdock_1PPE_rec.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_1PPE_lig.pdb",
            tmp_path / "lightdock_1PPE_lig.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "gso_1.out", tmp_path / "gso_1.out"
        )
        command = "lgd_generate_conformations.py %s %s %s %d > test.out" % (
            tmp_path / "1PPE_rec.pdb",
            tmp_path / "1PPE_lig.pdb",
            tmp_path / "gso_1.out",
            num_conformations,
        )
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "lightdock_0.pdb",
            tmp_path / "lightdock_0.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_1.pdb",
            tmp_path / "lightdock_1.pdb",
        )
