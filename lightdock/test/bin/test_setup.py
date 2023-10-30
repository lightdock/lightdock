"""Regression tests for testing simulation setup"""

import shutil
import os
import filecmp
from pathlib import Path
from lightdock.test.support import compare_two_files


class TestSetupWithoutRestraints:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_setup"

    def test_lightdock_setup_automatic(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(self.golden_data_path / "2UUY_rec.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "2UUY_lig.pdb", tmp_path)

        num_glowworms = 25

        command = f"lgd_setup.py 2UUY_rec.pdb 2UUY_lig.pdb -g {num_glowworms} -anm --noxt --noh"
        command += ">> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            tmp_path / "init" / "swarm_centers.pdb",
        )
        assert compare_two_files(
            tmp_path / "setup.json", self.golden_data_path / "setup.json",
            ignore=["setup_version", "start_time"]
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_0.dat",
            tmp_path / "init" / "initial_positions_0.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_45.dat",
            tmp_path / "init" / "initial_positions_45.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_2UUY_rec.pdb",
            tmp_path / "lightdock_2UUY_rec.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_2UUY_lig.pdb",
            tmp_path / "lightdock_2UUY_lig.pdb",
        )
