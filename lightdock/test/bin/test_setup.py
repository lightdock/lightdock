"""Regression tests for testing simulation setup"""

import shutil
import os
import filecmp
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest
from lightdock.test.support import compare_two_files


class TestSetupWithoutRestraints(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_setup_no_restraints"
        self.golden_data_path = self.path / "golden_data" / "regression_setup"

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "2UUY_rec.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "2UUY_lig.pdb", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_setup_automatic(self):
        os.chdir(self.test_path)

        num_glowworms = 25

        command = f"lightdock3_setup.py 2UUY_rec.pdb 2UUY_lig.pdb -g {num_glowworms} -anm --noxt --noh"
        command += ">> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            self.test_path / "init" / "swarm_centers.pdb",
        )
        assert compare_two_files(
            self.test_path / "setup.json", self.golden_data_path / "setup.json",
            ignore=["setup_version", "start_time"]
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_0.dat",
            self.test_path / "init" / "initial_positions_0.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_45.dat",
            self.test_path / "init" / "initial_positions_45.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_2UUY_rec.pdb",
            self.test_path / "lightdock_2UUY_rec.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_2UUY_lig.pdb",
            self.test_path / "lightdock_2UUY_lig.pdb",
        )
