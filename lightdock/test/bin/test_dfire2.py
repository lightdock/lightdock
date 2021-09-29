"""Regression tests for testing DFIRE2 scoring function"""

import shutil
import os
import filecmp
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestRegressionDFIRE2Short(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_dfire2_short"
        self.golden_data_path = self.path / "golden_data" / "regression_dfire2_short"

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "2UUY_rec.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "2UUY_lig.pdb", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_2uuy_5_steps_25_glowworms_100_swarms(self):
        os.chdir(self.test_path)
        num_swarms = 100
        num_glowworms = 25
        steps = 5

        command = f"lightdock3_setup.py 2UUY_rec.pdb 2UUY_lig.pdb -g {num_glowworms} -s {num_swarms} -anm"
        command += ">> test_lightdock.out"
        os.system(command)

        command = f"lightdock3.py -c 1 -s dfire2 setup.json {steps} -l 10 >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_10" / "gso_0.out",
            self.test_path / "swarm_10" / "gso_0.out",
        )
        assert (self.test_path / "swarm_10" / "gso_5.out").exists()


class TestRegressionDFIRE2Restraints(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_dfire2_restraints"
        self.golden_data_path = (
            self.path / "golden_data" / "regression_dfire2_restraints"
        )

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "2UUY_rec.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "2UUY_lig.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "restraints.list", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_2uuy_5_steps_25_glowworms_rst(self):
        os.chdir(self.test_path)
        num_glowworms = 25
        steps = 5

        command = (
            f"lightdock3_setup.py 2UUY_rec.pdb 2UUY_lig.pdb -g {num_glowworms} -anm "
        )
        command += "-rst restraints.list >> test_lightdock.out"
        os.system(command)

        command = f"lightdock3.py -c 1 -s dfire2 setup.json {steps} -l 0 >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_0.out",
            self.test_path / "swarm_0" / "gso_0.out",
        )
        assert (self.test_path / "swarm_0" / "gso_5.out").exists()


class TestRegressionDFIRE2Long(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_dfire2_long"
        self.golden_data_path = self.path / "golden_data" / "regression_dfire2_long"

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "2UUY_rec.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "2UUY_lig.pdb", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_2uuy_40_steps_50_glowworms(self):
        os.chdir(self.test_path)
        num_glowworms = 50
        steps = 40

        command = f"lightdock3_setup.py 2UUY_rec.pdb 2UUY_lig.pdb -g {num_glowworms} >> test_lightdock.out"
        os.system(command)

        command = f"lightdock3.py -c 1 -s dfire2 setup.json {steps} -l 100 >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_0.out",
            self.test_path / "swarm_100" / "gso_0.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            self.test_path / "init" / "swarm_centers.pdb",
        )
        assert (self.test_path / "swarm_100" / "gso_10.out").exists()
        assert (self.test_path / "swarm_100" / "gso_20.out").exists()
        assert (self.test_path / "swarm_100" / "gso_30.out").exists()
        assert (self.test_path / "swarm_100" / "gso_40.out").exists()
