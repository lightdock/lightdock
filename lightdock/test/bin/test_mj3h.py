"""Regression tests for testing Mj3h scoring function"""

import shutil
import os
import filecmp
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestRegressionMj3hShort(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_mj3h_short"
        self.golden_data_path = self.path / "golden_data" / "regression_mj3h_short"

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "2UUY_rec.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "2UUY_lig.pdb", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_2uuy_10_steps_25_glowworms_100_swarms(self):
        os.chdir(self.test_path)
        num_swarms = 100
        num_glowworms = 25
        steps = 10

        command = f"lightdock3_setup.py 2UUY_rec.pdb 2UUY_lig.pdb -g {num_glowworms} -s {num_swarms} "
        command += ">> test_lightdock.out"
        os.system(command)

        command = f"lightdock3.py -c 1 -f {self.golden_data_path / 'glowworm.conf'} "
        command += f"-s mj3h setup.json {steps} -l 0 >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_0.out",
            self.test_path / "swarm_0" / "gso_0.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_10.out",
            self.test_path / "swarm_0" / "gso_10.out",
        )


class TestRegressionMj3hLong(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_mj3h_long"
        self.golden_data_path = self.path / "golden_data" / "regression_mj3h_long"

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "1PPE_rec.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "1PPE_lig.pdb", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_1ppe_50_steps_20_glowworms_100_swarms(self):
        os.chdir(self.test_path)
        num_swarms = 100
        num_glowworms = 20
        steps = 30

        command = f"lightdock3_setup.py 1PPE_rec.pdb 1PPE_lig.pdb -g {num_glowworms} -s {num_swarms} "
        command += ">> test_lightdock.out"
        os.system(command)

        command = f"lightdock3.py -c 1 -f {self.golden_data_path / 'glowworm.conf'} "
        command += f"-s mj3h setup.json {steps} -l 0 >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_0.out",
            self.test_path / "swarm_0" / "gso_0.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_10.out",
            self.test_path / "swarm_0" / "gso_10.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_20.out",
            self.test_path / "swarm_0" / "gso_20.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_30.out",
            self.test_path / "swarm_0" / "gso_30.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            self.test_path / "init" / "swarm_centers.pdb",
        )
