"""Regression tests for testing Mj3h scoring function"""

import shutil
import os
import filecmp
from pathlib import Path


class TestRegressionMj3hShort:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_mj3h_short"

    def test_lightdock_2uuy_10_steps_25_glowworms_100_swarms(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(self.golden_data_path / "2UUY_rec.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "2UUY_lig.pdb", tmp_path)

        num_swarms = 100
        num_glowworms = 25
        steps = 10

        command = f"lgd_setup.py 2UUY_rec.pdb 2UUY_lig.pdb -g {num_glowworms} -s {num_swarms} "
        command += ">> test_lightdock.out"
        os.system(command)

        command = f"lgd_run.py -c 1 -f {self.golden_data_path / 'glowworm.conf'} "
        command += f"-s mj3h setup.json {steps} -l 0 >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_0.out",
            tmp_path / "swarm_0" / "gso_0.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_10.out",
            tmp_path / "swarm_0" / "gso_10.out",
        )


class TestRegressionMj3hLong:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_mj3h_long"

    def test_lightdock_1ppe_50_steps_20_glowworms_100_swarms(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(self.golden_data_path / "1PPE_rec.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "1PPE_lig.pdb", tmp_path)

        num_swarms = 100
        num_glowworms = 20
        steps = 30

        command = f"lgd_setup.py 1PPE_rec.pdb 1PPE_lig.pdb -g {num_glowworms} -s {num_swarms} "
        command += ">> test_lightdock.out"
        os.system(command)

        command = f"lgd_run.py -c 1 -f {self.golden_data_path / 'glowworm.conf'} "
        command += f"-s mj3h setup.json {steps} -l 0 >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_0.out",
            tmp_path / "swarm_0" / "gso_0.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_10.out",
            tmp_path / "swarm_0" / "gso_10.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_20.out",
            tmp_path / "swarm_0" / "gso_20.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_30.out",
            tmp_path / "swarm_0" / "gso_30.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            tmp_path / "init" / "swarm_centers.pdb",
        )
