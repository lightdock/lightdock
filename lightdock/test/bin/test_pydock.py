"""Regression tests for testing CPyDock scoring function"""

import shutil
import os
import filecmp
from pathlib import Path


class TestRegressionPyDockRestraints:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_pydock_restraints"

    def test_lightdock_1ppe_20_steps_25_glowworms_rst(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(self.golden_data_path / "1PPE_rec.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "1PPE_lig.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "restraints.list", tmp_path)

        num_glowworms = 25
        steps = 20

        command = f"lgd_setup.py 1PPE_rec.pdb 1PPE_lig.pdb -g {num_glowworms} -anm "
        command += "-rst restraints.list >> test_lightdock.out"
        os.system(command)

        command = f"lgd_run.py -c 1 -s cpydock setup.json {steps} -l 0 >> test_lightdock.out"
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


class TestRegressionCPyDockLong:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_pydock_long"

    def test_lightdock_1ppe_40_steps_50_glowworms(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(self.golden_data_path / "1PPE_rec.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "1PPE_lig.pdb", tmp_path)

        num_glowworms = 50
        steps = 40

        command = f"lgd_setup.py 1PPE_rec.pdb 1PPE_lig.pdb -g {num_glowworms} >> test_lightdock.out"
        os.system(command)

        command = f"lgd_run.py -c 1 -s cpydock setup.json {steps} -l 100 >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_0.out",
            tmp_path / "swarm_100" / "gso_0.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_10.out",
            tmp_path / "swarm_100" / "gso_10.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_20.out",
            tmp_path / "swarm_100" / "gso_20.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_30.out",
            tmp_path / "swarm_100" / "gso_30.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_40.out",
            tmp_path / "swarm_100" / "gso_40.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            tmp_path / "init" / "swarm_centers.pdb",
        )
