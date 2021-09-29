"""Regression tests for testing CPyDock scoring function"""

import shutil
import os
import filecmp
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestRegressionPyDockRestraints(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_pydock_restraints"
        self.golden_data_path = (
            self.path / "golden_data" / "regression_pydock_restraints"
        )

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "1PPE_rec.pdb.H", self.test_path)
        shutil.copy(self.golden_data_path / "1PPE_lig.pdb.H", self.test_path)
        shutil.copy(self.golden_data_path / "restraints.list", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_2uuy_20_steps_25_glowworms_rst(self):
        os.chdir(self.test_path)
        num_glowworms = 25
        steps = 20

        command = f"lightdock3_setup.py 1PPE_rec.pdb.H 1PPE_lig.pdb.H -g {num_glowworms} -anm "
        command += "-rst restraints.list >> test_lightdock.out"
        os.system(command)

        command = f"lightdock3.py -c 1 -s cpydock setup.json {steps} -l 0 >> test_lightdock.out"
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


class TestRegressionCPyDockLong(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_pydock_long"
        self.golden_data_path = self.path / "golden_data" / "regression_pydock_long"

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "1PPE_rec.pdb.H", self.test_path)
        shutil.copy(self.golden_data_path / "1PPE_lig.pdb.H", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_2uuy_40_steps_50_glowworms(self):
        os.chdir(self.test_path)
        num_glowworms = 50
        steps = 40

        command = f"lightdock3_setup.py 1PPE_rec.pdb.H 1PPE_lig.pdb.H -g {num_glowworms} >> test_lightdock.out"
        os.system(command)

        command = f"lightdock3.py -c 1 -s cpydock setup.json {steps} -l 100 >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_0.out",
            self.test_path / "swarm_100" / "gso_0.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_10.out",
            self.test_path / "swarm_100" / "gso_10.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_20.out",
            self.test_path / "swarm_100" / "gso_20.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_30.out",
            self.test_path / "swarm_100" / "gso_30.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_100" / "gso_40.out",
            self.test_path / "swarm_100" / "gso_40.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            self.test_path / "init" / "swarm_centers.pdb",
        )
