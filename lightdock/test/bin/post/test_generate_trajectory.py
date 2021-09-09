"""Test for lgd_generate_trajectory post script"""

import os
import filecmp
import shutil
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestGenerateTrajectory(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_generate_trajectory"
        self.golden_data_path = self.path / "golden_data" / "generate_trajectory"

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_generate_trajectory(self):
        # Prepare folder structure for this test
        os.chdir(self.test_path)
        swarm_dir = "swarm_0"
        os.mkdir(swarm_dir)
        shutil.copyfile(
            self.golden_data_path / swarm_dir / "gso_0.out",
            self.test_path / swarm_dir / "gso_0.out",
        )
        shutil.copyfile(
            self.golden_data_path / swarm_dir / "gso_10.out",
            self.test_path / swarm_dir / "gso_10.out",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_4IZ7_A_noh.pdb",
            self.test_path / "lightdock_4IZ7_A_noh.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_4IZ7_B_noh.pdb",
            self.test_path / "lightdock_4IZ7_B_noh.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "setup.json", self.test_path / "setup.json"
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_rec.nm.npy",
            self.test_path / "lightdock_rec.nm.npy",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_lig.nm.npy",
            self.test_path / "lightdock_lig.nm.npy",
        )

        os.chdir(self.test_path / swarm_dir)
        receptor = self.test_path / "lightdock_4IZ7_A_noh.pdb"
        ligand = self.test_path / "lightdock_4IZ7_B_noh.pdb"
        setup = self.test_path / "setup.json"
        command = (
            f"lgd_generate_trajectory.py 4 10 {receptor} {ligand} {setup} > test.out"
        )
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "trajectory_4_step_0.pdb",
            self.test_path / "swarm_0" / "trajectory_4_step_0.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "trajectory_4_step_10.pdb",
            self.test_path / "swarm_0" / "trajectory_4_step_10.pdb",
        )
