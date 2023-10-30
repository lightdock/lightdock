"""Test for lgd_generate_trajectory post script"""

import os
import filecmp
import shutil
from pathlib import Path


class TestGenerateTrajectory:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "generate_trajectory"

    def test_generate_trajectory(self, tmp_path):
        os.chdir(tmp_path)
        swarm_dir = "swarm_0"
        os.mkdir(swarm_dir)
        shutil.copyfile(
            self.golden_data_path / swarm_dir / "gso_0.out",
            tmp_path / swarm_dir / "gso_0.out",
        )
        shutil.copyfile(
            self.golden_data_path / swarm_dir / "gso_10.out",
            tmp_path / swarm_dir / "gso_10.out",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_4IZ7_A_noh.pdb",
            tmp_path / "lightdock_4IZ7_A_noh.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_4IZ7_B_noh.pdb",
            tmp_path / "lightdock_4IZ7_B_noh.pdb",
        )
        shutil.copyfile(
            self.golden_data_path / "setup.json", tmp_path / "setup.json"
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_rec.nm.npy",
            tmp_path / "lightdock_rec.nm.npy",
        )
        shutil.copyfile(
            self.golden_data_path / "lightdock_lig.nm.npy",
            tmp_path / "lightdock_lig.nm.npy",
        )

        os.chdir(tmp_path / swarm_dir)
        receptor = tmp_path / "lightdock_4IZ7_A_noh.pdb"
        ligand = tmp_path / "lightdock_4IZ7_B_noh.pdb"
        setup = tmp_path / "setup.json"
        command = (
            f"lgd_generate_trajectory.py 4 10 {receptor} {ligand} {setup} > test.out"
        )
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "trajectory_4_step_0.pdb",
            tmp_path / "swarm_0" / "trajectory_4_step_0.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "trajectory_4_step_10.pdb",
            tmp_path / "swarm_0" / "trajectory_4_step_10.pdb",
        )
