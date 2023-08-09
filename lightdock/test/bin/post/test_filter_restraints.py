"""Test for lgd_top post script"""

import os
import filecmp
import shutil
from pathlib import Path


class TestFilterRestraints:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "filter_restraints"

    def test_filter_restraints(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copyfile(
            self.golden_data_path / "rank_by_scoring.list",
            tmp_path / "rank_by_scoring.list",
        )
        shutil.copyfile(
            self.golden_data_path / "restraints.list",
            tmp_path / "restraints.list",
        )
        shutil.copytree(self.golden_data_path / "swarm_0", tmp_path / "swarm_0")
        shutil.copytree(self.golden_data_path / "swarm_1", tmp_path / "swarm_1")
        shutil.copytree(self.golden_data_path / "swarm_2", tmp_path / "swarm_2")
        shutil.copytree(self.golden_data_path / "swarm_3", tmp_path / "swarm_3")

        command = "lgd_filter_restraints.py -cutoff 5.0 -fnat 0.4 rank_by_scoring.list restraints.list A B > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "filtered" / "rank_filtered.list",
            tmp_path / "filtered" / "rank_filtered.list",
        )
        assert filecmp.cmp(
            self.golden_data_path / "filtered" / "swarm_1_34.pdb",
            tmp_path / "filtered" / "swarm_1_34.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "filtered" / "swarm_2_16.pdb",
            tmp_path / "filtered" / "swarm_2_16.pdb",
        )
