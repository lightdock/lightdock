"""Test for lgd_top post script"""

import os
import filecmp
import shutil
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestFilterRestraints(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_filter_restraints"
        self.golden_data_path = self.path / "golden_data" / "filter_restraints"

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_filter_restraints(self):

        # Prepare folder structure for this test
        os.chdir(self.test_path)
        shutil.copyfile(
            self.golden_data_path / "rank_by_scoring.list",
            self.test_path / "rank_by_scoring.list",
        )
        shutil.copyfile(
            self.golden_data_path / "restraints.list",
            self.test_path / "restraints.list",
        )
        shutil.copytree(self.golden_data_path / "swarm_0", self.test_path / "swarm_0")
        shutil.copytree(self.golden_data_path / "swarm_1", self.test_path / "swarm_1")
        shutil.copytree(self.golden_data_path / "swarm_2", self.test_path / "swarm_2")
        shutil.copytree(self.golden_data_path / "swarm_3", self.test_path / "swarm_3")

        command = "lgd_filter_restraints.py -cutoff 5.0 -fnat 0.4 rank_by_scoring.list restraints.list A B > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "filtered" / "rank_filtered.list",
            self.test_path / "filtered" / "rank_filtered.list",
        )
        assert filecmp.cmp(
            self.golden_data_path / "filtered" / "swarm_1_34.pdb",
            self.test_path / "filtered" / "swarm_1_34.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "filtered" / "swarm_2_16.pdb",
            self.test_path / "filtered" / "swarm_2_16.pdb",
        )
