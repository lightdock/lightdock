"""Test for lgd_rank post script"""

import os
from pathlib import Path
import filecmp
import shutil
from lightdock.test.bin.regression import RegressionTest


class TestGenerateRanking(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_lgd_rank"
        self.golden_data_path = self.path / "golden_data" / "4IZ7"

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_rank_with_clusters(self):
        num_swarms = 4
        num_steps = 10

        # Prepare folder structure for this test
        os.chdir(self.test_path)
        for i in range(num_swarms):
            swarm_dir = f"swarm_{i}"
            os.mkdir(swarm_dir)
            shutil.copyfile(
                self.golden_data_path / swarm_dir / f"gso_{num_steps}.out",
                self.test_path / swarm_dir / f"gso_{num_steps}.out",
            )
            shutil.copyfile(
                self.golden_data_path / swarm_dir / "cluster.repr",
                self.test_path / swarm_dir / "cluster.repr",
            )

        command = f"lgd_rank.py {num_swarms} {num_steps} > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "rank_by_scoring.list",
            self.test_path / "rank_by_scoring.list",
        )


class TestGenerateRankingWithoutClusters(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_lgd_rank_no_clust"
        self.golden_data_path = self.path / "golden_data" / "4IZ7"

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_rank_with_clusters(self):
        num_swarms = 4
        num_steps = 10

        # Prepare folder structure for this test
        os.chdir(self.test_path)
        for i in range(num_swarms):
            swarm_dir = f"swarm_{i}"
            os.mkdir(swarm_dir)
            shutil.copyfile(
                self.golden_data_path / swarm_dir / f"gso_{num_steps}.out",
                self.test_path / swarm_dir / f"gso_{num_steps}.out",
            )
            shutil.copyfile(
                self.golden_data_path / swarm_dir / "cluster.repr",
                self.test_path / swarm_dir / "cluster.repr",
            )

        command = f"lgd_rank.py {num_swarms} {num_steps} --ignore_clusters > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "rank_by_scoring_noclust.list",
            self.test_path / "rank_by_scoring.list",
        )
