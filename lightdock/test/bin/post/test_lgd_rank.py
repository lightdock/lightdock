"""Test for lgd_rank post script"""

import os
from pathlib import Path
import filecmp
import shutil


class TestGenerateRanking:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "4IZ7"

    def test_rank_with_clusters(self, tmp_path):
        num_swarms = 4
        num_steps = 10

        os.chdir(tmp_path)
        for i in range(num_swarms):
            swarm_dir = f"swarm_{i}"
            os.mkdir(swarm_dir)
            shutil.copyfile(
                self.golden_data_path / swarm_dir / f"gso_{num_steps}.out",
                tmp_path / swarm_dir / f"gso_{num_steps}.out",
            )
            shutil.copyfile(
                self.golden_data_path / swarm_dir / "cluster.repr",
                tmp_path / swarm_dir / "cluster.repr",
            )

        command = f"lgd_rank.py {num_swarms} {num_steps} > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "rank_by_scoring.list",
            tmp_path / "rank_by_scoring.list",
        )


class TestGenerateRankingWithoutClusters:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "4IZ7"

    def test_rank_with_clusters(self, tmp_path):
        num_swarms = 4
        num_steps = 10

        os.chdir(tmp_path)
        for i in range(num_swarms):
            swarm_dir = f"swarm_{i}"
            os.mkdir(swarm_dir)
            shutil.copyfile(
                self.golden_data_path / swarm_dir / f"gso_{num_steps}.out",
                tmp_path / swarm_dir / f"gso_{num_steps}.out",
            )
            shutil.copyfile(
                self.golden_data_path / swarm_dir / "cluster.repr",
                tmp_path / swarm_dir / "cluster.repr",
            )

        command = f"lgd_rank.py {num_swarms} {num_steps} --ignore_clusters > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "rank_by_scoring_noclust.list",
            tmp_path / "rank_by_scoring.list",
        )
