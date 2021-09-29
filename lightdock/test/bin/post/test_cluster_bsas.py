"""Test for lgd_cluster_bsas post script"""

import os
import filecmp
import shutil
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestClusterBSAS(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_cluster_bsas"
        self.golden_data_path = self.path / "golden_data" / "cluster_bsas"

    def setup(self):
        self.ini_path()

    def teardown(self):
        self.clean_path()

    def test_cluster_bsas(self):
        os.chdir(self.test_path)
        shutil.copyfile(
            self.golden_data_path / "gso_10.out", self.test_path / "gso_10.out"
        )
        for i in range(10):
            shutil.copyfile(
                self.golden_data_path / f"lightdock_{i}.pdb",
                self.test_path / f"lightdock_{i}.pdb",
            )

        command = f"lgd_cluster_bsas.py {self.test_path / 'gso_10.out'} > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "cluster.repr", self.test_path / "cluster.repr"
        )
