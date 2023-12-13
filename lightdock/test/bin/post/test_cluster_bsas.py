"""Test for lgd_cluster_bsas post script"""

import os
import filecmp
import shutil
from pathlib import Path


class TestClusterBSAS:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "cluster_bsas"

    def test_cluster_bsas_default(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copyfile(
            self.golden_data_path / "gso_10.out", tmp_path / "gso_10.out"
        )
        for i in range(10):
            shutil.copyfile(
                self.golden_data_path / f"lightdock_{i}.pdb",
                tmp_path / f"lightdock_{i}.pdb",
            )

        command = f"lgd_cluster_bsas.py {tmp_path / 'gso_10.out'} > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "cluster.repr", tmp_path / "cluster.repr"
        )

    def test_cluster_bsas_with_4_cutoff(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copyfile(
            self.golden_data_path / "gso_10.out", tmp_path / "gso_10.out"
        )
        for i in range(10):
            shutil.copyfile(
                self.golden_data_path / f"lightdock_{i}.pdb",
                tmp_path / f"lightdock_{i}.pdb",
            )

        command = f"lgd_cluster_bsas.py {tmp_path / 'gso_10.out'} -c 4.0 > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "cluster.repr", tmp_path / "cluster.repr"
        )

    def test_cluster_bsas_with_10_cutoff(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copyfile(
            self.golden_data_path / "gso_10.out", tmp_path / "gso_10.out"
        )
        for i in range(10):
            shutil.copyfile(
                self.golden_data_path / f"lightdock_{i}.pdb",
                tmp_path / f"lightdock_{i}.pdb",
            )

        command = f"lgd_cluster_bsas.py {tmp_path / 'gso_10.out'} -c 10.0 > test.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "cluster.repr.10", tmp_path / "cluster.repr"
        )