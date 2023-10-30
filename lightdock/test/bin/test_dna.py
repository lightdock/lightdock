"""Regression tests for testing DNA scoring function"""

import shutil
import os
import filecmp
from pathlib import Path


class TestRegressionDNAShort:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_dna_short"

    def test_lightdock_1diz_10_steps_20_glowworms(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(self.golden_data_path / "1DIZ_rec.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "1DIZ_lig.pdb", tmp_path)

        num_glowworms = 25
        steps = 10

        command = (
            f"lgd_setup.py 1DIZ_rec.pdb 1DIZ_lig.pdb -g {num_glowworms} -anm"
        )
        command += ">> test_lightdock.out"
        os.system(command)

        command = (
            f"lgd_run.py -c 1 -s dna setup.json {steps} -l 0 >> test_lightdock.out"
        )
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_0.out",
            tmp_path / "swarm_0" / "gso_0.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_10.out",
            tmp_path / "swarm_0" / "gso_10.out",
        )


class TestRegressionDNARestraints:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_dna_restraints"

    def test_lightdock_1diz_30_steps_50_glowworms(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(self.golden_data_path / "1DIZ_rec.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "1DIZ_lig.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "restraints.list", tmp_path)

        num_glowworms = 50
        steps = 30

        command = f"lgd_setup.py 1DIZ_rec.pdb 1DIZ_lig.pdb -g {num_glowworms} -anm "
        command += "-rst restraints.list >> test_lightdock.out"
        os.system(command)

        command = (
            f"lgd_run.py -c 1 -s dna setup.json {steps} -l 0 >> test_lightdock.out"
        )
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
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_30.out",
            tmp_path / "swarm_0" / "gso_30.out",
        )
