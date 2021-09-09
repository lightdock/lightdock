"""Regression tests for testing DNA scoring function"""

import shutil
import os
import filecmp
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestRegressionDNAShort(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_dna_short"
        self.golden_data_path = self.path / "golden_data" / "regression_dna_short"

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "1DIZ_rec.pdb.H", self.test_path)
        shutil.copy(self.golden_data_path / "1DIZ_lig.pdb.H", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_1diz_10_steps_20_glowworms(self):
        os.chdir(self.test_path)
        num_glowworms = 25
        steps = 10

        command = (
            f"lightdock3_setup.py 1DIZ_rec.pdb.H 1DIZ_lig.pdb.H -g {num_glowworms} -anm"
        )
        command += ">> test_lightdock.out"
        os.system(command)

        command = (
            f"lightdock3.py -c 1 -s dna setup.json {steps} -l 0 >> test_lightdock.out"
        )
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_0.out",
            self.test_path / "swarm_0" / "gso_0.out",
        )
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_10.out",
            self.test_path / "swarm_0" / "gso_10.out",
        )


class TestRegressionDNARestraints(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_dna_restraints"
        self.golden_data_path = self.path / "golden_data" / "regression_dna_restraints"

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "1DIZ_rec.pdb.H", self.test_path)
        shutil.copy(self.golden_data_path / "1DIZ_lig.pdb.H", self.test_path)
        shutil.copy(self.golden_data_path / "restraints.list", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_1diz_30_steps_50_glowworms(self):
        os.chdir(self.test_path)
        num_glowworms = 50
        steps = 30

        command = f"lightdock3_setup.py 1DIZ_rec.pdb.H 1DIZ_lig.pdb.H -g {num_glowworms} -anm "
        command += "-rst restraints.list >> test_lightdock.out"
        os.system(command)

        command = (
            f"lightdock3.py -c 1 -s dna setup.json {steps} -l 0 >> test_lightdock.out"
        )
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
        assert filecmp.cmp(
            self.golden_data_path / "swarm_0" / "gso_30.out",
            self.test_path / "swarm_0" / "gso_30.out",
        )
