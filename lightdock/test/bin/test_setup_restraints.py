"""Regression tests for testing simulation setup with restraints information"""

import shutil
import os
import filecmp
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest
from lightdock.test.support import compare_two_files


class TestSetupWithRestraints(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_setup_restraints"
        self.golden_data_path = self.path / "golden_data" / "regression_setup_rst"

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "2UUY_rec.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "2UUY_lig.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "restraints.list", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_setup_automatic(self):
        os.chdir(self.test_path)

        num_glowworms = 100

        command = f"lightdock3_setup.py 2UUY_rec.pdb 2UUY_lig.pdb -g {num_glowworms} -anm --noxt --noh --now "
        command += "-rst restraints.list -sp >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            self.test_path / "init" / "swarm_centers.pdb",
        )
        assert compare_two_files(
            self.test_path / "setup.json", self.golden_data_path / "setup.json",
            ignore=["setup_version", "start_time"]
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_0.dat",
            self.test_path / "init" / "initial_positions_0.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_18.dat",
            self.test_path / "init" / "initial_positions_18.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_2UUY_rec.pdb",
            self.test_path / "lightdock_2UUY_rec.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_2UUY_lig.pdb",
            self.test_path / "lightdock_2UUY_lig.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "starting_poses_15.bild",
            self.test_path / "init" / "starting_poses_15.bild",
        )


class TestSetupWithRestraintsAndInsertCodes(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_setup_restraints_insert"
        self.golden_data_path = (
            self.path / "golden_data" / "regression_setup_rst_insert"
        )

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "receptor.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "ligand.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "restraints.list", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_setup_automatic(self):
        os.chdir(self.test_path)

        command = "lightdock3_setup.py receptor.pdb ligand.pdb --noh --noxt --now -anm "
        command += "-rst restraints.list >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            self.test_path / "init" / "swarm_centers.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_0.dat",
            self.test_path / "init" / "initial_positions_0.dat",
        )
        assert compare_two_files(
            self.test_path / "setup.json", self.golden_data_path / "setup.json",
            ignore=["setup_version", "start_time"]
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_receptor.pdb",
            self.test_path / "lightdock_receptor.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_ligand.pdb",
            self.test_path / "lightdock_ligand.pdb",
        )


class TestSetupWithRestraintsAndFlipMode(RegressionTest):
    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_regression_setup_rst_flip"
        self.golden_data_path = self.path / "golden_data" / "regression_setup_rst_flip"

    def setup(self):
        self.ini_path()
        shutil.copy(self.golden_data_path / "3V6Z_A.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "3V6Z_B.pdb", self.test_path)
        shutil.copy(self.golden_data_path / "restraints.list", self.test_path)

    def teardown(self):
        self.clean_path()

    def test_lightdock_setup_automatic(self):
        os.chdir(self.test_path)

        command = "lightdock3_setup.py 3V6Z_A.pdb 3V6Z_B.pdb --noh --noxt -anm -flip "
        command += "-rst restraints.list >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            self.test_path / "init" / "swarm_centers.pdb",
        )
        assert compare_two_files(
            self.test_path / "setup.json", self.golden_data_path / "setup.json",
            ignore=["setup_version", "start_time"]
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_0.dat",
            self.test_path / "init" / "initial_positions_0.dat",
        )
