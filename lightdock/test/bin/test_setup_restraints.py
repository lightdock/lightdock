"""Regression tests for testing simulation setup with restraints information"""

import shutil
import os
import filecmp
from pathlib import Path
from lightdock.test.support import compare_two_files


class TestSetupWithRestraints:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_setup_rst"

    def test_lightdock_setup_automatic(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(self.golden_data_path / "2UUY_rec.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "2UUY_lig.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "restraints.list", tmp_path)

        num_glowworms = 100

        command = f"lgd_setup.py 2UUY_rec.pdb 2UUY_lig.pdb -g {num_glowworms} -anm --noxt --noh --now "
        command += "-rst restraints.list -sp >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            tmp_path / "init" / "swarm_centers.pdb",
        )
        assert compare_two_files(
            tmp_path / "setup.json", self.golden_data_path / "setup.json",
            ignore=["setup_version", "start_time"]
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_0.dat",
            tmp_path / "init" / "initial_positions_0.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_18.dat",
            tmp_path / "init" / "initial_positions_18.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_2UUY_rec.pdb",
            tmp_path / "lightdock_2UUY_rec.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_2UUY_lig.pdb",
            tmp_path / "lightdock_2UUY_lig.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "starting_poses_15.bild",
            tmp_path / "init" / "starting_poses_15.bild",
        )


class TestSetupWithRestraintsAndInsertCodes:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_setup_rst_insert"

    def test_lightdock_setup_automatic(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(self.golden_data_path / "receptor.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "ligand.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "restraints.list", tmp_path)

        command = "lgd_setup.py receptor.pdb ligand.pdb --noh --noxt --now -anm "
        command += "-rst restraints.list >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            tmp_path / "init" / "swarm_centers.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_0.dat",
            tmp_path / "init" / "initial_positions_0.dat",
        )
        assert compare_two_files(
            tmp_path / "setup.json", self.golden_data_path / "setup.json",
            ignore=["setup_version", "start_time"]
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_receptor.pdb",
            tmp_path / "lightdock_receptor.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_ligand.pdb",
            tmp_path / "lightdock_ligand.pdb",
        )


class TestSetupWithRestraintsAndFlipMode:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_setup_rst_flip"

    def test_lightdock_setup_automatic(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(self.golden_data_path / "3V6Z_A.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "3V6Z_B.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "restraints.list", tmp_path)

        command = "lgd_setup.py 3V6Z_A.pdb 3V6Z_B.pdb --noh --noxt -anm -flip "
        command += "-rst restraints.list >> test_lightdock.out"
        os.system(command)

        assert filecmp.cmp(
            self.golden_data_path / "init" / "swarm_centers.pdb",
            tmp_path / "init" / "swarm_centers.pdb",
        )
        assert compare_two_files(
            tmp_path / "setup.json", self.golden_data_path / "setup.json",
            ignore=["setup_version", "start_time"]
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_0.dat",
            tmp_path / "init" / "initial_positions_0.dat",
        )
