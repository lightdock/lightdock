"""Regression tests for testing simulation setup using membrane"""

import shutil
import os
import filecmp
from pathlib import Path
from lightdock.test.support import compare_two_files


class TestSetupWithMembrane:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_setup_membrane_all_anm"

    def test_lightdock_setup_with_membrane_automatic(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(
            self.golden_data_path / "3x29_receptor_membrane.pdb", tmp_path
        )
        shutil.copy(self.golden_data_path / "3x29_ligand.pdb", tmp_path)

        command = "lgd_setup.py 3x29_receptor_membrane.pdb 3x29_ligand.pdb "
        command += "--noxt --noh -membrane -anm >> test_lightdock.out"
        os.system(command)

        assert compare_two_files(
            tmp_path / "setup.json", self.golden_data_path / "setup.json",
            ignore=["setup_version", "start_time"]
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_3x29_receptor_membrane_mask.npy",
            tmp_path / "lightdock_3x29_receptor_membrane_mask.npy",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_3x29_ligand_mask.npy",
            tmp_path / "lightdock_3x29_ligand_mask.npy",
        )


class TestSetupWithMembraneANM:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_setup_membrane"

    def test_lightdock_setup_with_membrane_manual(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(
            self.golden_data_path / "3x29_receptor_membrane.pdb", tmp_path
        )
        shutil.copy(self.golden_data_path / "3x29_ligand.pdb", tmp_path)

        num_swarms = 400
        num_glowworms = 50

        command = f"lgd_setup.py 3x29_receptor_membrane.pdb 3x29_ligand.pdb -g {num_glowworms} "
        command += f" -s {num_swarms} --noxt --noh -membrane >> test_lightdock.out"
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
            self.golden_data_path / "init" / "initial_positions_35.dat",
            tmp_path / "init" / "initial_positions_35.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "init" / "initial_positions_111.dat",
            tmp_path / "init" / "initial_positions_111.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_3x29_receptor_membrane.pdb",
            tmp_path / "lightdock_3x29_receptor_membrane.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_3x29_ligand.pdb",
            tmp_path / "lightdock_3x29_ligand.pdb",
        )


class TestSetupWithMembraneAndRestraints:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "regression_setup_membrane_restraints"

    def test_lightdock_setup_with_membrane_automatic(self, tmp_path):
        os.chdir(tmp_path)
        shutil.copy(
            self.golden_data_path / "3x29_receptor_membrane.pdb", tmp_path
        )
        shutil.copy(self.golden_data_path / "3x29_ligand.pdb", tmp_path)
        shutil.copy(self.golden_data_path / "restraints.list", tmp_path)

        num_swarms = 400
        num_glowworms = 50

        command = f"lgd_setup.py 3x29_receptor_membrane.pdb 3x29_ligand.pdb -g {num_glowworms} "
        command += f" -s {num_swarms} --noxt --noh -membrane -rst restraints.list >> test_lightdock.out"
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
            self.golden_data_path / "init" / "initial_positions_5.dat",
            tmp_path / "init" / "initial_positions_5.dat",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_3x29_receptor_membrane.pdb",
            tmp_path / "lightdock_3x29_receptor_membrane.pdb",
        )
        assert filecmp.cmp(
            self.golden_data_path / "lightdock_3x29_ligand.pdb",
            tmp_path / "lightdock_3x29_ligand.pdb",
        )
