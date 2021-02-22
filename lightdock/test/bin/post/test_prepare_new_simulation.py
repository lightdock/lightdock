"""Test for lgd_prepare_new_simulation post script"""

import os
import filecmp
import shutil
from pathlib import Path
from lightdock.test.bin.regression import RegressionTest


class TestPrepareNewSimulation(RegressionTest):

    def __init__(self):
        super().__init__()
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / 'scratch_prepare_new_simulation'
        self.golden_data_path = self.path / 'golden_data' / 'prepare_new_simulation'

    def setup(self):
        self.ini_test_path()

    def teardown(self):
        self.clean_test_path()

    def test_prepare_new_simulation(self):
        # Prepare folder structure for this test
        os.chdir(self.test_path)
        swarm_dir = 'swarm_0'
        os.mkdir(swarm_dir)
        shutil.copyfile(self.golden_data_path / swarm_dir / 'gso_10.out',
                        self.test_path / swarm_dir / 'gso_10.out')
        shutil.copyfile(self.golden_data_path / swarm_dir / 'cluster.repr',
                        self.test_path / swarm_dir / 'cluster.repr')

        command = "lgd_prepare_new_simulation.py 0 10 new -nm > test.out"
        os.system(command)

        assert filecmp.cmp(self.golden_data_path / 'new' / 'init' / 'swarm_centers.pdb',
                           self.test_path / 'new' / 'init' / 'swarm_centers.pdb')
