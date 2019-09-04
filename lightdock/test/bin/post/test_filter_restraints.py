"""Test for lgd_top post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestFilterRestraints(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_filter_restraints/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/filter_restraints/'

    def teardown(self):
        self.clean_test_path()

    def test_filter_restraints(self):

        # Prepare folder structure for this test
        os.chdir(self.test_path)
        shutil.copyfile(os.path.join(self.golden_data_path, 'rank_by_scoring.list'),
                        os.path.join(self.test_path, 'rank_by_scoring.list'))
        shutil.copyfile(os.path.join(self.golden_data_path, 'restraints.list'),
                        os.path.join(self.test_path, 'restraints.list'))
        shutil.copytree(os.path.join(self.golden_data_path, 'swarm_0'),
                        os.path.join(self.test_path, 'swarm_0'))
        shutil.copytree(os.path.join(self.golden_data_path, 'swarm_1'),
                        os.path.join(self.test_path, 'swarm_1'))
        shutil.copytree(os.path.join(self.golden_data_path, 'swarm_2'),
                        os.path.join(self.test_path, 'swarm_2'))
        shutil.copytree(os.path.join(self.golden_data_path, 'swarm_3'),
                        os.path.join(self.test_path, 'swarm_3'))

        command = "lgd_filter_restraints.py -cutoff 5.0 -fnat 0.4 rank_by_scoring.list restraints.list A B > test.out"
        os.system(command)

        assert filecmp.cmp(os.path.join(self.golden_data_path, 'filtered', 'rank_filtered.list'), 
                           os.path.join(self.test_path, 'filtered', 'rank_filtered.list'))
        assert filecmp.cmp(os.path.join(self.golden_data_path, 'filtered', 'swarm_1_34.pdb'), 
                           os.path.join(self.test_path, 'filtered', 'swarm_1_34.pdb'))
        assert filecmp.cmp(os.path.join(self.golden_data_path, 'filtered', 'swarm_2_16.pdb'), 
                           os.path.join(self.test_path, 'filtered', 'swarm_2_16.pdb'))
 