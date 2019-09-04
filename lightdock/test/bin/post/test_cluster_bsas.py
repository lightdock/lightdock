"""Test for lgd_cluster_bsas post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestClusterBSAS(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_cluster_bsas/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) \
                                + '/golden_data/cluster_bsas/'

    def teardown(self):
        self.clean_test_path()

    def test_cluster_bsas(self):
        os.chdir(self.test_path)
        shutil.copyfile(self.golden_data_path + 'gso_10.out', self.test_path + 'gso_10.out')
        for i in range(10):
            shutil.copyfile(self.golden_data_path + 'lightdock_%d.pdb' % i, self.test_path + 'lightdock_%d.pdb' % i)

        command = "lgd_cluster_bsas.py %s > test.out" % (self.test_path + 'gso_10.out')
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'cluster.repr', self.test_path + 'cluster.repr')
