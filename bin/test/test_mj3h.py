"""Regression tests for testing MJ3H scoring function"""

import shutil
import os
import filecmp
from bin.test.regression import RegressionTest


class TestRegressionMJ3HShort(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_mj3h_short/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_mj3h_short/'
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_lig.pdb'), self.test_path)

    def teardown(self):
        self.clean_test_path()
    
    def test_lightdock_2uuy_10_steps_50_glowworms(self):
        os.chdir(self.test_path)
        num_clusters = 1
        num_glowworms = 25
        steps = 10
        command = "lightdock %s %s %d %d %d -f %s" \
                  " -s mj3h -c 1 > lightdock.out" % (self.test_path + '2UUY_rec.pdb',
                                                     self.test_path + '2UUY_lig.pdb',
                                                     num_clusters,
                                                     num_glowworms,
                                                     steps,
                                                     self.golden_data_path + 'glowworm.conf')
        os.system(command)
        
        assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_0.dat',
                           self.test_path + 'init/initial_positions_0.dat')
        assert filecmp.cmp(self.golden_data_path + 'init/starting_positions_0.pdb',
                           self.test_path + 'init/starting_positions_0.pdb')
        assert filecmp.cmp(self.golden_data_path + 'init/cluster_centers.pdb',
                           self.test_path + 'init/cluster_centers.pdb')
        assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_0.out',
                           self.test_path + 'cluster_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_10.out',
                           self.test_path + 'cluster_0/gso_10.out')


class TestRegressionMJ3HLong(RegressionTest):
    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_mj3h_long/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_mj3h_long/'
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_lig.pdb'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_2uuy_10_steps_50_glowworms(self):
        if os.environ.has_key('LIGHTDOCK_LONG_TEST') and os.environ['LIGHTDOCK_LONG_TEST'] == 'true':
            os.chdir(self.test_path)
            num_clusters = 1
            num_glowworms = 50
            steps = 20
            command = "lightdock %s %s %d %d %d -f %s" \
                      " -s mj3h -c 1 > lightdock.out" % (self.test_path + '2UUY_rec.pdb',
                                                         self.test_path + '2UUY_lig.pdb',
                                                         num_clusters,
                                                         num_glowworms,
                                                         steps,
                                                         self.golden_data_path + 'glowworm.conf')
            os.system(command)

            assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_0.dat',
                               self.test_path + 'init/initial_positions_0.dat')
            assert filecmp.cmp(self.golden_data_path + 'init/starting_positions_0.pdb',
                               self.test_path + 'init/starting_positions_0.pdb')
            assert filecmp.cmp(self.golden_data_path + 'init/cluster_centers.pdb',
                               self.test_path + 'init/cluster_centers.pdb')
            assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_0.out',
                               self.test_path + 'cluster_0/gso_0.out')
            assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_10.out',
                               self.test_path + 'cluster_0/gso_10.out')
            assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_20.out',
                               self.test_path + 'cluster_0/gso_20.out')