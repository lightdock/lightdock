"""Regression tests for testing CPyDock scoring function"""

import shutil
import os
import filecmp

from bin.test.regression import RegressionTest


class TestRegressionPyDockShort(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_pydock_short/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_pydock_short/'
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_rec.pdb.H'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_lig.pdb.H'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_rec.pdb.amber'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_lig.pdb.amber'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_1ppe_1_step_5_glowworms_1_cluster(self):
        os.chdir(self.test_path)
        num_clusters = 1
        num_glowworms = 5
        steps = 1
        command = "lightdock %s %s %d %d %d -c 1 -s " \
                  "cpydock > test_lightdock.out" % (self.test_path + '1PPE_rec.pdb.H',
                                                    self.test_path + '1PPE_lig.pdb.H',
                                                    num_clusters,
                                                    num_glowworms,
                                                    steps)
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_0.out',
                           self.test_path + 'cluster_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_1.out',
                           self.test_path + 'cluster_0/gso_1.out')


class TestRegressionPyDockLong(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_pydock_long/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_pydock_long/'
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_rec.pdb.H'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_lig.pdb.H'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_rec.pdb.amber'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_lig.pdb.amber'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_1ppe_10_steps_10_glowworms_1_cluster(self):
        if os.environ.has_key('LIGHTDOCK_LONG_TEST') and os.environ['LIGHTDOCK_LONG_TEST'] == 'true':
            os.chdir(self.test_path)
            num_clusters = 1
            num_glowworms = 10
            steps = 10
            command = "lightdock %s %s %d %d %d -c 1 -s " \
                      "cpydock > test_lightdock.out" % (self.test_path + '1PPE_rec.pdb.H',
                                                        self.test_path + '1PPE_lig.pdb.H',
                                                        num_clusters,
                                                        num_glowworms,
                                                        steps)
            os.system(command)

            assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_0.out',
                               self.test_path + 'cluster_0/gso_0.out')
            assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_10.out',
                               self.test_path + 'cluster_0/gso_10.out')
