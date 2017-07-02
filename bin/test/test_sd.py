"""Regression tests for testing SD scoring function"""

import shutil
import os
import filecmp

from bin.test.regression import RegressionTest


class TestRegressionSDShort(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_sd_short/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_sd_short/'
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_rec.pdb.H'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_lig.pdb.H'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_1ppe_10_step_10_glowworms_1_cluster(self):
        os.chdir(self.test_path)
        num_clusters = 1
        num_glowworms = 10
        steps = 10
        command = "lightdock %s %s %d %d %d -c 1 -nm --noxt -s " \
                  "sd > test_lightdock.out" % (self.test_path + '1PPE_rec.pdb.H',
                                               self.test_path + '1PPE_lig.pdb.H',
                                               num_clusters,
                                               num_glowworms,
                                               steps)
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_0.out',
                           self.test_path + 'cluster_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_10.out',
                           self.test_path + 'cluster_0/gso_10.out')


class TestRegressionSDLong(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_sd_long/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_sd_long/'
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_rec.pdb.H'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_lig.pdb.H'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_1ppe_30_steps_50_glowworms_1_cluster(self):
        if os.environ.has_key('LIGHTDOCK_LONG_TEST') and os.environ['LIGHTDOCK_LONG_TEST'] == 'true':
            os.chdir(self.test_path)
            num_clusters = 1
            num_glowworms = 50
            steps = 30
            command = "lightdock %s %s %d %d %d -c 1 -nm --noxt -s " \
                      "sd > test_lightdock.out" % (self.test_path + '1PPE_rec.pdb.H',
                                                   self.test_path + '1PPE_lig.pdb.H',
                                                   num_clusters,
                                                   num_glowworms,
                                                   steps)
            os.system(command)

            assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_0.out',
                               self.test_path + 'cluster_0/gso_0.out')
            # Uncomment to check, problems with number instability:
            # assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_10.out',
            #                   self.test_path + 'cluster_0/gso_10.out')
            # assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_20.out',
            #                   self.test_path + 'cluster_0/gso_20.out')
            # assert filecmp.cmp(self.golden_data_path + 'cluster_0/gso_30.out',
            #                   self.test_path + 'cluster_0/gso_30.out')
