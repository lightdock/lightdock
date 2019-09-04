"""Regression tests for testing DNA scoring function"""

import shutil
import os
import filecmp

from .regression import RegressionTest


class TestRegressionDNAShort(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_dna_short/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_dna_short/'
        shutil.copy(os.path.join(self.golden_data_path, '1DIZ_rec.pdb.H'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1DIZ_lig.pdb.H'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_1diz_10_steps_20_glowworms_1_swarm(self):
        os.chdir(self.test_path)
        num_swarms = 1
        num_glowworms = 20
        steps = 10

        command = "lightdock3_setup.py %s %s %d %d -anm > test_lightdock.out" % ('1DIZ_rec.pdb.H',
                                                                             '1DIZ_lig.pdb.H',
                                                                             num_swarms,
                                                                             num_glowworms
                                                                             )
        os.system(command)
        command = "lightdock3.py %s %d -c 1 -s " \
                  "dna >> test_lightdock.out" % (self.test_path + 'setup.json', steps)
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out',
                           self.test_path + 'swarm_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_10.out',
                           self.test_path + 'swarm_0/gso_10.out')


class TestRegressionDNARestraintsOnlyReceptor(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_dna_restraints/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_dna_restraints/'
        shutil.copy(os.path.join(self.golden_data_path, '1DIZ_rec.pdb.H'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1DIZ_lig.pdb.H'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_1diz_5_steps_10_glowworms_4_swarms(self):
        os.chdir(self.test_path)
        num_swarms = 4
        num_glowworms = 10
        steps = 5

        command = "lightdock3_setup.py %s %s %d %d -rst %s > test_lightdock.out" % ('1DIZ_rec.pdb.H',
                                                                                '1DIZ_lig.pdb.H',
                                                                                num_swarms,
                                                                                num_glowworms,
                                                                                self.golden_data_path + 'restraints.list'
                                                                                )
        os.system(command)
        command = "lightdock3.py %s %d -c 1 -s " \
                  "dna >> test_lightdock.out" % (self.test_path + 'setup.json', steps)
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out',
                           self.test_path + 'swarm_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_5.out',
                           self.test_path + 'swarm_0/gso_5.out')
