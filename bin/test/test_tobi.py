"""Regression tests for TOBI scoring function"""

import shutil
import os
import filecmp
from .regression import RegressionTest


class TestRegressionTOBIShort(RegressionTest):
    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_tobi_short/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_tobi_short/'
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_lig.pdb'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_1ppe_4_steps_5_glowworms(self):
        os.chdir(self.test_path)
        num_swarms = 1
        num_glowworms = 5
        steps = 4

        command = "lightdock3_setup.py %s %s %d %d " \
                  "-ft %s > test_lightdock.out" % ('1PPE_rec.pdb',
                                                   '1PPE_lig.pdb',
                                                   num_swarms,
                                                   num_glowworms,
                                                   self.golden_data_path + '1PPE.ftdock')
        os.system(command)
        command = "lightdock3.py %s %d -c 1 -f %s -s " \
                  "tobi >> test_lightdock.out" % (self.test_path + 'setup.json',
                                                  steps,
                                                  self.golden_data_path + 'glowworm.conf')
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_0.dat',
                           self.test_path + 'init/initial_positions_0.dat')
        assert filecmp.cmp(self.golden_data_path + 'init/starting_positions_0.pdb',
                           self.test_path + 'init/starting_positions_0.pdb')
        assert filecmp.cmp(self.golden_data_path + 'init/cluster_centers.pdb',
                           self.test_path + 'init/cluster_centers.pdb')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out',
                           self.test_path + 'swarm_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_4.out',
                           self.test_path + 'swarm_0/gso_4.out')


class TestRegressionTOBIRestraints(RegressionTest):
    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_tobi_restraints/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_tobi_restraints/'
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_lig.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, 'restraints.list'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_1ppe_10_steps_20_glowworms(self):
        os.chdir(self.test_path)
        num_swarms = 4
        num_glowworms = 20
        steps = 10

        command = "lightdock3_setup.py %s %s %d %d --noh --noxt " \
                  "-rst %s > test_lightdock.out" % ('2UUY_rec.pdb',
                                                    '2UUY_lig.pdb',
                                                    num_swarms,
                                                    num_glowworms,
                                                    self.golden_data_path + 'restraints.list')
        os.system(command)
        command = "lightdock3.py %s %d -c 1 -s " \
                  "tobi >> test_lightdock.out" % (self.test_path + 'setup.json',
                                                  steps)
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out',
                           self.test_path + 'swarm_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_10.out',
                           self.test_path + 'swarm_0/gso_10.out')


class TestRegressionTOBIMoreGlowworms(RegressionTest):
    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_tobi_long/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) \
                                + '/golden_data/regression_tobi_long/'
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_lig.pdb'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_1pee_2_steps_20_glowworms(self):
        if 'LIGHTDOCK_LONG_TEST' in os.environ and os.environ['LIGHTDOCK_LONG_TEST'] == 'true':
            os.chdir(self.test_path)
            num_swarms = 1
            num_glowworms = 20
            steps = 2
            
            command = "lightdock3_setup.py %s %s %d %d " \
                      "-ft %s > test_lightdock.out" % ('1PPE_rec.pdb',
                                                       '1PPE_lig.pdb',
                                                       num_swarms,
                                                       num_glowworms,
                                                       self.golden_data_path + '1PPE.ftdock')
            os.system(command)
            command = "lightdock3.py %s %d -c 1 -f %s -s " \
                      "tobi >> test_lightdock.out" % (self.test_path + 'setup.json',
                                                      steps,
                                                      self.golden_data_path + 'glowworm.conf')
            os.system(command)

            assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_0.dat',
                               self.test_path + 'init/initial_positions_0.dat')
            assert filecmp.cmp(self.golden_data_path + 'init/starting_positions_0.pdb',
                               self.test_path + 'init/starting_positions_0.pdb')
            assert filecmp.cmp(self.golden_data_path + 'init/cluster_centers.pdb',
                               self.test_path + 'init/cluster_centers.pdb')
            assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out',
                               self.test_path + 'swarm_0/gso_0.out')
            assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_2.out',
                               self.test_path + 'swarm_0/gso_2.out')
