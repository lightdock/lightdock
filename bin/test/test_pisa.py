"""Regression tests for PISA scoring function"""

import shutil
import os
import filecmp
from .regression import RegressionTest


class TestRegressionPISAShort(RegressionTest):
    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_pisa_short/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) \
                                + '/golden_data/regression_pisa_short/'
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_lig.pdb'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_2uuy_5_steps_5_glowworms(self):
        os.chdir(self.test_path)
        num_swarms = 1
        num_glowworms = 5
        steps = 5
        
        command = "lightdock3_setup.py %s %s %d %d > test_lightdock.out" % ('2UUY_rec.pdb',
                                                                        '2UUY_lig.pdb',
                                                                        num_swarms,
                                                                        num_glowworms
                                                                        )
        os.system(command)
        command = "lightdock3.py %s %d -c 1 -f %s -s " \
                  "pisa >> test_lightdock.out" % (self.test_path + 'setup.json',
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
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_5.out',
                           self.test_path + 'swarm_0/gso_5.out')


class TestRegressionPISALong(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_pisa_long/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) \
                                + '/golden_data/regression_pisa_long/'
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_lig.pdb'), self.test_path)

    def teardown(self):
        self.clean_test_path()
    
    def test_lightdock_2uuy_10_steps_50_glowworms(self):
        if 'LIGHTDOCK_LONG_TEST' in os.environ and os.environ['LIGHTDOCK_LONG_TEST'] == 'true':
            os.chdir(self.test_path)
            num_swarms = 1
            num_glowworms = 50
            steps = 10
            
            command = "lightdock3_setup.py %s %s %d %d > test_lightdock.out" % ('2UUY_rec.pdb',
                                                                            '2UUY_lig.pdb',
                                                                            num_swarms,
                                                                            num_glowworms
                                                                            )
            os.system(command)
            command = "lightdock3.py %s %d -c 1 -f %s -s " \
                      "pisa >> test_lightdock.out" % (self.test_path + 'setup.json',
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
            assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_10.out',
                               self.test_path + 'swarm_0/gso_10.out')


class TestRegressionPISARestraints(RegressionTest):
    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_pisa_restraints/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) \
                                + '/golden_data/regression_pisa_restraints/'
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_lig.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, 'restraints.list'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_2uuy_10_steps_5_glowworms(self):
        os.chdir(self.test_path)
        num_swarms = 4
        num_glowworms = 5
        steps = 10
        
        command = "lightdock3_setup.py %s %s %d %d -rst %s --noh > test_lightdock.out" % ('2UUY_rec.pdb',
                                                                                      '2UUY_lig.pdb',
                                                                                      num_swarms,
                                                                                      num_glowworms,
                                                                                      'restraints.list'
                                                                                      )
        os.system(command)
        command = "lightdock3.py %s %d -c 1 -s " \
                  "pisa >> test_lightdock.out" % (self.test_path + 'setup.json', 
                                                  steps)
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_0.dat',
                           self.test_path + 'init/initial_positions_0.dat')
        assert filecmp.cmp(self.golden_data_path + 'init/starting_positions_0.pdb',
                           self.test_path + 'init/starting_positions_0.pdb')
        assert filecmp.cmp(self.golden_data_path + 'init/cluster_centers.pdb',
                           self.test_path + 'init/cluster_centers.pdb')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out',
                           self.test_path + 'swarm_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_10.out',
                           self.test_path + 'swarm_0/gso_10.out')
