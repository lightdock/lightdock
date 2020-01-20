"""Regression tests for testing CPyDock scoring function"""

import shutil
import os
import filecmp

from .regression import RegressionTest


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

    def test_lightdock_1ppe_1_step_5_glowworms_1_swarm(self):
        os.chdir(self.test_path)
        num_swarms = 1
        num_glowworms = 5
        steps = 1

        command = "lightdock3_setup.py %s %s %d %d > test_lightdock.out" % ('1PPE_rec.pdb.H',
                                                                        '1PPE_lig.pdb.H',
                                                                        num_swarms,
                                                                        num_glowworms
                                                                        )
        os.system(command)
        command = "lightdock3.py %s %d -c 1 -s " \
                  "cpydock >> test_lightdock.out" % (self.test_path + 'setup.json',
                                                     steps)
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out',
                           self.test_path + 'swarm_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_1.out',
                           self.test_path + 'swarm_0/gso_1.out')


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

    def test_lightdock_1ppe_10_steps_10_glowworms_1_swarm(self):
        if 'LIGHTDOCK_LONG_TEST' in os.environ and os.environ['LIGHTDOCK_LONG_TEST'] == 'true':
            os.chdir(self.test_path)
            num_swarms = 1
            num_glowworms = 10
            steps = 10
            
            command = "lightdock3_setup.py %s %s %d %d > test_lightdock.out" % ('1PPE_rec.pdb.H',
                                                                            '1PPE_lig.pdb.H',
                                                                            num_swarms,
                                                                            num_glowworms
                                                                            )
            os.system(command)
            command = "lightdock3.py %s %d -c 1 -s " \
                      "cpydock >> test_lightdock.out" % (self.test_path + 'setup.json',
                                                         steps)
            os.system(command)

            assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out',
                               self.test_path + 'swarm_0/gso_0.out')
            assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_10.out',
                               self.test_path + 'swarm_0/gso_10.out')


class TestRegressionPyDockRestraints(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_pydock_restraints/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_pydock_restraints/'
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_lig.pdb'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_2uuy_20_steps_25_glowworms_1_swarm(self):
        os.chdir(self.test_path)
        num_swarms = 4
        num_glowworms = 25
        steps = 20

        command = "lightdock3_setup.py %s %s %d %d -rst %s > test_lightdock.out" % ('2UUY_rec.pdb',
                                                                                '2UUY_lig.pdb',
                                                                                num_swarms,
                                                                                num_glowworms,
                                                                                self.golden_data_path + 'restraints.list'
                                                                                )
        os.system(command)
        command = "lightdock3.py %s %d -c 1 -s " \
                  "cpydock >> test_lightdock.out" % (self.test_path + 'setup.json',
                                                     steps)
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out',
                           self.test_path + 'swarm_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_10.out',
                           self.test_path + 'swarm_0/gso_10.out')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_20.out',
                           self.test_path + 'swarm_0/gso_20.out')
