"""Regression tests for testing DFIRE2 scoring function"""

import shutil
import os
import filecmp
from .regression import RegressionTest


class TestRegressionDFIRE2Short(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_dfire2_short/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_dfire2_short/'
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_lig.pdb'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_2uuy_10_steps_25_glowworms_1_swarm(self):
        os.chdir(self.test_path)
        num_swarms = 1
        num_glowworms = 25
        steps = 10

        command = "lightdock3_setup.py %s %s %d %d > test_lightdock.out" % ('2UUY_rec.pdb',
                                                                        '2UUY_lig.pdb',
                                                                        num_swarms,
                                                                        num_glowworms
                                                                        )
        os.system(command)
        command = "lightdock3.py %s %d -c 1 -f %s -s " \
                  "dfire2 >> test_lightdock.out" % (self.test_path + 'setup.json',
                                                   steps,
                                                   self.golden_data_path + 'glowworm.conf')
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out', self.test_path + 'swarm_0/gso_0.out')
        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_10.out', self.test_path + 'swarm_0/gso_10.out')


class TestRegressionDFIRE2Restraints(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_dfire2_restraints/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_dfire2_restraints/'
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '2UUY_lig.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, 'restraints.list'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_2uuy_10_steps_25_glowworms_1_swarm(self):
        os.chdir(self.test_path)
        num_swarms = 4
        num_glowworms = 25
        steps = 5

        command = "lightdock3_setup.py %s %s %d %d -rst %s > test_lightdock.out" % ('2UUY_rec.pdb',
                                                                                '2UUY_lig.pdb',
                                                                                num_swarms,
                                                                                num_glowworms,
                                                                                'restraints.list'
                                                                                )
        os.system(command)
        
        command = "lightdock3.py %s %d -c 1 -s " \
                  "dfire2 >> test_lightdock.out" % (self.test_path + 'setup.json', steps)
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_0.out', self.test_path + 'swarm_0/gso_0.out')
        # Slight differences due to caching
        #assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_5.out', self.test_path + 'swarm_0/gso_5.out')


class TestRegressionDFIRE2Long(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_dfire2_long/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_dfire2_long/'
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_rec.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, '1PPE_lig.pdb'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_1ppe_10_steps_100_glowworms(self):
        if 'LIGHTDOCK_LONG_TEST' in os.environ and os.environ['LIGHTDOCK_LONG_TEST'] == 'true':
            os.chdir(self.test_path)
            num_swarms = 5
            num_glowworms = 50
            steps = 10

            command = "lightdock3_setup.py %s %s %d %d > test_lightdock.out" % ('1PPE_rec.pdb',
                                                                            '1PPE_lig.pdb',
                                                                            num_swarms,
                                                                            num_glowworms
                                                                            )
            os.system(command)
            command = "lightdock3.py %s %d -c 1 -f %s -s " \
                      "dfire2 >> test_lightdock.out" % (self.test_path + 'setup.json',
                                                       steps,
                                                       self.golden_data_path + 'glowworm.conf')
            os.system(command)

            assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_0.dat',
                               self.test_path + 'init/initial_positions_0.dat')
            assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_1.dat',
                               self.test_path + 'init/initial_positions_1.dat')
            assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_2.dat',
                               self.test_path + 'init/initial_positions_2.dat')
            assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_3.dat',
                               self.test_path + 'init/initial_positions_3.dat')
            assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_4.dat',
                               self.test_path + 'init/initial_positions_4.dat')
            assert filecmp.cmp(self.golden_data_path + 'init/starting_positions_0.pdb',
                               self.test_path + 'init/starting_positions_0.pdb')
            assert filecmp.cmp(self.golden_data_path + 'init/starting_positions_1.pdb',
                               self.test_path + 'init/starting_positions_1.pdb')
            assert filecmp.cmp(self.golden_data_path + 'init/starting_positions_2.pdb',
                               self.test_path + 'init/starting_positions_2.pdb')
            assert filecmp.cmp(self.golden_data_path + 'init/starting_positions_3.pdb',
                               self.test_path + 'init/starting_positions_3.pdb')
            assert filecmp.cmp(self.golden_data_path + 'init/starting_positions_4.pdb',
                               self.test_path + 'init/starting_positions_4.pdb')
            assert filecmp.cmp(self.golden_data_path + 'init/cluster_centers.pdb',
                               self.test_path + 'init/cluster_centers.pdb')
            assert filecmp.cmp(self.golden_data_path + 'swarm_0/gso_10.out', self.test_path + 'swarm_0/gso_10.out')
            assert filecmp.cmp(self.golden_data_path + 'swarm_1/gso_10.out', self.test_path + 'swarm_1/gso_10.out')
            assert filecmp.cmp(self.golden_data_path + 'swarm_2/gso_10.out', self.test_path + 'swarm_2/gso_10.out')
            assert filecmp.cmp(self.golden_data_path + 'swarm_3/gso_10.out', self.test_path + 'swarm_3/gso_10.out')
            assert filecmp.cmp(self.golden_data_path + 'swarm_4/gso_10.out', self.test_path + 'swarm_4/gso_10.out')
