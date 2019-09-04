"""Regression tests for testing simulation setup using membrane"""

import shutil
import os
import filecmp

from .regression import RegressionTest


class TestSetupWithMembrane(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_setup_membrane/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_setup_membrane/'
        shutil.copy(os.path.join(self.golden_data_path, 'receptor_membrane.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, 'ligand.pdb'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_setup_with_membrane(self):
        os.chdir(self.test_path)
        num_swarms = 50
        num_glowworms = 10

        command = "lightdock3_setup.py %s %s %d %d --noh --noxt -membrane > test_lightdock.out" % ('receptor_membrane.pdb',
                                                                                               'ligand.pdb',
                                                                                               num_swarms,
                                                                                               num_glowworms
                                                                                               )
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_0.dat',
                           self.test_path + 'init/initial_positions_0.dat')
        assert filecmp.cmp(self.golden_data_path + 'init/cluster_centers.pdb',
                           self.test_path + 'init/cluster_centers.pdb')
        assert filecmp.cmp(self.golden_data_path + 'setup.json',
                           self.test_path + 'setup.json')


class TestSetupWithMembraneAndRestraints(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_setup_membrane_restraints/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/regression_setup_membrane_restraints/'
        shutil.copy(os.path.join(self.golden_data_path, 'receptor_membrane.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, 'ligand.pdb'), self.test_path)
        shutil.copy(os.path.join(self.golden_data_path, 'restraints.list'), self.test_path)

    def teardown(self):
        self.clean_test_path()

    def test_lightdock_setup_with_membrane_and_restraints(self):
        os.chdir(self.test_path)
        num_swarms = 50
        num_glowworms = 10

        command = "lightdock3_setup.py %s %s %d %d --noh --noxt -membrane -rst %s > test_lightdock.out" % ('receptor_membrane.pdb',
                                                                                                       'ligand.pdb',
                                                                                                       num_swarms,
                                                                                                       num_glowworms,
                                                                                                       'restraints.list'
                                                                                                       )
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'init/initial_positions_0.dat',
                           self.test_path + 'init/initial_positions_0.dat')
        assert filecmp.cmp(self.golden_data_path + 'init/cluster_centers.pdb',
                           self.test_path + 'init/cluster_centers.pdb')
        assert filecmp.cmp(self.golden_data_path + 'setup.json',
                           self.test_path + 'setup.json')
