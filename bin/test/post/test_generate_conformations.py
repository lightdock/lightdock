"""Test for generate_conformations post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestGenerateConformations(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_generate_conformations/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/generate_conformations/'

    def teardown(self):
        self.clean_test_path()

    def test_generate_conformations(self):
        os.chdir(self.test_path)
        num_conformations = 2
        shutil.copyfile(self.golden_data_path + '1PPE_rec.pdb', self.test_path + '1PPE_rec.pdb')
        shutil.copyfile(self.golden_data_path + '1PPE_lig.pdb', self.test_path + '1PPE_lig.pdb')
        shutil.copyfile(self.golden_data_path + 'lightdock_1PPE_rec.pdb', self.test_path + 'lightdock_1PPE_rec.pdb')
        shutil.copyfile(self.golden_data_path + 'lightdock_1PPE_lig.pdb', self.test_path + 'lightdock_1PPE_lig.pdb')
        shutil.copyfile(self.golden_data_path + 'gso_1.out', self.test_path + 'gso_1.out')
        command = "lgd_generate_conformations.py %s %s %s %d > test.out" % (self.test_path + '1PPE_rec.pdb',
                                                                            self.test_path + '1PPE_lig.pdb',
                                                                            self.test_path + 'gso_1.out',
                                                                            num_conformations
                                                                            )
        os.system(command)

        assert filecmp.cmp(self.golden_data_path + 'lightdock_0.pdb', self.test_path + 'lightdock_0.pdb')
        assert filecmp.cmp(self.golden_data_path + 'lightdock_1.pdb', self.test_path + 'lightdock_1.pdb')

