"""Tests for GSOParameters class"""

from lightdock.gso.parameters import GSOParameters
from lightdock.error.lightdock_errors import GSOParameteresError

from nose.tools import assert_almost_equals #@UnresolvedImport
from nose.tools import raises

import os
import shutil


class TestGSOParameters:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch/'
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
        '/golden_data/parameters/'

    def tearDown(self):
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
    
    def test_read_gso_parameters(self):
        parameteres = GSOParameters()
        
        assert_almost_equals(0.4, parameteres.rho)
        assert_almost_equals(0.6, parameteres.gamma)
        assert_almost_equals(0.08, parameteres.beta)
        assert_almost_equals(5.0, parameteres.initial_luciferin)
        assert_almost_equals(0.2, parameteres.initial_vision_range)
        assert_almost_equals(5.0, parameteres.max_vision_range)
        assert 5 == parameteres.max_neighbors
        
    def test_read_gso_parameters_with_file(self):
        parameteres = GSOParameters(self.golden_data_path + 'glowworm.conf')
        
        assert_almost_equals(0.1, parameteres.rho)
        assert_almost_equals(0.2, parameteres.gamma)
        assert_almost_equals(0.6, parameteres.beta)
        assert_almost_equals(0.3, parameteres.initial_luciferin)
        assert_almost_equals(0.4, parameteres.initial_vision_range)
        assert_almost_equals(0.5, parameteres.max_vision_range)
        assert 7 == parameteres.max_neighbors
        
    @raises(GSOParameteresError)
    def test_read_gso_parameters_wrong_file(self):
        parameteres = GSOParameters(self.golden_data_path + 'no_file.conf')
        assert False
        
    @raises(GSOParameteresError)
    def test_read_gso_parameters_wrong_parameters(self):
        parameteres = GSOParameters(self.golden_data_path + 'wrong_glowworm.conf')
        assert False
