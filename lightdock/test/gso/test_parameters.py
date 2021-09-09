"""Tests for GSOParameters class"""

import os
import shutil
from pathlib import Path
from lightdock.gso.parameters import GSOParameters
from lightdock.error.lightdock_errors import GSOParameteresError
from nose.tools import assert_almost_equals
from nose.tools import raises


class TestGSOParameters:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_gso_parameters"
        self.golden_data_path = self.path / "golden_data" / "parameters"

    def setup(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.mkdir(self.test_path)

    def teardown(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass

    def test_read_gso_parameters(self):
        parameters = GSOParameters()

        assert_almost_equals(0.4, parameters.rho)
        assert_almost_equals(0.6, parameters.gamma)
        assert_almost_equals(0.08, parameters.beta)
        assert_almost_equals(5.0, parameters.initial_luciferin)
        assert_almost_equals(0.2, parameters.initial_vision_range)
        assert_almost_equals(5.0, parameters.max_vision_range)
        assert parameters.max_neighbors == 5

    def test_read_gso_parameters_with_file(self):
        parameters = GSOParameters(self.golden_data_path / "glowworm.conf")

        assert_almost_equals(0.1, parameters.rho)
        assert_almost_equals(0.2, parameters.gamma)
        assert_almost_equals(0.6, parameters.beta)
        assert_almost_equals(0.3, parameters.initial_luciferin)
        assert_almost_equals(0.4, parameters.initial_vision_range)
        assert_almost_equals(0.5, parameters.max_vision_range)
        assert parameters.max_neighbors == 7

    @raises(GSOParameteresError)
    def test_read_gso_parameters_wrong_file(self):
        parameters = GSOParameters(self.golden_data_path / "no_file.conf")
        assert parameters is not None

    @raises(GSOParameteresError)
    def test_read_gso_parameters_wrong_parameters(self):
        parameters = GSOParameters(self.golden_data_path / "wrong_glowworm.conf")
        assert parameters is not None
