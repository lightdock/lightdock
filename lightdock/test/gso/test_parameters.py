"""Tests for GSOParameters class"""

import pytest
from pathlib import Path
from lightdock.gso.parameters import GSOParameters
from lightdock.error.lightdock_errors import GSOParameteresError


class TestGSOParameters:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data" / "parameters"

    def test_read_gso_parameters(self):
        parameters = GSOParameters()

        assert 0.4 == pytest.approx(parameters.rho)
        assert 0.6 == pytest.approx(parameters.gamma)
        assert 0.08 == pytest.approx(parameters.beta)
        assert 5.0 == pytest.approx(parameters.initial_luciferin)
        assert 0.2 == pytest.approx(parameters.initial_vision_range)
        assert 5.0 == pytest.approx(parameters.max_vision_range)
        assert parameters.max_neighbors == 5

    def test_read_gso_parameters_with_file(self):
        parameters = GSOParameters(self.golden_data_path / "glowworm.conf")

        assert 0.1 == pytest.approx(parameters.rho)
        assert 0.2 == pytest.approx(parameters.gamma)
        assert 0.6 == pytest.approx(parameters.beta)
        assert 0.3 == pytest.approx(parameters.initial_luciferin)
        assert 0.4 == pytest.approx(parameters.initial_vision_range)
        assert 0.5 == pytest.approx(parameters.max_vision_range)
        assert parameters.max_neighbors == 7

    def test_read_gso_parameters_wrong_file(self):
        with pytest.raises(GSOParameteresError):
            parameters = GSOParameters(self.golden_data_path / "no_file.conf")
            assert parameters is not None

    def test_read_gso_parameters_wrong_parameters(self):
        with pytest.raises(GSOParameteresError):
            parameters = GSOParameters(self.golden_data_path / "wrong_glowworm.conf")
            assert parameters is not None
