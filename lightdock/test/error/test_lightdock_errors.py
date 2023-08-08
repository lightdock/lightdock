"""Tests for errors module"""

import pytest
from lightdock.error.lightdock_errors import LightDockError, GSOError


class TestLightDockError:
    def test_create_lightdock_exception(self):
        e = LightDockError("Testing")

        assert str(e) == "[LightDockError] Testing"

    def test_raising_lightdock_exception(self):
        with pytest.raises(LightDockError):
            raise LightDockError("Testing")

    def test_subclassing_base_exception(self):
        with pytest.raises(GSOError):
            e = GSOError("Testing")

            assert str(e) == "[GSOError] Testing"
            raise e
