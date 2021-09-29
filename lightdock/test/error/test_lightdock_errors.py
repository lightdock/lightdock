"""Tests for errors module"""

from lightdock.error.lightdock_errors import LightDockError, GSOError
from nose.tools import raises


class TestLightDockError:
    def test_create_lightdock_exception(self):
        e = LightDockError("Testing")

        assert str(e) == "[LightDockError] Testing"

    @raises(LightDockError)
    def test_raising_lightdock_exception(self):
        raise LightDockError("Testing")

    @raises(GSOError)
    def test_subclassing_base_exception(self):
        e = GSOError("Testing")

        assert str(e) == "[GSOError] Testing"
        raise e
