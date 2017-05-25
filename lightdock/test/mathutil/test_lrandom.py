"""Tests for RandomNumberGenerator interface"""

from lightdock.mathutil.lrandom import RandomNumberGenerator
from nose.tools import raises


class TestRandomNumberGenerator:

    @raises(NotImplementedError)
    def test_call_interface(self):
        gen = RandomNumberGenerator()
        gen()
