"""Tests for C-util class"""

from lightdock.mathutil.cython.cutil import sum_of_squares
from nose.tools import assert_almost_equals #@UnresolvedImport


class TestCutil:

    def setUp(self):
        self.values = list(range(1000))

    def tearDown(self):
        pass
    
    def test_sum_of_squares(self):
        assert_almost_equals(332833500.0, sum_of_squares(self.values))
