"""Tests for C-util class"""

import pytest
from lightdock.mathutil.cython.cutil import sum_of_squares


class TestCutil:
    def setup_class(self):
        self.values = list(range(1000))

    def test_sum_of_squares(self):
        assert 332833500.0 == pytest.approx(sum_of_squares(self.values))
