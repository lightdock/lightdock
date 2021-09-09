"""Tests for J2 function"""

from lightdock.gso.searchspace.benchmark_ofunctions import J2
from lightdock.gso.coordinates import Coordinates
from nose.tools import assert_almost_equals


class TestJ2:
    def __init__(self):
        self.expected_values = [
            [50, 51.25, 25, 51.25, 50],
            [51.25, 52.5, 26.25, 52.5, 51.25],
            [25, 26.25, 0, 26.25, 25],
            [51.25, 52.5, 26.25, 52.5, 51.25],
            [50, 51.25, 25, 51.25, 50],
        ]

    def test_compute_J1_matrix(self):
        j2 = J2()
        for i in range(5):
            for j in range(5):
                assert_almost_equals(
                    self.expected_values[i][j],
                    j2(Coordinates([-5.0 + j * 2.5, -5.0 + i * 2.5])),
                )
