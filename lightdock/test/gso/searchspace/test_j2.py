"""Tests for J2 function"""

import pytest
from lightdock.gso.searchspace.benchmark_ofunctions import J2
from lightdock.gso.coordinates import Coordinates


class TestJ2:
    def setup_class(self):
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
                assert self.expected_values[i][j] == pytest.approx(j2(Coordinates([-5.0 + j * 2.5, -5.0 + i * 2.5])))
