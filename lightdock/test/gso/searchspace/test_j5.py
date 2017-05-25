"""Tests for J5 function"""

from lightdock.gso.searchspace.benchmark_ofunctions import J5
from nose.tools import assert_equals
from lightdock.gso.coordinates import Coordinates


class TestJ5:

    def setUp(self):
        self.expected_values = [
             [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
             [1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1],
             [1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1],
             [1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1],
             [1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1],
             [1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1],
             [1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1],
             [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
             [1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1],
             [1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1],
             [1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1],
             [1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1],
             [1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1],
             [1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1],
             [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            ]

    def tearDown(self):
        pass
    
    def test_compute_J1_matrix(self):
        j5 = J5()
        for i in range(15):
            for j in range(15):
                assert_equals(self.expected_values[i][j], j5(Coordinates([-6.24+j*0.89143, -6.24+i*0.89143])))
