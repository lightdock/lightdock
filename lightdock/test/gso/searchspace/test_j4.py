"""Tests for J4 function"""

from lightdock.gso.searchspace.benchmark_ofunctions import J4
from lightdock.gso.coordinates import Coordinates
from nose.tools import assert_almost_equals


class TestJ4:

    def setUp(self):
        self.expected_values = [
            [ 29, 28, 27, 26, 25, 25 ],
            [ 28, 27, 26, 25, 24, 24 ],
            [ 27, 26, 25, 24, 23, 23 ],
            [ 26, 25, 24, 23, 22, 22 ],
            [ 25, 24, 23, 22, 21, 21 ],
            [ 25, 24, 23, 22, 21, 21 ]
            ]

    def tearDown(self):
        pass
    
    def test_compute_J1_matrix(self):
        j4 = J4()
        for i in range(6):
            for j in range(6):
                assert_almost_equals(self.expected_values[i][j], j4(Coordinates([-2.0+j*0.8, -2.0+i*0.8])))
