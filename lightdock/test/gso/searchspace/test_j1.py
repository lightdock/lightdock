"""Tests for J1 function"""

from lightdock.gso.searchspace.benchmark_ofunctions import J1
from lightdock.gso.coordinates import Coordinates
from nose.tools import assert_almost_equals


class TestJ1:

    def setUp(self):
        self.expected_values = [
            [6.671280296717448e-05, 0.004156270134668844, -0.2449540405743497, -0.02975999770667323, -5.864187872589536e-06],
            [-0.0004517936594085711, 0.3265409898164098, -5.680276001896266, -0.4404918548935705, 0.003599520786098113],
            [-0.03650620461319554, -2.773610019456342, 0.9810118431238463, 3.269463326439044, 0.03312494992430833],
            [-0.003078234252739988, 0.4784411467828077, 7.996620241631349, 1.185275846660814, 0.0044245231361739],
            [3.223535961269275e-05, 0.0311759440751708, 0.2998710282262348, 0.03200763718576478, 4.102972745826762e-05]
            ]

    def tearDown(self):
        pass
    
    def test_compute_J1_matrix(self):
        j1 = J1()
        for i in range(5):
            for j in range(5):
                assert_almost_equals(self.expected_values[i][j], j1(Coordinates([-3.0+j*1.5, -3.0+i*1.5])))
