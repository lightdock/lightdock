"""Tests for LandscapePosition class"""

from math import sqrt
from lightdock.gso.coordinates import Coordinates
from lightdock.gso.searchspace.landscape import LandscapePosition
from lightdock.gso.searchspace.benchmark_ofunctions import J1
from nose.tools import assert_almost_equals


class TestLandscapePosition:
    def __init__(self):
        self.objective_function = J1()
        self.coordinates1 = Coordinates([0.0, 0.0])
        self.coordinates2 = Coordinates([1.0, 1.0])
        self.coordinates3 = Coordinates([5.0, -5.0])
        self.coordinates4 = Coordinates([2.0, 2.0])
        self.coordinates5 = Coordinates([-1.0, -1.0])

    def test_LandscapePosition_with_J1(self):
        landscape_point = LandscapePosition(self.objective_function, self.coordinates1)

        assert_almost_equals(
            0.9810118431238463, landscape_point.evaluate_objective_function()
        )

    def test_LandscapePosition_with_J1_matrix(self):
        expected_values = [
            [
                6.671280296717448e-05,
                0.004156270134668844,
                -0.2449540405743497,
                -0.02975999770667323,
                -5.864187872589536e-06,
            ],
            [
                -0.0004517936594085711,
                0.3265409898164098,
                -5.680276001896266,
                -0.4404918548935705,
                0.003599520786098113,
            ],
            [
                -0.03650620461319554,
                -2.773610019456342,
                0.9810118431238463,
                3.269463326439044,
                0.03312494992430833,
            ],
            [
                -0.003078234252739988,
                0.4784411467828077,
                7.996620241631349,
                1.185275846660814,
                0.0044245231361739,
            ],
            [
                3.223535961269275e-05,
                0.0311759440751708,
                0.2998710282262348,
                0.03200763718576478,
                4.102972745826762e-05,
            ],
        ]

        for i in range(5):
            for j in range(5):
                coord = Coordinates([-3.0 + i * 1.5, -3.0 + j * 1.5])
                landscape_point = LandscapePosition(self.objective_function, coord)
                assert_almost_equals(
                    expected_values[j][i], landscape_point.evaluate_objective_function()
                )

    def test_equals(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates1)
        landscape_point2 = LandscapePosition(self.objective_function, self.coordinates2)
        landscape_point3 = LandscapePosition(self.objective_function, self.coordinates2)

        assert landscape_point1 != landscape_point2
        assert landscape_point2 == landscape_point3

        assert landscape_point1 != landscape_point2
        assert landscape_point2 == landscape_point3

    def test_clone(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates1)

        landscape_point2 = landscape_point1.clone()

        expected = LandscapePosition(self.objective_function, self.coordinates1)

        assert expected == landscape_point2

        landscape_point2.coordinates[0] = 1.0

        assert landscape_point1 != landscape_point2

    def test_addition_and_assigment(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates1)
        landscape_point2 = LandscapePosition(self.objective_function, self.coordinates2)

        expected = LandscapePosition(self.objective_function, self.coordinates2)

        landscape_point1 += landscape_point2

        assert expected == landscape_point1

    def test_self_addition_and_assigment(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates2)

        expected = LandscapePosition(self.objective_function, self.coordinates4)

        landscape_point1 += landscape_point1

        assert expected == landscape_point1

    def test_addition(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates1)
        landscape_point2 = LandscapePosition(self.objective_function, self.coordinates2)

        expected = LandscapePosition(self.objective_function, self.coordinates2)

        assert expected == (landscape_point1 + landscape_point2)

    def test_subtraction_and_assigment(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates1)
        landscape_point2 = LandscapePosition(self.objective_function, self.coordinates2)

        expected = LandscapePosition(self.objective_function, self.coordinates5)

        landscape_point1 -= landscape_point2

        assert expected == landscape_point1

    def test_subtraction(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates1)
        landscape_point2 = LandscapePosition(self.objective_function, self.coordinates2)

        expected = LandscapePosition(self.objective_function, self.coordinates5)

        assert expected == (landscape_point1 - landscape_point2)

    def test_multiply_by_scalar(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates2)

        expected = LandscapePosition(self.objective_function, self.coordinates4)

        assert expected == (landscape_point1 * 2.0)

    def test_distance_J1(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates2)
        landscape_point2 = LandscapePosition(self.objective_function, self.coordinates1)

        assert_almost_equals(sqrt(2.0), landscape_point1.distance(landscape_point2))

    def test_distance2_J1(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates2)
        landscape_point2 = LandscapePosition(self.objective_function, self.coordinates1)

        assert_almost_equals(2.0, landscape_point1.distance2(landscape_point2))

    def test_norm(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates2)

        assert_almost_equals(sqrt(2.0), landscape_point1.norm())

    def test_move_using_default_step(self):
        landscape_point1 = LandscapePosition(self.objective_function, self.coordinates1)
        landscape_point2 = LandscapePosition(self.objective_function, self.coordinates2)

        expected = LandscapePosition(
            self.objective_function,
            Coordinates([0.021213203435596423, 0.021213203435596423]),
        )

        landscape_point1.move(landscape_point2)

        assert expected == landscape_point1

    def test_move_using_other_step(self):
        landscape_point1 = LandscapePosition(
            self.objective_function, self.coordinates1, 1.0
        )
        landscape_point2 = LandscapePosition(self.objective_function, self.coordinates2)

        expected = LandscapePosition(
            self.objective_function,
            Coordinates([0.70710678118654746, 0.70710678118654746]),
        )

        landscape_point1.move(landscape_point2)

        assert expected == landscape_point1
