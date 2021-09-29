"""Tests for Glowworm class using J1 function"""

from math import sqrt
from nose.tools import assert_almost_equals
from lightdock.gso.parameters import GSOParameters
from lightdock.gso.searchspace.benchmark_ofunctions import J1
from lightdock.gso.searchspace.landscape import LandscapePosition
from lightdock.gso.coordinates import Coordinates
from lightdock.gso.glowworm import Glowworm


class TestGlowwormWithJ1:
    def __init__(self):
        self.gso_parameters = GSOParameters()
        self.objective_function = J1()
        self.landscape_position1 = [
            LandscapePosition(self.objective_function, Coordinates([0.0, 0.0]))
        ]
        self.landscape_position2 = [
            LandscapePosition(self.objective_function, Coordinates([1.0, 1.0]))
        ]
        self.landscape_position3 = [
            LandscapePosition(self.objective_function, Coordinates([2.0, 2.0]))
        ]
        self.landscape_position4 = [
            LandscapePosition(self.objective_function, Coordinates([0.5, 0.5]))
        ]
        self.landscape_position5 = [
            LandscapePosition(self.objective_function, Coordinates([5.0, 5.0]))
        ]
        self.landscape_position6 = [
            LandscapePosition(self.objective_function, Coordinates([0.0, 1.0]))
        ]
        self.landscape_position7 = [
            LandscapePosition(self.objective_function, Coordinates([0.1, 0.0]))
        ]
        self.landscape_position8 = [
            LandscapePosition(self.objective_function, Coordinates([0.0, 0.1]))
        ]
        self.landscape_position9 = [
            LandscapePosition(self.objective_function, Coordinates([0.2, 0.0]))
        ]

    def test_create_glowworm(self):
        glowworm = Glowworm(self.landscape_position1, self.gso_parameters)
        assert_almost_equals(0.4, glowworm.rho)

    def test_compute_luciferin(self):
        glowworm = Glowworm(self.landscape_position1, self.gso_parameters)
        assert_almost_equals(
            (1.0 - 0.4) * 5.0 + 0.6 * 0.9810118431238463, glowworm.compute_luciferin()
        )

    def check_compare_glowworms(self):
        glowworm1 = Glowworm(self.landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(self.landscape_position2, self.gso_parameters)

        assert glowworm1 == glowworm1
        assert glowworm1 != glowworm2

    def test_distance(self):
        glowworm1 = Glowworm(self.landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(self.landscape_position3, self.gso_parameters)

        assert_almost_equals(sqrt(8.0), glowworm1.distance(glowworm2))

    def test_distance2(self):
        glowworm1 = Glowworm(self.landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(self.landscape_position3, self.gso_parameters)

        assert_almost_equals(8.0, glowworm1.distance2(glowworm2))

    def test_search_my_neighbors(self):
        glowworm1 = Glowworm(self.landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(self.landscape_position4, self.gso_parameters)
        glowworm3 = Glowworm(self.landscape_position5, self.gso_parameters)
        glowworm1.luciferin = 3.0
        glowworm1.vision_range = 1.0

        glowworms = [glowworm1, glowworm2, glowworm3]

        glowworm1.search_neighbors(glowworms)

        assert len(glowworm1.neighbors) == 1
        assert glowworm1.neighbors[0] == glowworm2

    def test_search_my_neighbors_when_distance_equal_range(self):
        glowworm1 = Glowworm(self.landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(self.landscape_position4, self.gso_parameters)
        glowworm3 = Glowworm(self.landscape_position5, self.gso_parameters)
        glowworm4 = Glowworm(self.landscape_position6, self.gso_parameters)
        glowworm1.luciferin = 3.0
        glowworm1.vision_range = 1.0

        glowworms = [glowworm1, glowworm2, glowworm3, glowworm4]

        glowworm1.search_neighbors(glowworms)

        assert len(glowworm1.neighbors) == 1
        assert glowworm1.neighbors[0] == glowworm2

    def test_search_my_neighbors_when_luciferin_equal_to_my_luciferin(self):
        glowworm1 = Glowworm(self.landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(self.landscape_position4, self.gso_parameters)
        glowworm3 = Glowworm(self.landscape_position5, self.gso_parameters)
        glowworm4 = Glowworm(self.landscape_position6, self.gso_parameters)
        glowworm1.luciferin = 3.0
        glowworm1.vision_range = 1.0
        glowworm4.luciferin = 3.0
        glowworm4.vision_range = 1.0

        glowworms = [glowworm1, glowworm2, glowworm3, glowworm4]

        glowworm1.search_neighbors(glowworms)

        assert len(glowworm1.neighbors) == 1
        assert glowworm1.neighbors[0] == glowworm2

    def test_compute_probability_moving_toward_neighbor(self):
        glowworm1 = Glowworm(self.landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(self.landscape_position7, self.gso_parameters)
        glowworm3 = Glowworm(self.landscape_position8, self.gso_parameters)
        glowworm4 = Glowworm(self.landscape_position9, self.gso_parameters)
        glowworm1.luciferin = 3.0
        glowworm1.vision_range = 1.0
        glowworm2.luciferin = 4.0
        glowworm2.vision_range = 1.0
        glowworm3.luciferin = 5.0
        glowworm3.vision_range = 1.0
        glowworm4.luciferin = 6.0
        glowworm4.vision_range = 1.0

        glowworm1.neighbors = [glowworm2, glowworm3, glowworm4]

        glowworm1.compute_probability_moving_toward_neighbor()

        assert_almost_equals(1.0 / 6.0, glowworm1.probabilities[0])
        assert_almost_equals(2.0 / 6.0, glowworm1.probabilities[1])
        assert_almost_equals(3.0 / 6.0, glowworm1.probabilities[2])

    def test_select_a_glowworm_by_random_number(self):
        glowworm1 = Glowworm(self.landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(self.landscape_position7, self.gso_parameters)
        glowworm3 = Glowworm(self.landscape_position8, self.gso_parameters)
        glowworm4 = Glowworm(self.landscape_position9, self.gso_parameters)
        glowworm1.luciferin = 3.0
        glowworm1.vision_range = 1.0
        glowworm2.luciferin = 4.0
        glowworm2.vision_range = 1.0
        glowworm3.luciferin = 5.0
        glowworm3.vision_range = 1.0
        glowworm4.luciferin = 6.0
        glowworm4.vision_range = 1.0

        glowworm1.neighbors = [glowworm2, glowworm3, glowworm4]

        glowworm1.probabilities = [1.0 / 6.0, 2.0 / 6.0, 3.0 / 6.0]

        assert glowworm4 == glowworm1.select_random_neighbor(0.55)
        assert glowworm2 == glowworm1.select_random_neighbor(0.1)
        assert glowworm3 == glowworm1.select_random_neighbor(0.2)

    def test_update_visio_range_when_number_neighbors_is_less_than_maximum_neighbors(
        self,
    ):
        glowworm1 = Glowworm(self.landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(self.landscape_position7, self.gso_parameters)
        glowworm3 = Glowworm(self.landscape_position8, self.gso_parameters)
        glowworm1.luciferin = 3.0
        glowworm1.vision_range = 2.0
        glowworm1.max_vision_range = 1.0
        glowworm2.luciferin = 4.0
        glowworm2.vision_range = 2.0
        glowworm2.max_vision_range = 1.0
        glowworm3.luciferin = 5.0
        glowworm3.vision_range = 2.0
        glowworm3.max_vision_range = 1.0

        glowworm1.neighbors = [glowworm2, glowworm3]

        glowworm1.update_vision_range()

        assert_almost_equals(1.0, glowworm1.vision_range)

    def test_update_visio_range_when_number_neighbors_is_greater_than_maximum_neighbors(
        self,
    ):
        glowworm1 = Glowworm(self.landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(self.landscape_position7, self.gso_parameters)
        glowworm1.luciferin = 5.0
        glowworm1.vision_range = 3.0
        glowworm1.max_vision_range = 3.0
        glowworm2.luciferin = 5.0
        glowworm2.vision_range = 3.0
        glowworm2.max_vision_range = 3.0

        for _ in range(23):
            glowworm1.neighbors.append(glowworm2)

        assert len(glowworm1.neighbors) == 23

        glowworm1.update_vision_range()

        assert_almost_equals(1.5600000000000001, glowworm1.vision_range)

    def test_move_different_coordinates(self):
        landscape_position1 = [
            LandscapePosition(self.objective_function, Coordinates([1.0, 2.0]))
        ]
        landscape_position2 = [
            LandscapePosition(self.objective_function, Coordinates([0.0, 1.0]))
        ]

        glowworm1 = Glowworm(landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(landscape_position2, self.gso_parameters)
        glowworm1.luciferin = 4.0
        glowworm1.vision_range = 2.0
        glowworm1.max_vision_range = 1.0
        glowworm2.luciferin = 4.0
        glowworm2.vision_range = 2.0
        glowworm2.max_vision_range = 1.0

        glowworm1.move(glowworm2)

        expected = LandscapePosition(
            self.objective_function,
            Coordinates([-0.03 / sqrt(2.0) + 1.0, -0.03 / sqrt(2.0) + 2.0]),
        )

        assert expected == glowworm1.landscape_positions[0]

    def test_move_same_coordinates(self):
        landscape_position1 = [
            LandscapePosition(self.objective_function, Coordinates([1.0, 2.0]))
        ]
        landscape_position2 = [
            LandscapePosition(self.objective_function, Coordinates([1.0, 2.0]))
        ]

        glowworm1 = Glowworm(landscape_position1, self.gso_parameters)
        glowworm2 = Glowworm(landscape_position2, self.gso_parameters)
        glowworm1.luciferin = 4.0
        glowworm1.vision_range = 2.0
        glowworm1.max_vision_range = 1.0
        glowworm2.luciferin = 4.0
        glowworm2.vision_range = 2.0
        glowworm2.max_vision_range = 1.0

        glowworm1.move(glowworm2)

        expected = LandscapePosition(self.objective_function, Coordinates([1.0, 2.0]))

        assert expected == glowworm1.landscape_positions[0]
