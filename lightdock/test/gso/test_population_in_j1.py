"""Tests for Population class using J1 function"""

from lightdock.gso.parameters import GSOParameters
from lightdock.gso.searchspace.benchmark_ofunctions import J1
from lightdock.gso.searchspace.landscape import LandscapePosition
from lightdock.gso.coordinates import Coordinates
from lightdock.gso.swarm import Swarm


class TestSwarmInJ1:
    def __init__(self):
        self.gso_parameters = GSOParameters()
        self.objective_function = J1()
        self.landscape_position1 = LandscapePosition(
            self.objective_function, Coordinates([0.0, 0.0])
        )
        self.landscape_position2 = LandscapePosition(
            self.objective_function, Coordinates([1.0, 1.0])
        )
        self.positions = [[self.landscape_position1, self.landscape_position2]]

    def test_create_swarm(self):
        swarm = Swarm(self.positions, self.gso_parameters)

        expected = "#Coordinates  Luciferin  Neighbor's number  Vision Range  Scoring\n(0.0, 0.0)   5.00000000  0 0.200   0.00000000\n(1.0, 1.0)   5.00000000  0 0.200   0.00000000\n"

        assert swarm.get_size() == 2
        assert expected == str(swarm)
