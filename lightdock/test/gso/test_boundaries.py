"""Tests for Boundaries module"""

from lightdock.gso.boundaries import Boundary
from lightdock.gso.boundaries import BoundingBox


class TestBoundingBox:
    def __init__(self):
        self.dimension1 = Boundary(0, 1)
        self.dimension2 = Boundary(0, 5)

        self.boundaries = [self.dimension1, self.dimension2]

        self.box = BoundingBox(self.boundaries)

    def test_check_dimension_bounding_box(self):
        assert self.box.dimension == 2

    def test_get_boundaries(self):
        expected1 = self.dimension1.clone()
        expected2 = self.dimension2.clone()

        assert expected1 == self.box.get_boundary_of_dimension(0)
        assert expected2 == self.box.get_boundary_of_dimension(1)
        assert expected1 != self.box.get_boundary_of_dimension(1)

    def test_boundary_representation(self):
        assert str(self.dimension1) == "[0, 1]"

    def test_bounding_box_representation(self):
        assert str(self.box) == "[0, 1] [0, 5]"
