
import numpy as np


class SpacePoints(object):
    """A collection of spatial points"""
    def __init__(self, coordinates):
        self.coordinates = np.array(coordinates)

    def clone(self):
        return SpacePoints(self.coordinates.copy())

    def translate(self, vector):
        """Translates coordinates based on vector"""
        self.coordinates += vector

    def rotate(self, q):
        """Rotates coordinates using a quaternion q"""
        for i in range(self.coordinates.shape[0]):
            self.coordinates[i, ...] = q.rotate(self.coordinates[i, ...])

    def __getitem__(self, item):
        return self.coordinates[item]

    def __setitem__(self, index, item):
        self.coordinates[index] = item

    def __iter__(self):
        for coordinates in self.coordinates:
            yield coordinates

    def __len__(self):
        return self.coordinates.shape[0]

    def __eq__(self, other):
        return np.allclose(self.coordinates, other.coordinates)

    def __ne__(self, other):
        return not (self == other)

    def __sub__(self, other):
        return self.coordinates - other.coordinates

    def __str__(self):
        return str(self.coordinates)

    def shape(self):
        return self.coordinates.shape
