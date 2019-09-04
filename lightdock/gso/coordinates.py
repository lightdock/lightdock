"""Glowworm's position in a given landscape."""

import os
from lightdock.mathutil.cython.cutil import float_equals as cfloat_equals
from lightdock.mathutil.cython.cutil import norm as cnorm
from lightdock.mathutil.cython.cutil import sum_of_squares as csum_of_squares
from lightdock.mathutil.cython.cutil import sum_of_square_difference as csum_of_square_difference
from lightdock.error.lightdock_errors import GSOCoordinatesError


class Coordinates(object):
    """A Coordinates object is an array of float numbers with size equal to the
    dimensions of the solution space used in the objective function."""
    def __init__(self, values):
        self._values = values
        self.dimension = len(self._values)
        
    def __getitem__(self, index):
        """Gets the item at index"""
        return self._values[index]
    
    def __setitem__(self, index, value):
        """Sets a value at index"""
        self._values[index] = value
    
    def __eq__(self, other):
        """Compares for equality"""
        if self.dimension == other.dimension:
            for c1, c2 in zip(self._values, other._values):
                if not cfloat_equals(c1, c2):
                    return False
            return True
        else:
            return False
    
    def __ne__(self, other):
        """Compares for unequality"""
        return not self.__eq__(other)
    
    def clone(self):
        """Get a copy of the current coordinate"""
        return Coordinates(self._values*1)
    
    def __add__(self, other):
        """Adds two coordinates"""
        return Coordinates([sum(pair) for pair in zip(self._values, other._values)])

    def __iadd__(self, other):
        """Adds and assigns another coordinate"""
        for i in range(self.dimension):
            self._values[i] += other._values[i]
        return self
    
    def __sub__(self, other):
        """Subtracts two coordinates"""
        return Coordinates([c1-c2 for c1, c2 in zip(self._values, other._values)])
  
    def __isub__(self, other):
        """Subtracts and assigns another coordinate"""
        for i in range(self.dimension):
            self._values[i] -= other._values[i]
        return self
  
    def __imul__(self, scalar):
        """Multiplies a coordinate by a scalar"""
        self._values = [v*scalar for v in self._values]
        return self

    def __mul__(self, scalar):
        """Multiplies a coordinate by a scalar"""
        values = [v*scalar for v in self._values]
        return Coordinates(values)
    
    def norm(self):
        """Calculates the norm of a coordinate"""
        return cnorm(self._values)
    
    def distance(self, other):
        """Distance between two coordinates"""
        return (self - other).norm()
    
    def distance2(self, other):
        """Square distance between two coordinates"""
        return csum_of_square_difference(self._values, other._values)
    
    def sum_of_squares(self):
        """Calculates the sum of squares of these coordinates"""
        return csum_of_squares(self._values)
    
    def move(self, other, step=1.0):
        """Move from one coordinate to another a given step"""
        if self != other:
            delta_x = other - self
            delta_x *= step/delta_x.norm()
            self += delta_x
        return self
    
    def __repr__(self):
        """Coordinate representation"""
        # Path for equivalent representation from Python 2.7
        coord = ["{:.12g}".format(f) for f in self._values]
        return "(%s)" % ', '.join(["{}.0".format(f) if '.' not in f else f for f in coord])

    def __len__(self):
        return self.dimension


class CoordinatesFileReader(object):
    """Reads spatial coordinates from a given file"""
    def __init__(self, dimension):
        self.dimension = dimension
        
    def get_coordinates_from_file(self, coordinates_file):
        """Parses and creates coordinates from coordinates_file"""
        coordinates = []
        try:
            input_file = open(coordinates_file)
            for line in input_file:
                values = line.rstrip(os.linesep).split()
                if len(values) == self.dimension:
                    coordinates.append(Coordinates([float(value) for value in values]))
                else:
                    raise Exception("dimension %d does not correspond with values in line %s" % (self.dimension,
                                                                                                 line))
        except Exception as e:
            raise GSOCoordinatesError("Error reading coordinates from file: %s" % str(e))
        
        return coordinates
