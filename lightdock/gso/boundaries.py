"""Module for representing dimension boundaries"""

from lightdock.mathutil.cython.cutil import float_equals as cfloat_equals


class Boundary(object):
    """Represents a boundary for a given dimension"""

    def __init__(self, lower_limit, upper_limit):
        """Creates a boundary of a dimension with lower and upper limits"""
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit

    def clone(self):
        """Gets a copy of this boundary"""
        return Boundary(self.lower_limit, self.upper_limit)

    def __eq__(self, other):
        """Compares for equality two boundaries"""
        return cfloat_equals(self.lower_limit, other.lower_limit) and cfloat_equals(
            self.upper_limit, other.upper_limit
        )

    def __ne__(self, other):
        """Compares for unequality two boundaries"""
        return not self.__eq__(other)

    def __repr__(self):
        return "[%s, %s]" % (self.lower_limit, self.upper_limit)


class BoundingBox(object):
    """Represents a set of boundaries to apply to each dimension of a given space"""

    def __init__(self, boundaries):
        """Creates a bounding box with a Boundary for each dimension"""
        self.boundaries = boundaries
        self.dimension = len(self.boundaries)

    def get_boundary_of_dimension(self, dimension_index):
        """Gets the Boundary of the dimension with dimension_index"""
        return self.boundaries[dimension_index]

    def __repr__(self):
        return " ".join([str(b) for b in self.boundaries])
