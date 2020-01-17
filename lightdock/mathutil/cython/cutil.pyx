#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

from libc.math cimport sqrt
from libc.math cimport abs

def sum_of_squares(values):
    """
    Calculates the sum of squares of each value in values
    """
    cdef double sum
    sum = 0.0
    for value in values:
        sum += value*value
    return sum

def norm(values):
    """
    Calculates the norm of the given values
    """
    return sqrt(sum_of_squares(values))

def float_equals(double a, double b, precision=1e-7):
    """
    Compares two floats for equality given precision
    """
    return abs(a-b) < precision

def sum_of_square_difference(v1, v2):
    """
    Calculates the sum of the square differences of the components of 
    v1 and v2 vectors.
    """
    cdef double sum = 0.0
    for c1, c2 in zip(v1, v2):
        t = c1-c2
        sum += t*t
    return sum

def distance(double x1, double y1, double z1, double x2, double y2, double z2):
    """
    Calculates the distance between point 1 and point 2.
    """
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))

def distance2(double x1, double y1, double z1, double x2, double y2, double z2):
    """
    Calculates the distance^2 between point 1 and point 2.
    """
    return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)
