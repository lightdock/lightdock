#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

"""J5: plateaus function"""

from libc.math cimport cos

def j5(double x, double y):
    cdef double tmp = cos(x) + cos(y)
    return (tmp > 0.0) - (tmp < 0.0)
