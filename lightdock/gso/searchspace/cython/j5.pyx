"""J5: plateaus function"""

#cython: boundscheck=False
#cython: wraparound=False

from libc.math cimport cos

def j5(double x, double y):
    cdef double tmp = cos(x) + cos(y)
    return (tmp > 0.0) - (tmp < 0.0)
