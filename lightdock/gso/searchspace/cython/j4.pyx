#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

"""J4: staircase function"""

from libc.math cimport ceil

def j4(double x, double y):
    return 25.0 - ceil(x) - ceil(y)
