"""J4: staircase function"""

#cython: boundscheck=False
#cython: wraparound=False

from libc.math cimport ceil

def j4(double x, double y):
    return 25.0 - ceil(x) - ceil(y)
