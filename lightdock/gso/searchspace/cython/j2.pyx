#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

"""J2: Rastrigin's function"""

from libc.math cimport cos, M_PI

def j2(double x, double y):
    return 20.0 + (x*x - 10.0*cos(2*M_PI*x) + y*y - 10.0*cos(2*M_PI*y))
