#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

"""J1: peaks function"""

from libc.math cimport sqrt, exp, pow

def j1(double x, double y):
    return 3.0 * (1-x)*(1-x) * exp(-(x*x + (y+1)*(y+1))) \
        - 10.0 * (x/5.0 - pow(x,3) - pow(y,5)) * exp(-(x*x + y*y)) \
        - 1/3.0 * exp(-((x+1)*(x+1) + y*y))
