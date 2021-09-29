"""Represents GSO search space paradigm"""

from lightdock.gso.searchspace.cython.j1 import j1
from lightdock.gso.searchspace.cython.j2 import j2
from lightdock.gso.searchspace.cython.j3 import j3
from lightdock.gso.searchspace.cython.j4 import j4
from lightdock.gso.searchspace.cython.j5 import j5
from lightdock.gso.searchspace.ofunction import ObjectiveFunction


class J1(ObjectiveFunction):
    """Peaks function

    The Peaks function is a function of two variables, obtained by
    translating and scaling Gaussian distributions (Reutskiy and Chen 2006).
    """

    def __call__(self, coordinates):
        return j1(coordinates[0], coordinates[1])


class J2(ObjectiveFunction):
    """Rastrigin's function

    Function with large number of local minima and maxima (Torn and
    Zilinskas 1989).
    """

    def __call__(self, coordinates):
        return j2(coordinates[0], coordinates[1])


class J3(ObjectiveFunction):
    """Circles function

    Function that contains multiple concentric circles as the regions of
    local maxima (Muller et al. 2002). Note the infinite-peaks case.
    """

    def __call__(self, coordinates):
        return j3(coordinates[0], coordinates[1])


class J4(ObjectiveFunction):
    """Staircase function

    The Staircase function contains a series of stairs (Dreo and Siarry 2004).
    """

    def __call__(self, coordinates):
        return j4(coordinates[0], coordinates[1])


class J5(ObjectiveFunction):
    """Plateaus function

    Function containing multiple plateaus with equal objective function values
    (Singh et al. 2004).
    """

    def __call__(self, coordinates):
        return j5(coordinates[0], coordinates[1])
