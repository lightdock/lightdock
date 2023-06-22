from typing import Iterable


def sum_of_squares(values: Iterable[float]) -> float:
    """
    Calculates the sum of squares of each value in values
    """
    ...


def norm(values: Iterable[float]) -> float:
    """
    Calculates the norm of the given values
    """
    ...


def float_equals(a: float, b: float, precision: float = 1e-7) -> bool:
    """
    Compares two floats for equality given precision
    """
    ...


def sum_of_square_difference(v1: Iterable[float], v2: Iterable[float]) -> float:
    """
    Calculates the sum of the square differences of the components of 
    v1 and v2 vectors.
    """
    ...


def distance(x1: float, y1: float, z1: float, x2: float, y2: float, z2: float) -> float:
    """
    Calculates the distance between point 1 and point 2.
    """
    ...


def distance2(x1: float, y1: float, z1: float, x2: float, y2: float, z2: float) -> float:
    """
    Calculates the distance^2 between point 1 and point 2.
    """
    ...
