from lightdock.constants import DEFAULT_ROTATION_STEP
from lightdock.mathutil.lrandom import RandomNumberGenerator
from typing import Optional


class Quaternion:
    """An abstract data type representing a quaternion."""

    def __init__(self, w: float = 1., x: float = 0., y: float = 0., z: float = 0.) -> None:
        """
        Builds a quaternion.

        If not parameters are defined, returns the identity quaternion
        """
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def clone(self) -> Quaternion:
        """
        Creates a new instance of this quaternion
        """
        ...

    def __eq__(self, other: Quaternion) -> bool:
        """
        Compares two quaternions for equality using their components
        """
        ...

    def __ne__(self, other: Quaternion) -> bool:
        """
        Negation of the __eq__ function
        """
        ...

    def __neg__(self) -> Quaternion:
        """
        Implements quaternion inverse
        """
        ...

    def __add__(self, other: Quaternion) -> Quaternion:
        """
        Implements quaternion addition
        """
        ...

    def __sub__(self, other: Quaternion) -> Quaternion:
        """
        Implements quaternion substract
        """
        ...

    def __rmul__(self, scalar: float) -> Quaternion:
        """
        Implements multiplication of the form scalar*quaternion
        """
        ...

    def conjugate(self) -> Quaternion:
        """
        Calculates the conjugate of this quaternion
        """
        ...

    def __mul__(self, other: Quaternion) -> Quaternion:
        """
        Calculates quaternion multiplication
        """
        ...

    def __truediv__(self, scalar: float) -> Quaternion:
        """
        Calculates division of quaternion by scalar
        """
        ...

    def dot(self, other: Quaternion) -> float:
        """
        Calculates the dot product of two quaternions
        """
        ...

    def norm(self) -> float:
        """
        Calculates quaternion norm
        """
        ...

    def norm2(self) -> float:
        """
        Calculates quaternion norm^2
        """
        ...

    def normalize(self) -> Quaternion:
        """
        Normalizes a given quaternion
        """
        ...

    def inverse(self) -> Quaternion:
        """
        Calculates the inverse of this quaternion
        """
        ...

    def rotate(self, vec3) -> list[float]:
        """
        Rotates vec3 using quaternion
        """
        ...

    def __repr__(self) -> str:
        """
        Vector representation of the quaternion
        """
        ...

    def lerp(self, other: Quaternion, t: float) -> Quaternion:
        """
        Calculates the linear interpolation between two quaternions
        """
        ...

    def slerp(self, other: Quaternion, t: float = DEFAULT_ROTATION_STEP) -> Quaternion:
        """
        Calculates the spherical linear interpolation of two quaternions given a t step
        """
        ...

    def distance(self, other: Quaternion) -> float:
        """
        Calculates the closeness of two orientations represented in quaternions space.

        Quaternions must be normalized. Distance is 0 when quaternions are equal and 1 when
        the orientations are 180 degrees apart.
        See http://math.stackexchange.com/questions/90081/quaternion-distance
        """
        ...

    @staticmethod
    def random(rng: Optional[RandomNumberGenerator] = None) -> Quaternion:
        """
        Generates a random quaternion uniformly distributed:
        http://planning.cs.uiuc.edu/node198.html
        """
        ...
