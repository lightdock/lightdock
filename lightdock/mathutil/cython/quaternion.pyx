#cython: boundscheck=False
#cython: wraparound=False
#cython: language_level=3

from lightdock.mathutil.cython.cutil import float_equals as cfloat_equals
from math import sqrt, acos, sin, cos, pi
import random
from lightdock.constants import DEFAULT_ROTATION_STEP
from lightdock.mathutil.constants import LINEAR_THRESHOLD

class Quaternion:

    def __init__(self, w=1., x=0., y=0., z=0.):
        """
        Builds a quaternion.

        If not parameters are defined, returns the identity quaternion
        """
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def clone(self):
        """
        Creates a new instance of this quaternion
        """
        return Quaternion(self.w, self.x, self.y, self.z)

    def __eq__(self, other):
        """
        Compares two quaternions for equality using their components
        """
        return cfloat_equals(self.w, other.w) and \
            cfloat_equals(self.x, other.x) and \
            cfloat_equals(self.y, other.y) and \
            cfloat_equals(self.z, other.z)

    def __ne__(self, other):
        """
        Negation of the __eq__ function
        """
        return not self.__eq__(other)

    def __neg__(self):
        """
        Implements quaternion inverse
        """
        return Quaternion(-self.w, -self.x, -self.y, -self.z)

    def __add__(self, other):
        """
        Implements quaternion addition
        """
        return Quaternion(self.w+other.w, self.x+other.x, self.y+other.y, self.z+other.z)

    def __sub__(self, other):
        """
        Implements quaternion substract
        """
        return Quaternion(self.w-other.w, self.x-other.x, self.y-other.y, self.z-other.z)

    def __rmul__(self, scalar):
        """
        Implements multiplication of the form scalar*quaternion
        """
        return Quaternion(scalar * self.w, scalar * self.x, scalar * self.y, scalar * self.z)

    def conjugate(self):
        """
        Calculates the conjugate of this quaternion
        """
        return Quaternion(self.w, -self.x, -self.y, -self.z)

    def __mul__(self, other):
        """
        Calculates quaternion multiplication
        """
        w = (self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z)
        x = (self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y)
        y = (self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x)
        z = (self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w)
        return Quaternion(w, x, y, z)

    def __truediv__(self, scalar):
        """
        Calculates division of quaternion by scalar
        """
        return Quaternion(self.w/scalar, self.x/scalar, self.y/scalar, self.z/scalar)

    def dot(self, other):
        """
        Calculates the dot product of two quaternions
        """
        return self.w*other.w + self.x*other.x + self.y*other.y + self.z*other.z

    def norm(self):
        """
        Calculates quaternion norm
        """
        return sqrt(self.norm2())

    def norm2(self):
        """
        Calculates quaternion norm^2
        """
        return (self.w*self.w) + (self.x*self.x) + (self.y*self.y) + (self.z*self.z)

    def normalize(self):
        """
        Normalizes a given quaternion
        """
        return self/self.norm()

    def inverse(self):
        """
        Calculates the inverse of this quaternion
        """
        return self.conjugate() / float(self.norm2())

    def rotate(self, vec3):
        """
        Rotates vec3 using quaternion
        """
        v = Quaternion(0., vec3[0], vec3[1], vec3[2])
        r = self*v*self.inverse()
        return [r.x, r.y, r.z]


    def __repr__(self):
        """
        Vector representation of the quaternion
        """
        return "(%10.8f, %10.8f, %10.8f, %10.8f)" % (self.w, self.x, self.y, self.z)

    def lerp(self, other, t):
        """
        Calculates the linear interpolation between two quaternions
        """
        return (1.0-t)*self + t*other

    def slerp(self, other, t=DEFAULT_ROTATION_STEP):
        """
        Calculates the spherical linear interpolation of two quaternions given a t step
        """
        self = self.normalize()
        other = other.normalize()
        q_dot = self.dot(other)
        # Patch to avoid the long path
        if q_dot < 0:
            self = -self
            q_dot *= -1.

        if q_dot > LINEAR_THRESHOLD:
            # Linear interpolation if quaternions are too close
            result = self + t*(other-self)
            result.normalize()
            return result
        else:
            q_dot = max(min(q_dot, 1.0), -1.0)
            omega = acos(q_dot)
            so = sin(omega)
            return (sin((1.0-t)*omega) / so) * self + (sin(t*omega)/so) * other


    def distance(self, other):
        """
        Calculates the closeness of two orientations represented in quaternions space.

        Quaternions must be normalized. Distance is 0 when quaternions are equal and 1 when
        the orientations are 180 degrees apart.
        See http://math.stackexchange.com/questions/90081/quaternion-distance
        """
        return 1-self.dot(other)**2

    @staticmethod
    def random(rng=None):
        """
        Generates a random quaternion uniformly distributed:
        http://planning.cs.uiuc.edu/node198.html
        """
        u1 = 0.
        u2 = 0.
        u3 = 0.
        if rng:
            u1 = rng()
            u2 = rng()
            u3 = rng()
        else:
            u1 = random.random()
            u2 = random.random()
            u3 = random.random()
        return Quaternion(sqrt(1-u1)*sin(2*pi*u2), sqrt(1-u1)*cos(2*pi*u2),
                            sqrt(u1)*sin(2*pi*u3), sqrt(u1)*cos(2*pi*u3))
