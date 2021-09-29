"""Tests for Quaternion class"""

from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.mathutil.lrandom import MTGenerator
from nose.tools import assert_almost_equals


class TestQuaternion:
    def test_create_identity_quaternion(self):
        q = Quaternion()
        assert_almost_equals(1.0, q.w)
        assert_almost_equals(0.0, q.x)
        assert_almost_equals(0.0, q.y)
        assert_almost_equals(0.0, q.z)

    def test_create_quaternion(self):
        q = Quaternion(1.0, 1.2345, 6.7890, 2.3456)
        assert_almost_equals(1.0, q.w)
        assert_almost_equals(1.2345, q.x)
        assert_almost_equals(6.7890, q.y)
        assert_almost_equals(2.3456, q.z)

    def test_equal_quaternions(self):
        q1 = Quaternion(1.0, 2.0, -1.0, 3.0)
        q2 = Quaternion(1.0, 2.0, -1.0, 3.0)

        assert q1 == q2

    def test_copy_quaternions(self):
        q1 = Quaternion(1.0, 2.0, -1.0, 3.0)
        q2 = q1.clone()
        assert q1 == q2
        q2.w = 3.0
        assert q1 != q2

    def test_add_quaternions(self):
        q1 = Quaternion(1.0, 1.0, -1.0, 2.0)
        q2 = Quaternion(1.0, 1.0, 3.0, 0.0)

        expected = Quaternion(2.0, 2.0, 2.0, 2.0)

        assert expected == (q1 + q2)

    def test_subtract_quaternions(self):
        q1 = Quaternion(3.0, 1.0, -1.0, 0.0)
        q2 = Quaternion(2.0, 1.0, -1.0, 0.0)

        expected = Quaternion()

        assert expected == (q1 - q2)

    def test_negate_quaternion(self):
        q1 = Quaternion(3.0, 1.0, -1.0, 0.0)

        expected = Quaternion(-3.0, -1.0, 1.0, 0.0)

        assert expected == -q1

    def test_rmul(self):
        q1 = Quaternion(1, 1, -1, 2)

        expected = Quaternion(3, 3, -3, 6)

        assert expected == 3 * q1

    def test_division_by_scalar(self):
        q1 = Quaternion(1, 1, -1, 2)

        expected = Quaternion(1 / 2.0, 1 / 2.0, -1 / 2.0, 1)

        assert expected == q1 / 2.0

    def test_conjugate(self):
        q1 = Quaternion(1, 3, 4, 3)
        q2 = Quaternion(1, -1, -1, 3)

        expected1 = Quaternion(1, -3, -4, -3)
        expected2 = Quaternion(1, 1, 1, -3)

        assert expected1 == q1.conjugate()
        assert expected2 == q2.conjugate()

    def test_multiply_quaternions(self):
        q1 = Quaternion(1, 0, 0, 2)
        q2 = Quaternion(3, -1, 4, 3)
        q3 = Quaternion(1 / 2.0, -3, 2, 9)

        expected1 = Quaternion(-3, -9, 2, 9)
        expected2 = Quaternion(-3, 7, 6, 9)
        expected3 = Quaternion(-147.0 / 2, 97.0 / 2, -93.0, 19.0 / 2)

        assert expected1 == q1 * q2
        assert expected2 == q2 * q1
        assert expected3 == q2 * q1 * q3

    def test_conjugate_and_multiplication(self):
        q1 = Quaternion(1, 0, 0, 2)
        q2 = Quaternion(3, -1, 4, 3)

        expected = Quaternion(35, 0, 0, 0)

        assert (q1 * q2).conjugate() == q2.conjugate() * q1.conjugate()
        assert expected == q2.conjugate() * q2

    def test_norm(self):
        q1 = Quaternion(1, -3, 4, 3)
        q2 = Quaternion(3, -1, 4, 3)

        assert_almost_equals(5.91607978, q1.norm())
        assert_almost_equals((q1 * q2).norm(), q1.norm() * q2.norm())

    def test_inverse(self):
        q1 = Quaternion(1, 0, 0, 2)
        q2 = Quaternion(3, -1, 4, 3)

        expected = Quaternion(-3.0 / 175, 9.0 / 175, -2.0 / 175, -9.0 / 175)

        assert expected == (q1 * q2).inverse()

    def test_rotation(self):
        q = Quaternion(0.707106781, 0.0, 0.707106781, 0.0)
        v = [1.0, 0.0, 0.0]

        v2 = q.rotate(v)

        assert_almost_equals(0.0, v2[0])
        assert_almost_equals(0.0, v2[1])
        assert_almost_equals(-1.0, v2[2])

    def test_dot_product(self):
        q = Quaternion(0.707106781, 0.0, 0.707106781, 0.0)

        assert_almost_equals(1.0, q.dot(q))

    def test_lerp_t_0(self):
        q1 = Quaternion(1, 0, 0, 2)
        q2 = Quaternion(3, -1, 4, 3)

        s = q1.lerp(q2, 0.0)

        assert s == q1

    def test_lerp_t_1(self):
        q1 = Quaternion(1, 0, 0, 2)
        q2 = Quaternion(3, -1, 4, 3)

        s = q1.lerp(q2, 1.0)

        assert s == q2

    def test_slerp_t_0(self):
        q1 = Quaternion(1, 0, 0, 2)
        q2 = Quaternion(3, -1, 4, 3)
        expected = Quaternion(0.44721360, 0.000000, 0.000000, 0.89442719)

        s = q1.slerp(q2, 0.0)

        assert expected == s

    def test_slerp_t_1(self):
        q1 = Quaternion(1, 0, 0, 2)
        q2 = Quaternion(3, -1, 4, 3)
        expected = Quaternion(0.50709255, -0.16903085, 0.67612340, 0.50709255)

        s = q1.slerp(q2, 1.0)

        assert expected == s

    def test_slerp_t_half_y(self):
        q1 = Quaternion(1, 0, 0, 0)
        q2 = Quaternion(0, 0, 1, 0)

        s = q1.slerp(q2, 0.5)

        expected = Quaternion(0.707106781, 0.0, 0.707106781, 0.0)

        assert expected == s

    def test_slerp_t_half_bank_zero(self):
        q1 = Quaternion(0.707106781, 0.0, 0.0, 0.707106781)
        q2 = Quaternion(0.0, 0.707106781, 0.707106781, 0.0)

        s = q1.slerp(q2, 0.5)

        expected = Quaternion(0.5, 0.5, 0.5, 0.5)

        assert expected == s

    def test_slerp_same_quaternion(self):
        q1 = Quaternion(0.707106781, 0.0, 0.0, 0.707106781)
        q2 = Quaternion(0.707106781, 0.0, 0.0, 0.707106781)

        s = q1.slerp(q2, 0.1)

        assert s == q2

    def test_distance_is_zero(self):
        q = Quaternion(0.707106781, 0.0, 0.707106781, 0.0)

        assert_almost_equals(0.0, q.distance(q))

    def test_distance_is_one(self):
        q1 = Quaternion(0.707106781, 0.0, 0.707106781, 0.0)
        q2 = Quaternion(0.707106781, 0.0, -0.707106781, 0.0)

        assert_almost_equals(1.0, q1.distance(q2))

    def test_distance_is_half(self):
        q1 = Quaternion(0.707106781, 0.0, 0.707106781, 0.0)
        q2 = Quaternion(0, 0, 1.0, 0)

        assert_almost_equals(0.5, q1.distance(q2))

    def test_distance_composite_rotation(self):
        q1 = Quaternion(1.0, 0, 0, 0)
        q2 = Quaternion(0.5, 0.5, 0.5, 0.5)

        assert_almost_equals(0.75, q1.distance(q2))

    def test_get_random(self):
        q1 = Quaternion.random()
        assert_almost_equals(1.0, q1.norm())

        rng = MTGenerator(1379)
        q2 = Quaternion.random(rng)
        expected = Quaternion(0.53224081, -0.00589472, 0.84280352, 0.07979471)

        assert expected == q2
