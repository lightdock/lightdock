"""Tests for poses related module"""

import os
import shutil
import filecmp
import numpy as np
from nose.tools import assert_almost_equal
from lightdock.prep.poses import normalize_vector, quaternion_from_vectors
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.structure.complex import Complex
from lightdock.pdbutil.PDBIO import parse_complex_from_file, create_pdb_from_points


class TestPoses:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_poses/'
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'

    def tearDown(self):
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
    
    def test_normalize_vector1(self):
        v = np.array([1., 0., 0.])
        e = np.array([1., 0., 0.])
        n = normalize_vector(v)
        assert np.allclose(n, e)

    def test_normalize_vector2(self):
        v = np.array([4., 0., 0.])
        e = np.array([1., 0., 0.])
        n = normalize_vector(v)
        assert np.allclose(n, e)

    def test_normalize_vector3(self):
        v = np.array([0., 0., 0.])
        e = np.array([0., 0., 0.])
        n = normalize_vector(v)
        assert np.allclose(n, e)

    def test_normalize_vector4(self):
        v = np.array([2., -2., 0.])
        e = np.array([0.70710678, -0.70710678, 0.])
        n = normalize_vector(v)
        assert np.allclose(n, e)

    def test_quaternion_from_vectors1(self):
        a = np.array([1., 0., 0.])
        b = np.array([1., 0., 0.])
        q = quaternion_from_vectors(a, b)
        e = Quaternion()
        assert e == q

    def test_quaternion_from_vectors2(self):
        a = np.array([1., 0., 0.])
        b = np.array([2., 0., 0.])
        q = quaternion_from_vectors(a, b)
        e = Quaternion()
        assert e == q

    def test_quaternion_from_vectors3(self):
        a = np.array([1., 0., 0.])
        b = np.array([0., 2., 0.])
        q = quaternion_from_vectors(a, b)
        # 90 degrees in Z axis
        e = Quaternion(w=0.70710678, x=0.00000000, y=0.00000000, z=0.70710678)
        assert e == q

    def test_quaternion_from_vectors4(self):
        a = np.array([1., 0., 0.])
        b = np.array([-1., 0., 0.])
        q = quaternion_from_vectors(a, b)
        # Orthogonal rotation
        e = Quaternion(w=0., x=0.00000000, y=-1.0, z=0.)
        assert e == q
