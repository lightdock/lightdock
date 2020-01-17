"""Tests for geometry module"""

import os
import shutil
import filecmp
from nose.tools import assert_almost_equal
from lightdock.prep.geometry import sphere, axis, create_bild_file


class TestGeometry:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch/'
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
    
    def test_sphere(self):
        center = [5., 5., 5.]
        radius = 0.3

        expected = ".sphere 5.000000 5.000000 5.000000 0.300000"

        assert expected == sphere(center, radius)

    def test_axis_origin_length_default(self):
        # pose is translation + quaternion
        pose = [0., 0., 0., 1., 0., 0., 0.]

        expected = """.color cornflower blue
.sphere 0.000000 0.000000 0.000000 0.300000
.color 1 0 0
.arrow 0.000000 0.000000 0.000000 2.000000 0.000000 0.000000
.color 1 1 0
.arrow 0.000000 0.000000 0.000000 0.000000 2.000000 0.000000
.color 0 0 1
.arrow 0.000000 0.000000 0.000000 0.000000 0.000000 2.000000
"""

        assert expected == axis(pose)
        
    def test_create_bild_file(self):
        poses = [[0., 0., 0., 0.5, 0.5, 0.5, 0.5],
                 [2., 2., 2., 0., 0.7071, 0.7071, 0.],
                 [-2., -2., -2., 0., 0., 0.7071, 0.7071]]
        
        generated_file = os.path.join(self.test_path, 'test.bild')
        expected_file = os.path.join(self.golden_data_path, 'test.bild')

        create_bild_file(generated_file, poses)
    
        assert filecmp.cmp(generated_file, expected_file)
