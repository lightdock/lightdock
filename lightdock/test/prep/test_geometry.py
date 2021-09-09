"""Tests for geometry module"""

import os
from pathlib import Path
import shutil
import filecmp
from lightdock.prep.geometry import sphere, axis, create_bild_file


class TestGeometry:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_geometry"
        self.golden_data_path = self.path / "golden_data"

    def setup(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.mkdir(self.test_path)

    def teardown(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass

    def test_sphere(self):
        center = [5.0, 5.0, 5.0]
        radius = 0.3

        expected = ".sphere 5.000000 5.000000 5.000000 0.300000"

        assert expected == sphere(center, radius)

    def test_axis_origin_length_default(self):
        # pose is translation + quaternion
        pose = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]

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
        poses = [
            [0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5],
            [2.0, 2.0, 2.0, 0.0, 0.7071, 0.7071, 0.0],
            [-2.0, -2.0, -2.0, 0.0, 0.0, 0.7071, 0.7071],
        ]

        generated_file = self.test_path / "test.bild"
        expected_file = self.golden_data_path / "test.bild"

        create_bild_file(generated_file, poses)

        assert filecmp.cmp(generated_file, expected_file)
