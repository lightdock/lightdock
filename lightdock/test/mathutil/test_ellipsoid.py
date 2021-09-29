"""Tests for Ellipsoid module"""

import os
import shutil
from pathlib import Path
import numpy as np
from nose.tools import assert_almost_equal, raises
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.mathutil.ellipsoid import MinimumVolumeEllipsoid
from lightdock.error.lightdock_errors import MinimumVolumeEllipsoidError
from lightdock.mathutil.constants import ERROR_TOLERANCE


class TestEllipsoid:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_ellipsoid"
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

    def test_calculate_min_volume_ellipsoid(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPE_l_u.pdb"
        )
        protein = Complex(chains, atoms)

        ellipsoid = MinimumVolumeEllipsoid(protein.atom_coordinates[0].coordinates)

        assert_almost_equal(5.79979144, ellipsoid.center[0])
        assert_almost_equal(13.30609275, ellipsoid.center[1])
        assert_almost_equal(6.28378695, ellipsoid.center[2])

        assert_almost_equal(11.51000999, ellipsoid.radii[0])
        assert_almost_equal(17.41300089, ellipsoid.radii[1])
        assert_almost_equal(25.1317681, ellipsoid.radii[2])

        assert_almost_equal(-0.64868458, ellipsoid.rotation[0][0])
        assert_almost_equal(-0.43420895, ellipsoid.rotation[0][1])
        assert_almost_equal(0.62503673, ellipsoid.rotation[0][2])
        assert_almost_equal(0.75208829, ellipsoid.rotation[1][0])
        assert_almost_equal(-0.49144928, ellipsoid.rotation[1][1])
        assert_almost_equal(0.43913643, ellipsoid.rotation[1][2])
        assert_almost_equal(0.11649688, ellipsoid.rotation[2][0])
        assert_almost_equal(0.75494384, ellipsoid.rotation[2][1])
        assert_almost_equal(0.64535903, ellipsoid.rotation[2][2])

        expected_poles = [
            [13.266157381855532, 18.303842059830465, -0.91039204503235993],
            [-1.6665744987943629, 8.3083434397316651, 13.477965942387469],
            [-7.296322648355452, 21.863699498711949, -1.3628961564119457],
            [18.89590553141662, 4.7484860008501819, 13.930470053767054],
            [2.8720188105117521, -5.6669806736857815, -9.9352265853089641],
            [8.7275640725494164, 32.279166173247916, 22.502800482664075],
        ]

        assert np.allclose(expected_poles, ellipsoid.poles, ERROR_TOLERANCE)

    @raises(MinimumVolumeEllipsoidError)
    def test_exception_singular_matrix(self):
        coordinates = np.array([[2.0, 2.0, 2.0], [0.0, 0.0, 0.0]])

        ellipsoid = MinimumVolumeEllipsoid(coordinates)

        assert len(ellipsoid.poles) > 0
