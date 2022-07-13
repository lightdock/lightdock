"""Tests for starting_points module"""

import os
import shutil
import filecmp
from pathlib import Path
from nose.tools import assert_almost_equal
from lightdock.prep.starting_points import points_on_sphere, calculate_surface_points
from lightdock.structure.complex import Complex
from lightdock.pdbutil.PDBIO import parse_complex_from_file, create_pdb_from_points


class TestStartingPoints:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_starting_points"
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

    def test_points_on_sphere(self):
        points = points_on_sphere(5)
        expected = [
            [0.5999999999999999, -0.8, 0.0],
            [-0.6758097397797129, -0.39999999999999997, 0.6190970809322855],
            [0.08742572471695988, 5.551115123125783e-17, -0.9961710408648278],
            [0.5576434272376701, 0.4000000000000002, 0.7273471028736039],
            [-0.590828091189257, 0.8, -0.10450917022758696],
        ]
        for correct, point in zip(expected, points):
            assert_almost_equal(correct[0], point[0])
            assert_almost_equal(correct[1], point[1])
            assert_almost_equal(correct[2], point[2])

    def test_create_pdb_from_points(self):
        points = points_on_sphere(100)
        create_pdb_from_points(self.test_path / "points.pdb", points)
        assert filecmp.cmp(
            self.golden_data_path / "100_points.pdb", self.test_path / "points.pdb"
        )

    def test_calculate_starting_points(self):
        # Receptor
        file_name = self.golden_data_path / "1PPE_rec.pdb"
        atoms, _, chains = parse_complex_from_file(file_name)
        receptor = Complex(chains, atoms, structure_file_name=file_name)

        # Ligand
        file_name = self.golden_data_path / "1PPE_lig.pdb"
        atoms, _, chains = parse_complex_from_file(file_name)
        ligand = Complex(chains, atoms, structure_file_name=file_name)

        starting_points, rec_diameter, lig_diameter = calculate_surface_points(
            receptor, ligand, 100, [0.0, 0.0, 0.0], 50.0, [], [],
        )

        assert_almost_equal(rec_diameter, 50.213210831413676)
        assert_almost_equal(lig_diameter, 27.855559534857672)

        create_pdb_from_points(self.test_path / "points.pdb", starting_points)
        assert filecmp.cmp(
            self.golden_data_path / "starting_points.pdb", self.test_path / "points.pdb"
        )
