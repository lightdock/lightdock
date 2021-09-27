"""Tests for Model class"""

import numpy as np
from lightdock.structure.model import DockingModel
from lightdock.structure.residue import Residue
from lightdock.structure.atom import Atom
from lightdock.structure.space import SpacePoints
from lightdock.mathutil.cython.quaternion import Quaternion


class TestDockingModel:
    def __init__(self):
        self.atoms1 = [
            Atom(1, "CA", "", "A", "ALA", x=1.0, y=1.0, z=1.0),
            Atom(2, "N", "", "A", "ALA", x=2.0, y=2.0, z=2.0),
        ]
        self.atoms2 = [
            Atom(3, "CA", "", "A", "HIS", x=1.1, y=1.2, z=1.3),
            Atom(4, "N", "", "A", "HIS", x=2.9, y=2.8, z=2.7),
        ]
        self.residues = [
            Residue("ALA", 1, "", self.atoms1),
            Residue("HIS", 2, "", self.atoms2),
        ]

    def test_create_model(self):
        docking_model = DockingModel(
            objects=self.residues, coordinates=SpacePoints([[1, 1, 1], [1.1, 1.2, 1.3]])
        )

        expected_coordinates = SpacePoints([[1, 1, 1], [1.1, 1.2, 1.3]])

        assert len(docking_model.objects) == 2
        assert expected_coordinates == docking_model.coordinates[0]

    def test_create_full_model(self):
        restraints = {"active": [], "passive": []}
        membrane = []
        reference_points = [[0, 0, 0]]
        n_modes = [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
        nm_mask = [True, True, True, True]
        docking_model = DockingModel(
            objects=self.residues,
            coordinates=SpacePoints([[1, 1, 1], [1.1, 1.2, 1.3]]),
            restraints=restraints,
            membrane=membrane,
            reference_points=reference_points,
            n_modes=n_modes,
            nm_mask=nm_mask,
        )

        assert "active" in docking_model.restraints
        assert len(docking_model.n_modes) == 4
        assert len(docking_model.nm_mask) == 4 and all(docking_model.nm_mask)

    def test_translate(self):
        atom1 = Atom(1, "CA", "", "A", "ALA", x=2.0, y=2.0, z=2.0)
        atom2 = Atom(2, "CB", "", "A", "ALA", x=0.0, y=0.0, z=0.0)
        residues = [Residue("ALA", 1, "", [atom1, atom2])]
        docking_model = DockingModel(
            residues, SpacePoints([[2.0, 2.0, 2.0], [0, 0, 0]])
        )

        docking_model.translate([-2, -2, -2])

        expected_coordinates = SpacePoints([[0, 0, 0], [-2, -2, -2]])
        assert expected_coordinates == docking_model.coordinates[0]

        expected_coordinates = SpacePoints([[-1, -1, -1]])
        assert np.allclose(expected_coordinates, docking_model.reference_points)

    def test_null_rotation_one_atom(self):
        atom1 = Atom(1, "CA", "", "A", "ALA", x=2.0, y=2.0, z=2.0)
        atoms = [atom1]
        docking_model = DockingModel(atoms, SpacePoints([[2.0, 2.0, 2.0]]))
        q = Quaternion()
        docking_model.rotate(q)

        expected_coordinates = np.array([[2.0, 2.0, 2.0]])
        assert (expected_coordinates == docking_model.coordinates).all()

    def test_null_rotation(self):
        atom1 = Atom(1, "CA", "", "A", "ALA", x=2.0, y=2.0, z=2.0)
        atom2 = Atom(2, "CB", "", "A", "ALA", x=0.0, y=0.0, z=0.0)
        atom3 = Atom(3, "C", "", "A", "ALA", x=0.5, y=0.5, z=0.5)
        atoms = [atom1, atom2, atom3]
        docking_model = DockingModel(
            atoms, SpacePoints([[2.0, 2.0, 2.0], [0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
        )
        q = Quaternion()
        docking_model.rotate(q)

        expected_coordinates = SpacePoints(
            [[2.0, 2.0, 2.0], [0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
        )
        assert expected_coordinates == docking_model.coordinates[0]

        expected_coordinates = SpacePoints(
            [[0.833333333333, 0.833333333333, 0.833333333333]]
        )
        assert np.allclose(expected_coordinates, docking_model.reference_points)

    def test_reference_points_singular_matrix(self):
        atom1 = Atom(1, "CA", "", "A", "ALA", x=2.0, y=2.0, z=2.0)
        atom2 = Atom(2, "CB", "", "A", "ALA", x=0.0, y=0.0, z=0.0)
        atom3 = Atom(3, "C", "", "A", "ALA", x=1.0, y=1.0, z=1.0)
        atoms = [atom1, atom2, atom3]
        docking_model = DockingModel(
            atoms, SpacePoints([[2.0, 2.0, 2.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        )

        expected_coordinates = SpacePoints([[1.0, 1.0, 1.0]])

        assert expected_coordinates == docking_model.reference_points

    def test_reference_points_minimum_volume_ellipsoid(self):
        atom1 = Atom(1, "CA", "", "A", "ALA", x=1.2, y=-1.0, z=2.0)
        atom2 = Atom(2, "CB", "", "A", "ALA", x=0.0, y=0.0, z=0.0)
        atom3 = Atom(3, "C", "", "A", "ALA", x=0.5, y=3.0, z=0.5)
        atom4 = Atom(4, "N", "", "A", "ALA", x=0.85, y=-2.5, z=0.4)
        atoms = [atom1, atom2, atom3, atom4]
        docking_model = DockingModel(
            atoms,
            SpacePoints(
                [[1.2, -1.0, 2.0], [0.0, 0.0, 0.0], [0.5, 3.0, 0.5], [0.85, -2.5, 0.4]]
            ),
        )
        # Only center is used now
        expected_coordinates = SpacePoints([[0.6375, -0.125, 0.725]])
        assert expected_coordinates == docking_model.reference_points
