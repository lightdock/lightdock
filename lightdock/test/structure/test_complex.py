"""Tests for Complex class"""

import filecmp
import shutil
import os
from pathlib import Path
import numpy as np
from nose.tools import assert_almost_equals
from lightdock.structure.complex import Complex
from lightdock.structure.chain import Chain
from lightdock.structure.residue import Residue
from lightdock.structure.atom import Atom
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.pdbutil.PDBIO import write_pdb_to_file


class TestComplex:
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
        self.chains = [Chain("A", self.residues)]

        self.atoms3 = [
            Atom(1, "N", "", "A", "CYS", 3, x=2.496, y=13.096, z=10.611),
            Atom(2, "CA", "", "A", "CYS", 3, x=2.787, y=12.161, z=9.557),
            Atom(3, "C", "", "A", "CYS", 3, x=4.052, y=11.431, z=9.896),
            Atom(4, "O", "", "A", "CYS", 3, x=5.120, y=12.044, z=9.912),
            Atom(5, "CB", "", "A", "CYS", 3, x=2.879, y=12.853, z=8.246),
            Atom(6, "SG", "", "A", "CYS", 3, x=3.492, y=11.911, z=6.838),
        ]

        self.atoms4 = [
            Atom(7, "N", "", "A", "PRO", 4, x=3.987, y=10.183, z=10.286),
            Atom(8, "CA", "", "A", "PRO", 4, x=5.219, y=9.484, z=10.695),
            Atom(9, "C", "", "A", "PRO", 4, x=6.277, y=9.424, z=9.662),
            Atom(10, "O", "", "A", "PRO", 4, x=5.993, y=9.385, z=8.431),
            Atom(11, "CB", "", "A", "PRO", 4, x=4.664, y=8.140, z=11.092),
            Atom(12, "CG", "", "A", "PRO", 4, x=3.295, y=8.087, z=10.812),
            Atom(13, "CD", "", "A", "PRO", 4, x=2.797, y=9.351, z=10.427),
        ]
        self.residues2 = [
            Residue("CYS", 1, "", self.atoms3),
            Residue("PRO", 2, "", self.atoms4),
        ]
        self.chains2 = [Chain("A", self.residues2)]

        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_complex"
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

    def test_create_empty_complex(self):
        protein = Complex([])

        assert len(protein.chains) == 0
        assert protein.num_atoms == 0
        assert protein.num_residues == 0

    def test_create_complex_only_with_chains(self):
        protein = Complex(chains=self.chains)

        expected_coordinates = np.array(
            [[1, 1, 1], [2, 2, 2], [1.1, 1.2, 1.3], [2.9, 2.8, 2.7]]
        )

        assert protein.num_atoms == 4
        assert protein.num_residues == 2
        assert (expected_coordinates == protein.atom_coordinates).all()
        for expected_index, atom in enumerate(protein.atoms):
            assert expected_index == atom.index

    def test_clone_complex(self):
        protein1 = Complex(chains=self.chains)
        protein2 = protein1.clone()

        assert protein1.num_atoms == protein2.num_atoms
        assert protein1.num_residues == protein2.num_residues
        for atom1, atom2 in zip(protein1.atoms, protein2.atoms):
            assert atom1 == atom2
        assert protein1.atom_coordinates == protein2.atom_coordinates

        for residue1, residue2 in zip(protein1.residues, protein2.residues):
            assert residue1 == residue2

        protein2.atoms[0].name = "C"
        for atom1, atom2 in zip(protein1.atoms[1:], protein2.atoms[1:]):
            assert atom1 == atom2
        assert protein1.atoms[0] != protein2.atoms[0]

    def test_center_of_mass_empty_complex(self):
        protein = Complex([])
        com = protein.center_of_mass()
        assert_almost_equals(0.0, com[0])
        assert_almost_equals(0.0, com[1])
        assert_almost_equals(0.0, com[2])

    def test_center_of_mass(self):
        protein = Complex(chains=self.chains)
        com = protein.center_of_mass()
        assert_almost_equals(1.8037228011238935, com[0])
        assert_almost_equals(1.7998854581864723, com[1])
        assert_almost_equals(1.7960481152490517, com[2])

    def test_center_of_coordinates(self):
        protein = Complex(chains=self.chains)
        cc = protein.center_of_coordinates()
        assert_almost_equals(1.75, cc[0])
        assert_almost_equals(1.75, cc[1])
        assert_almost_equals(1.75, cc[2])

    def test_center_of_coordinates_zero_atoms(self):
        protein = Complex(chains=[])
        cc = protein.center_of_coordinates()
        assert_almost_equals(0.0, cc[0])
        assert_almost_equals(0.0, cc[1])
        assert_almost_equals(0.0, cc[2])

    def test_translate(self):
        atom1 = Atom(2, "C", "", "A", "ALA", x=2.0, y=2.0, z=2.0)
        atom2 = Atom(2, "C", "", "A", "ALA", x=0.0, y=0.0, z=0.0)
        residue = Residue("ALA", 1, "", [atom1, atom2])
        chains = [Chain("A", [residue])]
        protein = Complex(chains)

        com = protein.center_of_mass()
        protein.translate([c * -1 for c in com])

        expected_coordinates = np.array([[1, 1, 1], [-1, -1, -1]])

        assert (expected_coordinates == protein.atom_coordinates).all()

    def test_move_to_origin(self):
        atom1 = Atom(2, "C", "", "A", "ALA", x=0.0, y=0.0, z=0.0)
        atom2 = Atom(2, "C", "", "A", "ALA", x=2.0, y=2.0, z=2.0)
        atom3 = Atom(2, "C", "", "A", "ALA", x=4.0, y=4.0, z=4.0)
        residue = Residue("ALA", 1, "", [atom1, atom2, atom3])
        chains = [Chain("A", [residue])]
        protein = Complex(chains)

        protein.move_to_origin()

        expected_coordinates = np.array(
            [[-2.0, -2.0, -2.0], [0, 0, 0], [2.0, 2.0, 2.0]]
        )

        assert (expected_coordinates == protein.atom_coordinates).all()

    def test_null_rotation(self):
        protein = Complex(chains=self.chains2)
        q = Quaternion()
        protein.rotate(q)
        write_pdb_to_file(
            protein, self.test_path / "rotated.pdb", protein.atom_coordinates[0]
        )
        assert filecmp.cmp(
            self.golden_data_path / "two_residues.pdb", self.test_path / "rotated.pdb"
        )

    def test_rotation_180_degrees_y_axis_origin_is_0(self):
        """Expected file has been generated with Chimera fixing the rotation to the
        center of coordinates and modifying the column of atom name to have the
        same padding as the write_pdb_file function.
        """
        protein = Complex(chains=self.chains2)
        q = Quaternion(0, 0.0, 1.0, 0)

        protein.rotate(q)

        write_pdb_to_file(
            protein, self.test_path / "rotated.pdb", protein.atom_coordinates[0]
        )
        assert filecmp.cmp(
            self.golden_data_path / "two_residues_y_180.pdb",
            self.test_path / "rotated.pdb",
        )

    def test_rotation_90_degrees_y_axis_90_degrees_x_axis_origin_is_0(self):
        """Expected file has been generated with Chimera fixing the rotation to the
        center of coordinates and modifying the column of atom name to have the
        same padding as the write_pdb_file function.
        """
        protein = Complex(chains=self.chains2)
        # Heading 90degrees (Y axis)
        q1 = Quaternion(0.7071, 0.0, 0.7071, 0.0)
        # Attitude 90degrees (X axis)
        q2 = Quaternion(0.7071, 0.0, 0.0, 0.7071)
        q = q1 * q2

        protein.rotate(q)

        write_pdb_to_file(
            protein, self.test_path / "rotated.pdb", protein.atom_coordinates[0]
        )
        assert filecmp.cmp(
            self.golden_data_path / "two_residues_y_90_x_90.pdb",
            self.test_path / "rotated.pdb",
        )
