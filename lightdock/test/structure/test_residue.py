"""Tests for Residue class"""

from lightdock.structure.residue import Residue, AminoAcid, Cofactor, Ion
from lightdock.structure.atom import Atom


class TestResidue:
    def test_create_residue_empty_atoms(self):
        residue = Residue("ALA", 1)
        assert residue.name == "ALA" and residue.number == 1 and len(residue.atoms) == 0

    def test_create_residue_with_atoms(self):
        residue = Residue("ALA", 1, "", [Atom(), Atom()])
        assert residue.name == "ALA" and residue.number == 1 and len(residue.atoms) == 2

    def test_clone(self):
        residue1 = Residue("ALA", 1)
        residue2 = residue1.clone()

        assert (
            residue1.name == residue2.name
            and residue1.number == residue2.number
            and residue1.insertion == residue2.insertion
        )
        residue2.name = "MET"
        assert (
            residue1.name != residue2.name
            and residue1.number == residue2.number
            and residue1.insertion == residue2.insertion
        )

    def test_is_standard(self):
        residue1 = Residue("ALA", 1)
        residue2 = Residue("MG2", 2)

        assert residue1.is_standard()
        assert not residue2.is_standard()

    def test_is_protein(self):
        residue1 = Residue("ALA", 1)
        residue2 = Residue("DT", 2)

        assert residue1.is_protein()
        assert not residue2.is_protein()

    def test_is_nucleic(self):
        residue1 = Residue("DT", 1)
        residue2 = Residue("I", 2)

        assert residue1.is_nucleic()
        assert residue2.is_nucleic()

    def test_to_string(self):
        atoms = [Atom(1, "CA", "", "A", "ALA"), Atom(2, "N", "", "A", "ALA")]
        residue = Residue("ALA", 1, "", atoms)
        assert (
            str(residue)
            == "ALA.1    CA   0.000   0.000   0.000\nALA.1     N   0.000   0.000   0.000"
        )

    def test_get_atom(self):
        atoms = [Atom(1, "CA", "", "A", "ALA"), Atom(2, "N", "", "A", "ALA")]
        residue = Residue("ALA", 1, "", atoms)

        atom1 = residue.get_atom("CA")
        atom2 = residue.get_atom("N")
        atom3 = residue.get_atom("CG")

        assert atom1.name == "CA"
        assert atom2.name == "N"
        assert atom3 is None

    def test_get_calpha(self):
        atoms = [Atom(1, "CA", "", "A", "ALA"), Atom(2, "N", "", "A", "ALA")]
        residue = Residue("ALA", 1, "", atoms)

        atom = residue.get_calpha()

        assert atom.name == "CA"

        atoms = [Atom(1, "CB", "", "A", "ALA"), Atom(2, "N", "", "A", "ALA")]
        residue = Residue("ALA", 1, "", atoms)

        atom = residue.get_calpha()

        assert atom is None

    def test_get_non_hydrogen_atoms(self):
        atoms = [
            Atom(1, "CA", "", "A", "ALA"),
            Atom(2, "N", "", "A", "ALA"),
            Atom(3, "H", "", "A", "ALA"),
        ]
        residue = Residue("ALA", 1, "", atoms)

        no_h = residue.get_non_hydrogen_atoms()

        assert len(no_h) == 2
        assert no_h[0].name == "CA"
        assert no_h[1].name == "N"

    def test_backbone_sidechain(self):
        atoms = [
            Atom(1, "CA", "", "A", "ALA"),
            Atom(2, "N", "", "A", "ALA"),
            Atom(3, "C", "", "A", "ALA"),
            Atom(4, "O", "", "A", "ALA"),
            Atom(5, "CB", "", "A", "ALA"),
        ]
        residue = Residue("ALA", 1, "", atoms)

        assert len(residue.sidechain) == 1
        assert residue.sidechain[0].name == "CB"
        assert len(residue.backbone) == 4

    def test_dummy_residue(self):
        dummy = Residue.dummy(1.0, 2.0, 3.0)

        assert len(dummy.atoms) == 1
        assert dummy.atoms[0].name == "CA"
        assert str(dummy.atoms[0]) == "  CA   1.000   2.000   3.000"


class TestAminoAcid:
    def test_create_amino_acid(self):
        residue = AminoAcid("ALA", 1)
        assert residue.__class__.__name__ == "AminoAcid"


class TestCofactor:
    def test_create_cofactor(self):
        residue = Cofactor("ATP", 1)
        assert residue.__class__.__name__ == "Cofactor"


class TestIon:
    def test_create_ion(self):
        residue = Ion("FE2", 1)
        assert residue.__class__.__name__ == "Ion"
