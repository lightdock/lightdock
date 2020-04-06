"""Tests for Residue class"""

from lightdock.structure.residue import Residue, AminoAcid, Cofactor, Ion
from lightdock.structure.atom import Atom


class TestResidue:

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_create_residue_empty_atoms(self):
        residue = Residue('ALA', 1)
        assert residue.name == 'ALA' and residue.number == 1 and len(residue.atoms) == 0

    def test_create_residue_with_atoms(self):
        residue = Residue('ALA', 1, [Atom(), Atom()])
        assert residue.name == 'ALA' and residue.number == 1 and len(residue.atoms) == 2
        
    def test_clone(self):
        residue1 = Residue('ALA', 1)
        residue2 = residue1.clone()
        
        assert residue1.name == residue2.name and residue2.number == residue2.number
        residue2.name = "MET"
        assert residue1.name != residue2.name and residue2.number == residue2.number
        
    def test_is_standard(self):
        residue1 = Residue('ALA', 1)
        residue2 = Residue('MG2', 2)
        
        assert residue1.is_standard()
        assert not residue2.is_standard()

    def test_is_nucleic(self):
        residue1 = Residue('DT', 1)
        residue2 = Residue('I', 2)
        
        assert residue1.is_nucleic()
        assert residue2.is_nucleic()
    
    def test_to_string(self):
        atoms = [Atom(1, 'CA', '', 'A', 'ALA'), Atom(2, 'N', '', 'A', 'ALA')] 
        residue = Residue('ALA', 1, atoms)
        assert "ALA.1    CA   0.000   0.000   0.000\nALA.1     N   0.000   0.000   0.000" == str(residue)

    def test_get_atom(self):
        atoms = [Atom(1, 'CA', '', 'A', 'ALA'), Atom(2, 'N', '', 'A', 'ALA')]
        residue = Residue('ALA', 1, atoms)

        atom1 = residue.get_atom('CA')
        atom2 = residue.get_atom('N')
        atom3 = residue.get_atom('CG')

        assert 'CA' == atom1.name
        assert 'N' == atom2.name
        assert None == atom3

    def test_backbone_sidechain(self):
        atoms = [Atom(1, 'CA', '', 'A', 'ALA'),
                 Atom(2, 'N', '', 'A', 'ALA'),
                 Atom(3, 'C', '', 'A', 'ALA'),
                 Atom(4, 'O', '', 'A', 'ALA'),
                 Atom(5, 'CB', '', 'A', 'ALA')]
        residue = Residue('ALA', 1, atoms)

        assert 1 == len(residue.sidechain)
        assert 'CB' == residue.sidechain[0].name
        assert 4 == len(residue.backbone)

    def test_dummy_residue(self):
        dummy = Residue.dummy(1.0, 2.0, 3.0)

        assert 1 == len(dummy.atoms)
        assert 'CA' == dummy.atoms[0].name
        assert "  CA   1.000   2.000   3.000" == str(dummy.atoms[0])


class TestAminoAcid:
    
    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_create_amino_acid(self):
        residue = AminoAcid('ALA', 1)
        assert residue.__class__.__name__ == "AminoAcid"


class TestCofactor:
    
    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_create_cofactor(self):
        residue = Cofactor('ATP',1)
        assert residue.__class__.__name__ == "Cofactor"
        

class TestIon:
    
    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_create_ion(self):
        residue = Ion('FE2',1)
        assert residue.__class__.__name__ == "Ion"
