"""Tests for Atom class"""

from lightdock.structure.atom import Atom, HetAtom
from lightdock.error.lightdock_errors import AtomError
from nose.tools import assert_almost_equals #@UnresolvedImport
from nose.tools import raises


class TestAtom:

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_create_empty_atom(self):
        atom = Atom()
        assert_almost_equals(0.0, atom.x)
    
    def test_assign_element_and_mass(self):
        atom = Atom(1, 'CA', '', 'A', 'ALA', 1, '', x=0., y=0., z=0., 
                    occupancy=1., b_factor=0., element=None, mass=None)
        assert 'C' == atom.element
        assert_almost_equals(Atom.MASSES[atom.element], atom.mass)

    def test_create_atom_not_given_element(self):
        atom = Atom(1, 'PB', '', 'A', 'ALA', 1, '', x=0., y=0., z=0., 
                    occupancy=1., b_factor=0., element=None, mass=None)
        assert 'PB' == atom.element
        assert_almost_equals(Atom.MASSES[atom.element], atom.mass)
  
    @raises(AtomError)
    def test_assign_element_and_mass_of_not_recognized_element(self):
        atom = Atom(1, 'Ty', '', 'A', 'BSG', 1, '', x=0., y=0., z=0., 
                    occupancy=1., b_factor=0., element=None, mass=None)
        assert atom.name != 'Ty'
    
    @raises(AtomError)
    def test_not_recognized_element(self):
        """Test if Tylium is recognized"""
        atom = Atom(1, 'Ty', '', 'A', 'BSG', 1, '', x=0., y=0., z=0., 
                    occupancy=1., b_factor=0., element='Ty', mass=None)
        assert atom.name != 'Ty'
        
    def test_is_hydrogen(self):
        atom1 = Atom(1, 'CA', '', 'A', 'ALA', 1, '', x=1., y=2., z=-3., 
                    occupancy=1., b_factor=0., element=None, mass=None)
        atom2 = Atom()
        
        assert not atom1.is_hydrogen() and atom2.is_hydrogen()
    
    def test_is_backbone(self):
        atom1 = Atom(1, 'CA', '', 'A', 'ALA', 1, '', x=1., y=2., z=-3., 
                    occupancy=1., b_factor=0., element=None, mass=None)
        atom2 = Atom()
        
        assert atom1.is_backbone() and not atom2.is_backbone()
    
    def test_distance(self):
        atom1 = Atom(1, 'CA', '', 'A', 'ALA', 1, '', x=1., y=2., z=2., 
                    occupancy=1., b_factor=0., element=None, mass=None)
        atom2 = Atom()
        
        assert_almost_equals(3.0, atom1.distance(atom2))
        assert_almost_equals(3.0, atom2.distance(atom1))
        
    def test_clone(self):
        atom1 = Atom()
        atom2 = atom1.clone()
        
        assert atom1 == atom2
        assert not atom1 != atom2
        
        atom1.name = 'C'
        
        assert atom1 != atom2
    
    def test_to_string(self):
        atom1 = Atom(1, 'CA', '', 'A', 'ALA', 1, '', x=1., y=2., z=2., 
                    occupancy=1., b_factor=0., element=None, mass=None)
        assert "  CA   1.000   2.000   2.000" == str(atom1)


class TestHetAtom:
    
    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_create_atom(self):
        atom = HetAtom()
        assert atom.__class__.__name__ == "HetAtom"
