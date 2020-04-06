"""Module to package a complex residue representation and operations"""

from lightdock.error.lightdock_errors import ResidueNonStandardError, SideChainError, BackboneError
from lightdock.structure.atom import Atom


backbone = ["N", "CA", "C", "O"]

sidechain = {"ALA": ['CB'],
             "ARG": ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
             "ASN": ['CB', 'CG', 'OD1', 'ND2'],
             "ASP": ['CB', 'CG', 'OD1', 'OD2'],
             "CYS": ['CB', 'SG'],
             "GLU": ['CB', 'CG', 'CD', 'OE1', 'OE2'],
             "GLN": ['CB', 'CG', 'CD', 'OE1', 'NE2'],
             "GLY": [],
             "HIS": ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
             "ILE": ['CB', 'CG1', 'CG2', 'CD1'],
             "LEU": ['CB', 'CG', 'CD1', 'CD2'],
             "LYS": ['CB', 'CG', 'CD', 'CE', 'NZ'],
             "MET": ['CB', 'CG', 'SD', 'CE'],
             "PHE": ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
             "PRO": ['CB', 'CG', 'CD'],
             "SER": ['CB', 'OG'],
             "THR": ['CB', 'OG1', 'CG2'],
             "TRP": ['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
             "TYR": ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
             "VAL": ['CB', 'CG1', 'CG2']}

bond_length = {('CA', 'N'): 1.46, ('CA', 'C'): 1.52, ('C', 'O'): 1.23,
               ('N', 'CA'): 1.46, ('C', 'CA'): 1.52, ('O', 'C'): 1.23}


class Residue(object):
    """Represents a chemical residue in a complex"""
    
    STANDARD_TYPES = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                      'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                      'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                      'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
    
    DNA_STANDARD_TYPES = ['DA', 'DC', 'DI', 'DG', 'DT']
    RNA_STANDARD_TYPES = ['A', 'C', 'G', 'U', 'I']

    MODIFIED_TYPES = {'CYX': 'C', 'HIP': 'H', 'HID': 'H', 'HIE': 'H'}

    DUMMY_TYPES = ['MMB', 'DUM']
    
    def __init__(self, residue_name, residue_number, atoms=None, residue_index=0):
        """Creates a new residue"""
        self.name = residue_name.upper()
        self.number = residue_number
        if atoms:
            self.atoms = atoms
            self.set_backbone_and_sidechain()
        else:
            self.atoms = []
            self.backbone = []
            self.sidechain = []
        self.index = residue_index

    def clone(self):
        """Creates a copy of the current residue"""
        return Residue(self.name,
                       self.number,
                       [atom.clone() for atom in self.atoms],
                       self.index)

    def is_standard(self):
        """Checks if residue is standard"""
        return self.name in list(Residue.STANDARD_TYPES.keys())

    def is_nucleic(self):
        """Check if residue is Deoxyribonucleotide or Ribonucleotide"""
        return self.name in Residue.DNA_STANDARD_TYPES+Residue.RNA_STANDARD_TYPES

    def is_dummy(self):
        """Checks if residue is a dummy bead"""
        return self.name in Residue.DUMMY_TYPES

    def set_backbone_and_sidechain(self):
        """Classifies the atoms in backbone or side-chain"""
        if self.is_standard():
            self.backbone = [atom for atom in self.atoms if atom.name in backbone]
            self.sidechain = [atom for atom in self.atoms if atom.name in sidechain[self.name]]
        else:
            self.backbone = []
            self.sidechain = []

    def check(self):
        """Check if the residue has all the backbone and sidechain atoms, ignore dummy beads"""
        if self.is_standard():
            backbone_correct = sorted([a.name for a in self.backbone]) == sorted(backbone)
            if not backbone_correct:
                raise BackboneError("Incomplete backbone for residue %s.%s" % (self.name, self.number))
            
            sd_correct = sorted([a.name for a in self.sidechain]) == sorted(sidechain[self.name])
            if not sd_correct:
                raise SideChainError("Incomplete sidechain for residue %s.%s" % (self.name, self.number))
            
            return True
        else:
            if not (self.is_dummy() or self.is_nucleic()):
                raise ResidueNonStandardError("Can not check non-standard residue %s.%s" % (self.name, self.number))

    def __eq__(self, other):
        """Compares two residues for equality."""
        return self.number == other.number and self.name == other.name

    def __ne__(self, other):
        """Compares two residues for unequality"""
        return not self.__eq__(other)

    def get_atom(self, atom_name):
        """Gets the atom identified by atom_name"""
        for atom in self.atoms:
            if atom.name == atom_name:
                return atom
        return None

    def get_calpha(self):
        """Get the Calpha atom"""
        return self.get_atom('CA')

    def get_non_hydrogen_atoms(self):
        return [atom for atom in self.atoms if not atom.is_hydrogen()]

    def mutate_side_chain(self, rotamer):
        """Moves this residue's side chain using the rotamer angles"""
        pass

    @staticmethod
    def dummy(x=0., y=0., z=0.):
        atom = Atom(atom_name='CA', residue_name='DUM', x=x, y=y, z=z)
        return Residue(residue_name='DUM', residue_number=0, atoms=[atom])

    def __str__ (self):
        if len(self.atoms):
            representation = []
            for atom in self.atoms:
                representation.append("%s.%s  %s" % (self.name, self.number, str(atom)))
            return '\n'.join(representation)
        else:    
            return "%s.%s" % (self.name, self.number) 


class AminoAcid(Residue):
    """Amino acid residue type"""
    pass


class Cofactor(Residue):
    """Non-protein chemical compound type"""
    pass


class Ion(Residue):
    """Charged chemical compound type"""
    pass
