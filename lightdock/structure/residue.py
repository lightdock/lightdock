"""Module to package a complex residue representation and operations"""

import numpy as np
from scipy.spatial import cKDTree
from lightdock.error.lightdock_errors import (
    ResidueNonStandardError,
    SideChainError,
    BackboneError,
)
from lightdock.structure.atom import Atom


backbone = ["N", "CA", "C", "O"]

sidechain = {
    "ALA": ["CB"],
    "ARG": ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["CB", "CG", "OD1", "ND2"],
    "ASP": ["CB", "CG", "OD1", "OD2"],
    "CYS": ["CB", "SG"],
    "GLU": ["CB", "CG", "CD", "OE1", "OE2"],
    "GLN": ["CB", "CG", "CD", "OE1", "NE2"],
    "GLY": [],
    "HIS": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["CB", "CG1", "CG2", "CD1"],
    "LEU": ["CB", "CG", "CD1", "CD2"],
    "LYS": ["CB", "CG", "CD", "CE", "NZ"],
    "MET": ["CB", "CG", "SD", "CE"],
    "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["CB", "CG", "CD"],
    "SER": ["CB", "OG"],
    "THR": ["CB", "OG1", "CG2"],
    "TRP": ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["CB", "CG1", "CG2"],
}

bond_length = {
    ("CA", "N"): 1.46,
    ("CA", "C"): 1.52,
    ("C", "O"): 1.23,
    ("N", "CA"): 1.46,
    ("C", "CA"): 1.52,
    ("O", "C"): 1.23,
}


class Residue(object):
    """Represents a chemical residue in a complex"""

    STANDARD_TYPES = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
    }

    DNA_STANDARD_TYPES = ["DA", "DC", "DI", "DG", "DT"]
    RNA_STANDARD_TYPES = ["A", "C", "G", "U", "I"]

    MODIFIED_TYPES = {"CYX": "C", "HIP": "H", "HID": "H", "HIE": "H"}

    DUMMY_TYPES = ["MMB", "DUM"]

    def __init__(
        self,
        residue_name,
        residue_number,
        residue_insertion="",
        atoms=None,
        residue_index=0,
    ):
        """Creates a new residue"""
        self.name = residue_name.upper()
        self.number = residue_number
        self.insertion = residue_insertion.upper().strip()
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
        return Residue(
            self.name,
            self.number,
            self.insertion,
            [atom.clone() for atom in self.atoms],
            self.index,
        )

    def is_standard(self):
        """Checks if residue is standard"""
        return self.name in list(Residue.STANDARD_TYPES.keys())

    def is_protein(self):
        """Checks if residue is protein"""
        return self.name in list(Residue.STANDARD_TYPES.keys()) or self.name in list(
            Residue.MODIFIED_TYPES.keys()
        )

    def is_nucleic(self):
        """Check if residue is Deoxyribonucleotide or Ribonucleotide"""
        return self.name in Residue.DNA_STANDARD_TYPES + Residue.RNA_STANDARD_TYPES

    def is_dummy(self):
        """Checks if residue is a dummy bead"""
        return self.name in Residue.DUMMY_TYPES

    def set_backbone_and_sidechain(self):
        """Classifies the atoms in backbone or side-chain"""
        if self.is_standard():
            self.backbone = [atom for atom in self.atoms if atom.name in backbone]
            self.sidechain = [
                atom for atom in self.atoms if atom.name in sidechain[self.name]
            ]
        else:
            self.backbone = []
            self.sidechain = []

    def check(self):
        """Check if the residue has all the backbone and sidechain atoms, ignore dummy beads"""
        if self.is_standard():
            backbone_correct = set([a.name for a in self.backbone]) == set(backbone)
            if not backbone_correct:
                raise BackboneError(
                    f"Incomplete backbone for residue {self.name}.{self.number}{self.insertion}"
                )

            sd_correct = set([a.name for a in self.sidechain]) == set(
                sidechain[self.name]
            )
            if not sd_correct:
                raise SideChainError(
                    f"Incomplete sidechain for residue {self.name}.{self.number}{self.insertion}"
                )

            return True
        else:
            if not (self.is_dummy() or self.is_nucleic()):
                raise ResidueNonStandardError(
                    f"Can not check non-standard residue {self.name}.{self.number}{self.insertion}"
                )

    def __eq__(self, other):
        """Compares two residues for equality."""
        return (
            self.number == other.number
            and self.name == other.name
            and self.insertion == other.insertion
        )

    def __ne__(self, other):
        """Compares two residues for unequality"""
        return not self.__eq__(other)

    def get_atom(self, atom_name):
        """Gets the atom identified by atom_name"""
        for atom in self.atoms:
            if atom.name == atom_name:
                return atom
        return None

    def get_chain(self):
        """Gets the chain ID"""
        if self.atoms:
            return self.atoms[0].chain_id
        else:
            return None

    def get_calpha(self):
        """Get the Calpha atom"""
        return self.get_atom("CA")

    def get_central_atom(self):
        """Calculates center of coordiantes of residue and returns closest atom"""
        coordinates = np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])
        centroid = coordinates.mean(axis=0)
        min_dist, min_dist_idx = cKDTree(coordinates).query(centroid, 1)
        return self.atoms[min_dist_idx]

    def get_non_hydrogen_atoms(self):
        return [atom for atom in self.atoms if not atom.is_hydrogen()]

    @staticmethod
    def dummy(x=0.0, y=0.0, z=0.0):
        """Creates a dummy residue with DUM residue name"""
        atom = Atom(atom_name="CA", residue_name="DUM", x=x, y=y, z=z)
        return Residue(residue_name="DUM", residue_number=0, atoms=[atom])

    def __str__(self):
        if len(self.atoms):
            representation = []
            for atom in self.atoms:
                representation.append(
                    f"{self.name}.{self.number}{self.insertion}  {str(atom)}"
                )
            return "\n".join(representation)
        else:
            return f"{self.name}.{self.number}{self.insertion}"


    def full_name(self):
        """Get the full id of this residue"""
        return f"{self.name}.{self.number}{self.insertion}"


class AminoAcid(Residue):
    """Amino acid residue type"""

    pass


class Cofactor(Residue):
    """Non-protein chemical compound type"""

    pass


class Ion(Residue):
    """Charged chemical compound type"""

    pass
