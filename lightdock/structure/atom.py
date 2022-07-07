"""Module to package atom representation and operations"""

from lightdock.mathutil.cython.cutil import distance as cdistance
from lightdock.error.lightdock_errors import AtomError


class Atom(object):
    """Represents a chemical atom"""

    BACKBONE_ATOMS = ["CA", "C", "N", "O"]
    RECOGNIZED_ELEMENTS = ["C", "N", "O", "H", "S", "P", "CL", "MG", "FE", "PB", "SE", "F"]
    MASSES = {
        "H": 1.007825,
        "C": 12.01,
        "F": 18.9984,
        "O": 15.9994,
        "N": 14.0067,
        "S": 31.972071,
        "P": 30.973762,
        "CL": 35.45,
        "MG": 24.3050,
        "FE": 55.845,
        "PB": 207.2,
        "SE": 78.96,
    }

    def __init__(
        self,
        atom_number=99999,
        atom_name="H",
        atom_alternative="",
        chain_id="",
        residue_name="",
        residue_number=9999,
        residue_insertion="",
        x=0.0,
        y=0.0,
        z=0.0,
        occupancy=1.0,
        b_factor=0.0,
        element=None,
        mass=None,
        atom_index=0,
    ):
        """Creates a new atom.

        Mass will be selected depending upon atom element. By default, creates a
        regular hydrogen atom. Index can be used to quickly identify an atom or to use it
        as an index in an external data structure, i.e., a coordinates matrix.
        """
        self.number = atom_number
        self.name = atom_name
        self.alternative = atom_alternative
        self.chain_id = chain_id
        self.residue_name = residue_name
        self.residue_number = residue_number
        self.residue_insertion = residue_insertion.strip()
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.b_factor = b_factor
        if element:
            try:
                if element not in Atom.RECOGNIZED_ELEMENTS:
                    raise AtomError(
                        "Not recognized element '%s' for atom %s."
                        % (element, self.name)
                    )
                self.element = element
            except AtomError:
                self._assign_element()
        else:
            self._assign_element()
        if mass:
            self.mass = mass
        else:
            self.mass = Atom.MASSES[self.element]
        self.index = atom_index

    def _assign_element(self):
        """Assigns an element to an atom depending on its name"""
        atom_element = self.name[:2]
        if atom_element in Atom.RECOGNIZED_ELEMENTS:
            self.element = atom_element
        else:
            atom_element = self.name[0]
            if atom_element in Atom.RECOGNIZED_ELEMENTS:
                self.element = atom_element
            else:
                raise AtomError(
                    "Not recognized element '%s' for atom %s."
                    % (atom_element, self.name)
                )

    def distance(self, other):
        return cdistance(self.x, self.y, self.z, other.x, other.y, other.z)

    def is_hydrogen(self):
        """Checks if this atom is of hydrogen type"""
        return self.element == "H"

    def is_backbone(self):
        """Checks if this atom belongs to the residue backbone"""
        return self.name in Atom.BACKBONE_ATOMS

    def get_coordinates(self):
        """Gets the coordinates vector"""
        return [self.x, self.y, self.z]

    def clone(self):
        """Creates a copy of the current atom"""
        return Atom(
            self.number,
            self.name,
            self.alternative,
            self.chain_id,
            self.residue_name,
            self.residue_number,
            self.residue_insertion,
            self.x,
            self.y,
            self.z,
            self.occupancy,
            self.b_factor,
            self.element,
            self.mass,
            self.index,
        )

    def __eq__(self, other):
        """Compares two atoms for equality.

        Compare by number should be enough, but pdb files usually contain errors
        or numeration can be affected by chain id.
        """
        return (
            self.number == other.number
            and self.name == other.name
            and self.chain_id == other.chain_id
        )

    def __ne__(self, other):
        """Compares two atoms for unequality"""
        return not self.__eq__(other)

    def __str__(self):
        return "%4s%8.3f%8.3f%8.3f" % (self.name, self.x, self.y, self.z)


class HetAtom(Atom):
    """Represents an heterogeneous atom"""

    def __init__(
        self,
        atom_number=99999,
        atom_name="H",
        atom_alternative="",
        chain_id="",
        residue_name=None,
        residue_number=9999,
        residue_insertion="",
        x=0.0,
        y=0.0,
        z=0.0,
        occupancy=1.0,
        b_factor=0.0,
        element=None,
        mass=None,
        atom_index=0,
    ):
        """Creates a new hetatom"""
        super(HetAtom, self).__init__(
            atom_number,
            atom_name,
            atom_alternative,
            chain_id,
            residue_name,
            residue_number,
            residue_insertion,
            x,
            y,
            z,
            occupancy,
            b_factor,
            element,
            mass,
            atom_index,
        )
