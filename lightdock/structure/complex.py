"""Module to package a protein complex"""

from lightdock.structure.space import SpacePoints


class Complex(object):
    """Represents a molecular complex"""
    def __init__(self, chains, atoms=None, residues=None, structure_file_name='', structures=None, representative_id=0):
        """Creates a new complex that can deal with multiple coordinates for a given atom"""
        self.chains = chains
        # Set atoms at the upper level for fast indexing
        if atoms:
            self.atoms = atoms
        else:
            self.atoms = [atom for chain in self.chains
                          for residue in chain.residues for atom in residue.atoms]
        for atom_index, atom in enumerate(self.atoms):
            atom.index = atom_index

        # Same for residues
        if residues:
            self.residues = residues
        else:
            self.residues = [residue for chain in self.chains for residue in chain.residues]
        for residue_index, residue in enumerate(self.residues):
            residue.index = residue_index

        if structures:
            self.num_structures = len(structures)
            self.structure_file_names = [structure['file_name'] for structure in structures]
            self.atom_coordinates = [SpacePoints([[atom.x, atom.y, atom.z]
                                                  for atom in structure['atoms']]) for structure in structures]
        else:
            self.num_structures = 1
            self.structure_file_names = [structure_file_name]
            self.atom_coordinates = [SpacePoints([[atom.x, atom.y, atom.z] for atom in self.atoms])]

        self.num_atoms = len(self.atoms)
        self.num_residues = len(self.residues)
        self.representative_id = representative_id

    @staticmethod
    def from_structures(structures, representative_id=0):
        return Complex(structures[representative_id]['chains'],
                       structures[representative_id]['atoms'],
                       structures[representative_id]['residues'],
                       structures[representative_id]['file_name'],
                       structures, representative_id)

    def clone(self):
        """Creates a copy of the current complex"""
        molecule = Complex([chain.clone() for chain in self.chains])
        molecule.num_structures = self.num_structures
        molecule.structure_file_names = self.structure_file_names
        molecule.atom_coordinates = self.copy_coordinates()
        return molecule

    def copy_coordinates(self):
        return [coordinates.clone() for coordinates in self.atom_coordinates]

    def center_of_mass(self, structure=None):
        """Calculates the center of mass"""
        if not structure:
            structure = self.representative_id
        if len(self.atoms):
            total_x = 0.
            total_y = 0.
            total_z = 0.
            total_mass = 0.
            for atom in self.atoms:
                total_x += self.atom_coordinates[structure][atom.index][0] * atom.mass
                total_y += self.atom_coordinates[structure][atom.index][1] * atom.mass
                total_z += self.atom_coordinates[structure][atom.index][2] * atom.mass
                total_mass += atom.mass
            return [total_x/total_mass, total_y/total_mass, total_z/total_mass]
        else:
            return [0., 0., 0.]

    def center_of_coordinates(self, structure=None):
        """Calculates the center of coordinates"""
        if not structure:
            structure = self.representative_id
        atoms = [atom for atom in self.atoms if not atom.is_hydrogen()]
        dimension = len(atoms)
        if dimension:
            total_x = 0.
            total_y = 0.
            total_z = 0.
            for atom in atoms:
                total_x += self.atom_coordinates[structure][atom.index][0]
                total_y += self.atom_coordinates[structure][atom.index][1]
                total_z += self.atom_coordinates[structure][atom.index][2]
            return [total_x/dimension, total_y/dimension, total_z/dimension]
        else:
            return [0., 0., 0.]

    def translate(self, vector):
        """Translates atom coordinates based on vector"""
        for coordinates in self.atom_coordinates:
            coordinates.translate(vector)
        
    def rotate(self, q):
        """Rotates this complex using a quaternion q"""
        for coordinates in self.atom_coordinates:
            coordinates.rotate(q)

    def move_to_origin(self):
        """Moves the structure to the origin of coordinates"""
        translation = [-1*c for c in self.center_of_coordinates()]
        self.translate(translation)
        return translation

    def get_residue(self, chain_id, residue_name, residue_number):
        for chain in self.chains:
            if chain_id == chain.cid:
                for residue in chain.residues:
                    if residue.name == residue_name and int(residue.number) == int(residue_number):
                        return residue
        return None

    def __getitem__(self, item):
        return self.atom_coordinates[item]

    def __setitem__(self, index, item):
        self.atom_coordinates[index] = item

    def __iter__(self):
        for coordinates in self.atom_coordinates:
            yield coordinates

    def __len__(self):
        return self.atom_coordinates.shape[0]

    def representative(self, is_membrane=False):
        coordinates = self.atom_coordinates[self.representative_id]
        if is_membrane:
            transmembrane = []
            for atom_id, atom in enumerate(self.atoms):
                if atom.residue_name != 'MMB':
                    transmembrane.append(coordinates[atom_id])
            return transmembrane
        else:
            return coordinates
