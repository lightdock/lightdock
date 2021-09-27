"""DFIRE2 potential scoring function

Yuedong Yang, Yaoqi Zhou. Ab initio folding of terminal segments with secondary structures
reveals the fine difference between two closely related all-atom statistical energy functions.
Protein Science,17:1212-1219(2008)
"""

import os
import numpy as np
from lightdock.structure.model import DockingModel
from lightdock.scoring.functions import ModelAdapter, ScoringFunction
from lightdock.structure.space import SpacePoints
from lightdock.scoring.dfire2.c.cdfire2 import calculate_dfire2
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF

# Potential constants
atom_type_number = 167
bin_number = 30

DFIRE2_ATOM_TYPES = {
    "GLY CA": 40,
    "HIS C": 45,
    "VAL O": 137,
    "GLY O": 42,
    "GLY N": 39,
    "HIS O": 46,
    "HIS N": 43,
    "TRP CE3": 151,
    "GLY C": 41,
    "TRP CE2": 150,
    "LYS NZ": 69,
    "MET C": 80,
    "VAL N": 134,
    "PRO CA": 95,
    "MET O": 81,
    "MET N": 78,
    "SER OG": 126,
    "ARG NH2": 120,
    "VAL C": 136,
    "THR CG2": 133,
    "ALA CB": 4,
    "ALA CA": 1,
    "TRP CG": 146,
    "TRP CA": 142,
    "TRP CB": 145,
    "ALA N": 0,
    "ILE CB": 57,
    "ILE CA": 54,
    "TRP CH2": 154,
    "GLU CA": 20,
    "GLU CB": 23,
    "GLU CD": 25,
    "GLU CG": 24,
    "HIS CG": 48,
    "ASP OD1": 17,
    "HIS CA": 44,
    "CYS N": 5,
    "CYS O": 8,
    "HIS CE1": 51,
    "TYR CG": 160,
    "TYR CA": 156,
    "TYR CB": 159,
    "CYS C": 7,
    "ARG CB": 114,
    "LYS C": 63,
    "ARG CG": 115,
    "ARG CD": 116,
    "THR OG1": 132,
    "LYS O": 64,
    "LYS N": 61,
    "SER C": 123,
    "ILE CD1": 60,
    "PRO CB": 98,
    "PRO CD": 100,
    "PRO CG": 99,
    "ARG CZ": 118,
    "SER O": 124,
    "SER N": 121,
    "PHE CD1": 34,
    "PHE CD2": 35,
    "THR CA": 128,
    "HIS CD2": 50,
    "THR CB": 131,
    "PRO C": 96,
    "PRO N": 94,
    "PRO O": 97,
    "PHE CA": 29,
    "MET CE": 85,
    "MET CG": 83,
    "MET CA": 79,
    "ILE C": 55,
    "MET CB": 82,
    "TRP CD2": 148,
    "TRP CD1": 147,
    "GLN CD": 107,
    "ILE CG1": 58,
    "ILE CG2": 59,
    "PHE CE2": 37,
    "PHE CE1": 36,
    "GLU OE1": 26,
    "GLU OE2": 27,
    "ASP CG": 16,
    "ASP CB": 15,
    "ASP CA": 12,
    "THR O": 130,
    "THR N": 127,
    "SER CA": 122,
    "SER CB": 125,
    "PHE CG": 33,
    "GLU O": 22,
    "GLU N": 19,
    "PHE CB": 32,
    "VAL CG1": 139,
    "GLU C": 21,
    "ILE O": 56,
    "ILE N": 53,
    "GLN CA": 102,
    "GLN CB": 105,
    "ASN C": 88,
    "VAL CG2": 140,
    "TRP CZ2": 152,
    "TRP CZ3": 153,
    "PHE CZ": 38,
    "TRP O": 144,
    "TRP N": 141,
    "LEU CB": 74,
    "GLN N": 101,
    "GLN O": 104,
    "LEU O": 73,
    "GLN C": 103,
    "TRP C": 143,
    "HIS CB": 47,
    "GLN NE2": 109,
    "LEU CD2": 77,
    "ASP OD2": 18,
    "LEU CD1": 76,
    "VAL CA": 135,
    "ASN OD1": 92,
    "ALA O": 3,
    "MET SD": 84,
    "ALA C": 2,
    "THR C": 129,
    "TYR CD1": 161,
    "ARG NH1": 119,
    "TYR CD2": 162,
    "ASN ND2": 93,
    "TRP NE1": 149,
    "HIS ND1": 49,
    "LEU C": 72,
    "ASN O": 89,
    "ASN N": 86,
    "ASP C": 13,
    "LEU CA": 71,
    "ASP O": 14,
    "ASP N": 11,
    "CYS CB": 9,
    "LEU N": 70,
    "LEU CG": 75,
    "CYS CA": 6,
    "TYR OH": 166,
    "ASN CA": 87,
    "ASN CB": 90,
    "ASN CG": 91,
    "TYR CE2": 164,
    "ARG C": 112,
    "TYR CE1": 163,
    "HIS NE2": 52,
    "ARG O": 113,
    "ARG N": 110,
    "TYR C": 157,
    "GLN CG": 106,
    "ARG CA": 111,
    "TYR N": 155,
    "TYR O": 158,
    "CYS SG": 10,
    "TYR CZ": 165,
    "ARG NE": 117,
    "VAL CB": 138,
    "LYS CB": 65,
    "LYS CA": 62,
    "PHE C": 30,
    "LYS CG": 66,
    "LYS CE": 68,
    "LYS CD": 67,
    "GLN OE1": 108,
    "PHE N": 28,
    "PHE O": 31,
}


class DFIRE2Potential(object):
    """Loads DFIRE2 potentials information"""

    def __init__(self):
        data_path = os.path.dirname(os.path.realpath(__file__)) + "/data/"
        self.energy = np.load(data_path + "dfire2_energies.npy").ravel()


class DFIRE2Object(object):
    def __init__(self, residue_index, atom_index):
        self.residue_index = residue_index
        self.atom_index = atom_index


class DFIRE2Adapter(ModelAdapter, DFIRE2Potential):
    """Adapts a given Complex to a DockingModel object suitable for this
    DFIRE2 scoring function.
    """

    def _get_docking_model(self, molecule, restraints):
        """Builds a suitable docking model for this scoring function"""
        objects = []
        coordinates = []
        parsed_restraints = {}
        atom_index = 0
        for residue in molecule.residues:
            for rec_atom in residue.atoms:
                rec_atom_type = rec_atom.residue_name + " " + rec_atom.name
                if rec_atom_type in DFIRE2_ATOM_TYPES:
                    objects.append(
                        DFIRE2Object(residue.number, DFIRE2_ATOM_TYPES[rec_atom_type])
                    )
                    coordinates.append([rec_atom.x, rec_atom.y, rec_atom.z])
                    # Restraints support
                    res_id = f"{rec_atom.chain_id}.{residue.name}.{residue.number}{residue.insertion}"
                    if restraints and res_id in restraints:
                        try:
                            parsed_restraints[res_id].append(atom_index)
                        except:
                            parsed_restraints[res_id] = [atom_index]
                    atom_index += 1
        try:
            return DockingModel(
                objects,
                SpacePoints(coordinates),
                parsed_restraints,
                n_modes=molecule.n_modes.copy(),
            )
        except AttributeError:
            return DockingModel(objects, SpacePoints(coordinates), parsed_restraints)


class DFIRE2(ScoringFunction):
    """Implements DFIRE2 potential"""

    def __init__(self, weight=1.0):
        super(DFIRE2, self).__init__(weight)
        self.cached = False
        self.potential = DFIRE2Potential()

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        if not self.cached:
            self.res_index = []
            self.atom_index = []
            for o in receptor.objects:
                self.res_index.append(o.residue_index)
                self.atom_index.append(o.atom_index)
            last = self.res_index[-1]
            for o in ligand.objects:
                self.res_index.append(o.residue_index + last)
                self.atom_index.append(o.atom_index)
            self.res_index = np.array(self.res_index, dtype=np.int32)
            self.atom_index = np.array(self.atom_index, dtype=np.int32)
            self.molecule_length = len(self.res_index)
            self.cached = True
        return self.evaluate_energy(
            receptor, receptor_coordinates, ligand, ligand_coordinates
        )

    def evaluate_energy(
        self, receptor, receptor_coordinates, ligand, ligand_coordinates
    ):
        coordinates = np.append(
            receptor_coordinates.coordinates, ligand_coordinates.coordinates
        ).reshape((-1, 3))
        energy, interface_receptor, interface_ligand = calculate_dfire2(
            self.res_index,
            self.atom_index,
            coordinates,
            self.potential.energy,
            self.molecule_length,
            DEFAULT_CONTACT_RESTRAINTS_CUTOFF,
        )

        # Code to consider contacts in the interface
        perc_receptor_restraints = ScoringFunction.restraints_satisfied(
            receptor.restraints, set(interface_receptor)
        )
        perc_ligand_restraints = ScoringFunction.restraints_satisfied(
            ligand.restraints, set(interface_ligand)
        )
        return (
            energy + perc_receptor_restraints * energy + perc_ligand_restraints * energy
        ) * self.weight


# Needed to dynamically load the scoring functions from command line
DefinedScoringFunction = DFIRE2
DefinedModelAdapter = DFIRE2Adapter
