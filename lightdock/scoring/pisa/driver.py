"""PISA potential
Shruthi Viswanath, D. V. S. Ravikant, and Ron Elber.
Improving ranking of models for protein complexes with side chain modeling and atomic potentials.
Proteins 81(4) 2012

Original binary downloaded from http://clsb.ices.utexas.edu/web/dock_details.html
"""

import os
from lightdock.structure.model import DockingModel
from lightdock.scoring.functions import ModelAdapter, ScoringFunction
from lightdock.scoring.pisa.cython.cpisa import calculate_pisa
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF


class PISAPotential(object):
    """Loads PISA potential information"""
    # PISA constants
    num_atom_types = 32
    max_distance = 8.
    min_distance = 2.
    num_bins = 3

    def __init__(self):
        data_path = os.path.dirname(os.path.realpath(__file__)) + '/data/'
        self.pisa_energy = PISAPotential._read_parameters(data_path + 'pisa.params')

    @staticmethod
    def _read_parameters(data_file_name):
        """Reads PISA beans"""
        params_atom_potentials = [[[0 for k in range(PISAPotential.num_bins)]
                                   for j in range(PISAPotential.num_atom_types)]
                                  for i in range(PISAPotential.num_atom_types)]
        params_file = open(data_file_name)
        lines = params_file.readlines()
        params_file.close()
        line_count = 0
        for i in range(PISAPotential.num_atom_types):
            for j in range(i, PISAPotential.num_atom_types):
                for r in range(PISAPotential.num_bins):
                    params_atom_potentials[i][j][r] = float(lines[line_count])
                    line_count += 1
        return params_atom_potentials

    @staticmethod
    def get_atom_type(atom_name, res_name):
        if res_name == "LYS" and atom_name == "NZ":
            return 1
        if atom_name == "N":
            return 2
        if atom_name == "C" or (res_name == "ASN" and atom_name == "CG") or (res_name == "GLN" and atom_name == "CD"):
            return 3
        if atom_name == "O" or (res_name == "ASN" and atom_name == "OD1") or (
                        res_name == "GLN" and atom_name == "OE1") or (res_name == "GLN" and atom_name == "OE"):
            return 4
        if atom_name == "CA" and res_name != "PRO":
            return 5
        if ((res_name == "ALA" and atom_name == "CB")
            or (res_name == "ILE" and (atom_name == "CG2" or atom_name == "CD1"))
            or (res_name == "LEU" and (atom_name == "CD1" or atom_name == "CD2"))
            or (res_name == "LEU" and (atom_name == "CD" or atom_name == "CE"))
            or (res_name == "THR" and atom_name == "CG2")
            or (res_name == "VAL" and (atom_name == "CG1" or atom_name == "CG2"))):
            return 6
        if ((res_name == "ARG" and (atom_name == "CB" or atom_name == "CG"))
            or (res_name == "ASN" and atom_name == "CB")
            or (res_name == "GLN" and (atom_name == "CB" or atom_name == "CG"))
            or (res_name == "GLU" and atom_name == "CB")
            or (res_name == "HIS" and atom_name == "CB")
            or (res_name == "ILE" and (atom_name == "CB" or atom_name == "CG1"))
            or (res_name == "LEU" and (atom_name == "CB" or atom_name == "CG"))
            or (res_name == "LYS" and (atom_name == "CB" or atom_name == "CG" or atom_name == "CD"))
            or (res_name == "MET" and atom_name == "CB")
            or (res_name == "PHE" and atom_name == "CB")
            or (res_name == "PRO" and (atom_name == "CB" or atom_name == "CG"))
            or (res_name == "TRP" and atom_name == "CB")
            or (res_name == "TYR" and atom_name == "CB")
            or (res_name == "VAL" and atom_name == "CB")):
            return 7
        if ((res_name == "PHE" and (
                                        atom_name == "CG" or atom_name == "CD1" or atom_name == "CD2" or atom_name == "CE1" or atom_name == "CE2" or atom_name == "CZ"))
            or (res_name == "TRP" and (
                                    atom_name == "CE3" or atom_name == "CZ2" or atom_name == "CZ3" or atom_name == "CH2"))
            or (res_name == "TYR" and (
                                        atom_name == "CG" or atom_name == "CD1" or atom_name == "CD2" or atom_name == "CE1" or atom_name == "CE2"))
            or (res_name == "TYR" and (
                                        atom_name == "CC" or atom_name == "CD" or atom_name == "CE" or atom_name == "CG" or atom_name == "CH"))):
            return 8
        if res_name == "TYR" and (atom_name == "CZ" or atom_name == "CF"):
            return 9
        if ((res_name == "SER" and atom_name == "OG")
            or (res_name == "THR" and atom_name == "OG1")
            or (res_name == "TYR" and atom_name == "OH")):
            return 10
        if res_name == "TRP" and (atom_name == "CG" or atom_name == "CD2"):
            return 11
        if res_name == "TRP" and (atom_name == "CD1" or atom_name == "CE2"):
            return 12
        if res_name == "TRP" and atom_name == "NE1":
            return 13
        if res_name == "MET" and atom_name == "CG":
            return 14
        if res_name == "MET" and (atom_name == "SD" or atom_name == "S" or atom_name == "SE"):
            return 15
        if res_name == "MET" and atom_name == "CE":
            return 16
        if (res_name == "LYS" and atom_name == "CE") or (res_name == "LYS" and atom_name == "CZ"):
            return 17
        if (res_name == "SER" and atom_name == "CB") or (res_name == "THR" and atom_name == "CB"):
            return 18
        if res_name == "PRO" and (atom_name == "CD" or atom_name == "CA"):
            return 19
        if res_name == "CYS" and atom_name == "CB":
            return 20
        if res_name == "CYS" and (atom_name == "SG" or atom_name == "S" or atom_name == "SE"):
            return 21
        if res_name == "HIS" and (atom_name == "CG" or atom_name == "CD2"):
            return 22
        if res_name == "HIS" and (atom_name == "ND1" or atom_name == "NE2"):
            return 23
        if res_name == "HIS" and atom_name == "CE1":
            return 24
        if res_name == "ARG" and atom_name == "CD":
            return 25
        if res_name == "ARG" and atom_name == "NE":
            return 26
        if res_name == "ARG" and atom_name == "CZ":
            return 27
        if res_name == "ARG" and (atom_name == "NH1" or atom_name == "NH2"):
            return 28
        if (res_name == "ASN" and atom_name == "ND2") or (res_name == "GLN" and atom_name == "NE2"):
            return 29
        if (res_name == "ASP" and atom_name == "CB") or (res_name == "GLU" and atom_name == "CG"):
            return 30
        if ((res_name == "ASP" and atom_name == "CG")
            or (res_name == "GLU" and atom_name == "CD")
            or (res_name == "GLU" and atom_name == "CD1")):
            return 31
        if (atom_name == "OXT"
            or (res_name == "ASP" and (atom_name == "OD1" or atom_name == "OD2"))
            or (res_name == "GLU" and (atom_name == "OE1" or atom_name == "OE2"))
            or (res_name == "GLU" and (atom_name == "OE11" or atom_name == "OE21"))
            or (res_name == "GLU" and atom_name == "OE") or atom_name == "OX2"):
            return 32
        return -1


class PISAAdapter(ModelAdapter):
    """Adapts a given Complex to a DockingModel object suitable for this
    PISA scoring function.
    """
    def _get_docking_model(self, molecule, restraints):
        """Builds a suitable docking model for this scoring function"""
        pisa_objects = []
        coordinates = []
        parsed_restraints = {}
        for atom_index, atom in enumerate(molecule.atoms):
            atom_type = PISAPotential.get_atom_type(atom.name, atom.residue_name)
            if atom_type != -1:
                atom.pisa_type = atom_type
                pisa_objects.append(atom)
                coordinates.append([atom.x, atom.y, atom.z])
                res_id = "%s.%s.%s" % (atom.chain_id, atom.residue_name, str(atom.residue_number))
                if restraints and res_id in restraints:
                    try:
                        parsed_restraints[res_id].append(atom_index)
                    except:
                        parsed_restraints[res_id] = [atom_index]
        try:
            return DockingModel(pisa_objects, molecule.copy_coordinates(), parsed_restraints, n_modes=molecule.n_modes.copy())
        except AttributeError:
            return DockingModel(pisa_objects, molecule.copy_coordinates(), parsed_restraints)


class PISA(ScoringFunction):
    """Implements PISA scoring function"""
    def __init__(self, weight=1.0):
        super(PISA, self).__init__(weight)
        self.potential = PISAPotential()

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        energy, interface_receptor, interface_ligand = calculate_pisa(receptor, receptor_coordinates, 
                                                                      ligand, ligand_coordinates, 
                                                                      self.potential.pisa_energy,
                                                                      DEFAULT_CONTACT_RESTRAINTS_CUTOFF)
        perc_receptor_restraints = ScoringFunction.restraints_satisfied(receptor.restraints, interface_receptor)
        perc_ligand_restraints = ScoringFunction.restraints_satisfied(ligand.restraints, interface_ligand)
        return (energy + perc_receptor_restraints * energy + perc_ligand_restraints * energy) * self.weight


# Needed to dynamically load the scoring functions from command line
DefinedScoringFunction = PISA
DefinedModelAdapter = PISAAdapter
