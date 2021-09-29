"""MJ potentials scoring functions

Miyazawa S, Jernigan RL. Estimation of effective interresidue contact energies from protein crystal
structures--quasi-chemical approximation. Macromolecules 1985;18:534-552.
Miyazawa S, Jernigan RL. Residue-residue potentials with a favorable contact pair term and an unfavorable
high packing density term, for simulation and threading. J Mol Biol 1996;256:623-644.
Miyazawa S, Jernigan RL. Self-consistent estimation of interresidue protein contact energies based on an
equilibrium mixture approximation of residues. Proteins 1999;34:49-68.
"""

import os
from lightdock.error.lightdock_errors import PotentialsParsingError
from lightdock.structure.model import DockingModel
from lightdock.structure.space import SpacePoints
from lightdock.scoring.functions import ModelAdapter, ScoringFunction
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF


class MJPotential(object):
    """Loads MJ potentials information"""

    residues = {
        "LEU": 0,
        "PHE": 1,
        "ILE": 2,
        "MET": 3,
        "VAL": 4,
        "TRP": 5,
        "CYS": 6,
        "TYR": 7,
        "HIS": 8,
        "ALA": 9,
        "THR": 10,
        "GLY": 11,
        "PRO": 12,
        "ARG": 13,
        "GLN": 14,
        "SER": 15,
        "ASN": 16,
        "GLU": 17,
        "ASP": 18,
        "LYS": 19,
    }

    def __init__(self, data_file="MJ_potentials.dat"):
        data_path = os.path.dirname(os.path.realpath(__file__)) + "/data/"
        self.potentials = MJPotential._read_potentials(data_path + data_file)

    @staticmethod
    def _read_potentials(data_file_name):
        """Reads a given potentials file.

        Identifies each potential matrix by '>' and '<'
        """
        data_file = open(data_file_name)
        data = data_file.readlines()
        data_file.close()
        potentials = {}

        try:
            current_potential = ""
            for line in data:
                line = line.rstrip()
                if line:
                    if line[0] == "<":
                        current_potential = ""

                    if current_potential:
                        potentials[current_potential].append(
                            [float(x) for x in line.split()]
                        )

                    if line[0] == ">":
                        current_potential = line[1:]
                        potentials[current_potential] = []

        except Exception as e:
            raise PotentialsParsingError(
                "Error parsing %s file. Details: %s" % (data_file_name, str(e))
            )

        return potentials


class MJ3hAdapter(ModelAdapter):
    """Adapts a given Complex to a DockingModel object suitable for this
    MJ3h scoring function.
    """

    def _get_docking_model(self, molecule, restraints):
        """Builds a suitable docking model for this scoring function"""
        list_of_coordinates = []
        not_considered_atoms = ["O", "C", "N", "H"]
        residues_to_remove = []
        residues = [residue for chain in molecule.chains for residue in chain.residues]
        parsed_restraints = {}
        for structure in range(molecule.num_structures):
            coordinates = []
            for res_index, residue in enumerate(residues):
                c_x = 0.0
                c_y = 0.0
                c_z = 0.0
                count = 0
                chain_id = ""
                for atom in residue.atoms:
                    if atom.name not in not_considered_atoms:
                        c_x += molecule.atom_coordinates[structure][atom.index][0]
                        c_y += molecule.atom_coordinates[structure][atom.index][1]
                        c_z += molecule.atom_coordinates[structure][atom.index][2]
                        count += 1
                        chain_id = atom.chain_id
                if count:
                    count = float(count)
                    coordinates.append([c_x / count, c_y / count, c_z / count])
                    res_id = (
                        f"{chain_id}.{residue.name}.{residue.number}{residue.insertion}"
                    )
                    if restraints and res_id in restraints:
                        parsed_restraints[res_id] = []
                else:
                    residues_to_remove.append(res_index)

            if len(residues_to_remove):
                residues = [
                    residue
                    for res_index, residue in enumerate(residues)
                    if res_index not in residues_to_remove
                ]

            list_of_coordinates.append(SpacePoints(coordinates))

        return DockingModel(residues, list_of_coordinates, parsed_restraints)


class MJ3h(ScoringFunction):
    """Implements MJ3h potential"""

    potentials_dict = {
        "LEU": 0,
        "PHE": 1,
        "ILE": 2,
        "MET": 3,
        "VAL": 4,
        "TRP": 5,
        "CYS": 6,
        "TYR": 7,
        "HIS": 8,
        "ALA": 9,
        "THR": 10,
        "GLY": 11,
        "PRO": 12,
        "ARG": 13,
        "GLN": 14,
        "SER": 15,
        "ASN": 16,
        "GLU": 17,
        "ASP": 18,
        "LYS": 19,
    }

    max_distance_cutoff = 42.25  # 6.5^2
    min_distance_cutoff = 6.2  # 2.5^2

    def __init__(self, weight=1.0, penalization=3.0):
        super(MJ3h, self).__init__(weight, anm_support=False)
        self.penalization = penalization
        self.potential = MJPotential()
        self.potentials = self.potential.potentials["MJ3h"]
        self.cutoff = (
            DEFAULT_CONTACT_RESTRAINTS_CUTOFF * DEFAULT_CONTACT_RESTRAINTS_CUTOFF
        )

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        """Calculates the MJ3h potential taking into account the contacts between receptor
        and ligand. Receptor and ligand are DockingModel objects.
        """
        energy = 0.0
        interface_receptor = []
        interface_ligand = []
        for index_rec, res1 in enumerate(receptor.objects):
            for index_lig, res2 in enumerate(ligand.objects):
                [x1, y1, z1] = receptor_coordinates[index_rec]
                [x2, y2, z2] = ligand_coordinates[index_lig]
                distance = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                if distance < MJ3h.max_distance_cutoff:
                    if distance < MJ3h.min_distance_cutoff:
                        energy += self.penalization
                    else:
                        try:
                            i_rec = MJ3h.potentials_dict[res1.name]
                            i_lig = MJ3h.potentials_dict[res2.name]
                            energy += self.potentials[i_rec][i_lig]
                        except:
                            pass
                    if distance <= self.cutoff:
                        interface_receptor.append(index_rec)
                        interface_ligand.append(index_lig)
        interface_receptor = set(interface_receptor)
        interface_ligand = set(interface_ligand)
        energy *= -1.0
        perc_receptor_restraints = ScoringFunction.restraints_satisfied(
            receptor.restraints, interface_receptor
        )
        perc_ligand_restraints = ScoringFunction.restraints_satisfied(
            ligand.restraints, interface_ligand
        )
        return (
            energy + perc_receptor_restraints * energy + perc_ligand_restraints * energy
        ) * self.weight


# Needed to dynamically load the scoring functions from command line
DefinedScoringFunction = MJ3h
DefinedModelAdapter = MJ3hAdapter
