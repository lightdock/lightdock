"""DDNA Scoring function

C. Zhang, S. Liu, Q. Zhu, and Y. Zhou, 'A knowledge-based
energy function for protein-ligand, protein-protein and
protein-DNA complexes', J. Med. Chem. 48, 2325-2335 (2005).
"""

import os

from lightdock.structure.model import DockingModel
from lightdock.scoring.functions import ModelAdapter, ScoringFunction
from lightdock.scoring.ddna.cython.cddna import calculate_ddna
from lightdock.util.logger import LoggingManager
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF
from lightdock.error.lightdock_errors import NotSupportedInScoringError

log = LoggingManager.get_logger("ddna")


atom_map = {
    "CYS": {"N": "N.am", "CA": "C.3", "C": "C.2", "O": "O.2", "CB": "C.3", "SG": "S.3"},
    "MET": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.3",
        "SD": "S.3",
        "CE": "C.3",
    },
    "PHE": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.ar",
        "CD1": "C.ar",
        "CD2": "C.ar",
        "CE1": "C.ar",
        "CE2": "C.ar",
        "CZ": "C.ar",
    },
    "ILE": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG1": "C.3",
        "CG2": "C.3",
        "CD": "C.3",
        "CD1": "C.3",
    },
    "LEU": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.3",
        "CD1": "C.3",
        "CD2": "C.3",
    },
    "VAL": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG1": "C.3",
        "CG2": "C.3",
    },
    "TRP": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.2",
        "CD1": "C.2",
        "CD2": "C.ar",
        "NE1": "N.pl3",
        "CE2": "C.ar",
        "CE3": "C.ar",
        "CZ2": "C.ar",
        "CZ3": "C.ar",
        "CH2": "C.ar",
    },
    "TYR": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.ar",
        "CD1": "C.ar",
        "CD2": "C.ar",
        "CE1": "C.ar",
        "CE2": "C.ar",
        "CZ": "C.ar",
        "OH": "O.3",
    },
    "ALA": {"N": "N.am", "CA": "C.3", "C": "C.2", "O": "O.2", "CB": "C.3"},
    "GLY": {"N": "N.am", "CA": "C.3", "C": "C.2", "O": "O.2"},
    "THR": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "OG1": "O.3",
        "CG2": "C.3",
    },
    "SER": {"N": "N.am", "CA": "C.3", "C": "C.2", "O": "O.2", "CB": "C.3", "OG": "O.3"},
    "GLN": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.3",
        "CD": "C.2",
        "OE1": "O.2",
        "NE2": "N.am",
    },
    "ASN": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.2",
        "OD1": "O.2",
        "ND2": "N.am",
    },
    "GLU": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.3",
        "CD": "C.2",
        "OE1": "O.co2",
        "OE2": "O.co2",
    },
    "ASP": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.2",
        "OD1": "O.co2",
        "OD2": "O.co2",
    },
    "HIS": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.2",
        "ND1": "N.pl3",
        "CD2": "C.2",
        "CE1": "C.2",
        "NE2": "N.2",
    },
    "ARG": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.3",
        "CD": "C.3",
        "NE": "N.pl3",
        "CZ": "C.cat",
        "NH1": "N.pl3",
        "NH2": "N.pl3",
    },
    "LYS": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.3",
        "CD": "C.3",
        "CE": "C.3",
        "NZ": "N.4",
    },
    "PRO": {
        "N": "N.am",
        "CA": "C.3",
        "C": "C.2",
        "O": "O.2",
        "CB": "C.3",
        "CG": "C.3",
        "CD": "C.3",
    },
    "T": {
        "P": "P.3",
        "O1P": "O.co2",
        "O2P": "O.co2",
        "O5*": "O.3",
        "C5*": "C.3",
        "C4*": "C.3",
        "O4*": "O.3",
        "C3*": "C.3",
        "O3*": "O.3",
        "C2*": "C.3",
        "C1*": "C.3",
        "N1": "N.am",
        "C2": "C.2",
        "O2": "O.2",
        "N3": "N.am",
        "C4": "C.2",
        "O4": "O.2",
        "C5": "C.2",
        "C5M": "C.3",
        "C6": "C.2",
    },
    "DT": {
        "P": "P.3",
        "O1P": "O.co2",
        "O2P": "O.co2",
        "O5'": "O.3",
        "C5'": "C.3",
        "C4'": "C.3",
        "O4'": "O.3",
        "C3'": "C.3",
        "O3'": "O.3",
        "C2'": "C.3",
        "C1'": "C.3",
        "N1": "N.am",
        "C2": "C.2",
        "O2": "O.2",
        "N3": "N.am",
        "C4": "C.2",
        "O4": "O.2",
        "C5": "C.2",
        "C7": "C.3",
        "C6": "C.2",
    },
    "A": {
        "P": "P.3",
        "O1P": "O.co2",
        "O2P": "O.co2",
        "O5*": "O.3",
        "C5*": "C.3",
        "C4*": "C.3",
        "O4*": "O.3",
        "C3*": "C.3",
        "O3*": "O.3",
        "C2*": "C.3",
        "O2*": "O.3",
        "C1*": "C.3",
        "N9": "N.pl3",
        "C8": "C.2",
        "N7": "N.2",
        "C5": "C.ar",
        "C6": "C.ar",
        "N6": "N.pl3",
        "N1": "N.ar",
        "C2": "C.ar",
        "N3": "N.ar",
        "C4": "C.ar",
    },
    "DA": {
        "P": "P.3",
        "O1P": "O.co2",
        "O2P": "O.co2",
        "O5'": "O.3",
        "C5'": "C.3",
        "C4'": "C.3",
        "O4'": "O.3",
        "C3'": "C.3",
        "O3'": "O.3",
        "C2'": "C.3",
        "O2'": "O.3",
        "C1'": "C.3",
        "N9": "N.pl3",
        "C8": "C.2",
        "N7": "N.2",
        "C5": "C.ar",
        "C6": "C.ar",
        "N6": "N.pl3",
        "N1": "N.ar",
        "C2": "C.ar",
        "N3": "N.ar",
        "C4": "C.ar",
    },
    "G": {
        "P": "P.3",
        "O1P": "O.co2",
        "O2P": "O.co2",
        "O5*": "O.3",
        "C5*": "C.3",
        "C4*": "C.3",
        "O4*": "O.3",
        "C3*": "C.3",
        "O3*": "O.3",
        "C2*": "C.3",
        "O2*": "O.3",
        "C1*": "C.3",
        "N9": "N.pl3",
        "C8": "C.2",
        "N7": "N.2",
        "C5": "C.2",
        "C6": "C.2",
        "O6": "O.2",
        "N1": "N.am",
        "C2": "C.2",
        "N2": "N.pl3",
        "N3": "N.2",
        "C4": "C.2",
    },
    "DG": {
        "P": "P.3",
        "O1P": "O.co2",
        "O2P": "O.co2",
        "O5'": "O.3",
        "C5'": "C.3",
        "C4'": "C.3",
        "O4'": "O.3",
        "C3'": "C.3",
        "O3'": "O.3",
        "C2'": "C.3",
        "O2'": "O.3",
        "C1'": "C.3",
        "N9": "N.pl3",
        "C8": "C.2",
        "N7": "N.2",
        "C5": "C.2",
        "C6": "C.2",
        "O6": "O.2",
        "N1": "N.am",
        "C2": "C.2",
        "N2": "N.pl3",
        "N3": "N.2",
        "C4": "C.2",
    },
    "C": {
        "P": "P.3",
        "O1P": "O.co2",
        "O2P": "O.co2",
        "O5*": "O.3",
        "C5*": "C.3",
        "C4*": "C.3",
        "O4*": "O.3",
        "C3*": "C.3",
        "O3*": "O.3",
        "C2*": "C.3",
        "O2*": "O.3",
        "C1*": "C.3",
        "N1": "N.am",
        "C2": "C.2",
        "O2": "O.2",
        "N3": "N.2",
        "C4": "C.2",
        "N4": "N.pl3",
        "C5": "C.2",
        "C6": "C.2",
    },
    "DC": {
        "P": "P.3",
        "O1P": "O.co2",
        "O2P": "O.co2",
        "O5'": "O.3",
        "C5'": "C.3",
        "C4'": "C.3",
        "O4'": "O.3",
        "C3'": "C.3",
        "O3'": "O.3",
        "C2'": "C.3",
        "O2'": "O.3",
        "C1'": "C.3",
        "N1": "N.am",
        "C2": "C.2",
        "O2": "O.2",
        "N3": "N.2",
        "C4": "C.2",
        "N4": "N.pl3",
        "C5": "C.2",
        "C6": "C.2",
    },
}

atom_types = [
    "C.3",
    "C.2",
    "C.ar",
    "C.cat",
    "N.4",
    "N.2",
    "N.ar",
    "N.am",
    "N.pl3",
    "O.3",
    "O.2",
    "O.co2",
    "S.3",
    "S.o2",
    "P.3",
    "F",
    "Cl",
    "Br",
    "Met",
]


class DDNAPotential(object):
    """Loads DDNA potentials information"""

    def __init__(self):
        data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
        self.potentials = DDNAPotential._read_potentials(
            os.path.join(data_path, "fort.21_xscore_noH_Met")
        )
        self.map = DDNAPotential._createmap()

    def _read_potentials(data_file_name):
        """Reads DDNA data potentials"""
        potentials = [0.0] * (21 * 20 * 20)
        D1 = 20 * 20
        D2 = 20
        with open(data_file_name) as input_handle:
            for line in input_handle:
                line = line.rstrip(os.linesep)
                fields = line.split()
                if len(fields) == 6:
                    energy = float(fields[2])
                    m = int(fields[3])
                    i = int(fields[4])
                    j = int(fields[5])
                    potentials[m * D1 + i * D2 + j] = potentials[
                        m * D1 + j * D2 + i
                    ] = energy
        return potentials

    def _createmap():
        map = [-1] * 700
        for i in range(1, 50):
            if i < 4:
                map[i] = 1
            if i >= 4 and i < 16:
                map[i] = i - 2
            if i >= 16 and i < 50:
                map[i] = int(i * 0.5) + 6
        return map


class DDNAAdapter(ModelAdapter):
    """Adapts a given Complex to a DockingModel object suitable for this
    DDNA scoring function.
    """

    def _get_docking_model(self, molecule, restraints):
        """Builds a suitable docking model for this scoring function"""
        parsed_restraints = {}
        ddna_objects = []
        atom_index = 0
        for chain in molecule.chains:
            for residue in chain.residues:
                res_id = (
                    f"{chain.cid}.{residue.name}.{residue.number}{residue.insertion}"
                )
                in_restraint = False
                if restraints and res_id in restraints:
                    parsed_restraints[res_id] = []
                    in_restraint = True

                for rec_atom in residue.atoms:
                    try:
                        ddna_atom_type = atom_types.index(
                            atom_map[rec_atom.residue_name][rec_atom.name]
                        )
                    except KeyError:
                        raise NotSupportedInScoringError(
                            f"Atom {rec_atom.name} in residue {rec_atom.residue_name} not supported"
                        )
                    ddna_objects.append(ddna_atom_type)
                    if in_restraint:
                        parsed_restraints[res_id].append(atom_index)
                    atom_index += 1
        try:
            return DockingModel(
                ddna_objects,
                molecule.copy_coordinates(),
                parsed_restraints,
                n_modes=molecule.n_modes.copy(),
            )
        except AttributeError:
            return DockingModel(
                ddna_objects, molecule.copy_coordinates(), parsed_restraints
            )


class DDNA(ScoringFunction):
    """Implements the DDNA potential"""

    def __init__(self, weight=1.0):
        super(DDNA, self).__init__(weight)
        self.potential = DDNAPotential()
        self.cutoff = DEFAULT_CONTACT_RESTRAINTS_CUTOFF

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        energy, interface_receptor, interface_ligand = calculate_ddna(
            receptor,
            receptor_coordinates,
            ligand,
            ligand_coordinates,
            self.potential.potentials,
            self.potential.map,
            interface_cutoff=self.cutoff,
        )

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
DefinedScoringFunction = DDNA
DefinedModelAdapter = DDNAAdapter
