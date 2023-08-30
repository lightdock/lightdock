"""
TOBI atom-based docking potentials-I (ADPs-I), step 1 & 2 as described in:

Tobi, D. Designing coarse grained-and atom based-potentials for protein-protein docking.
BMC Struct Biol 10, 40 (2010). https://doi.org/10.1186/1472-6807-10-40
"""

from pathlib import Path
import numpy as np
import scipy.spatial

from lightdock.structure.model import DockingModel
from lightdock.scoring.functions import ModelAdapter, ScoringFunction
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF


class TOBIA1Potential(object):
    """Loads TOBIA1 potentials information"""

    recognized_residues = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]

    # N, CA, C, O and GLYCA are special
    atom_types = [
        [res + "N" for res in recognized_residues],
        [res + "CA" for res in recognized_residues if res != "GLY"],
        [res + "C" for res in recognized_residues],
        [res + "O" for res in recognized_residues],
        ["GLYCA"],
        [
            "ALACB",
            "ARGCB",
            "ASNCB",
            "ASPCB",
            "CYSCB",
            "GLNCB",
            "GLUCB",
            "HISCB",
            "ILECB",
            "LEUCB",
            "LYSCB",
            "METCB",
            "PHECB",
            "PROCB",
            "PROCG",
            "PROCD",
            "THRCB",
            "TRPCB",
            "TYRCB",
            "VALCB",
        ],
        ["LYSCE", "LYSNZ"],
        ["LYSCD"],
        ["ASPCG", "ASPOD1", "ASPOD2", "GLUCD", "GLUOE1", "GLUOE2"],
        ["ARGCZ", "ARGNH1", "ARGNH2"],
        ["ASNCG", "ASNOD1", "ASNND2", "GLNCD", "GLNOE1", "GLNNE2"],
        ["ARGCD", "ARGNE"],
        ["SERCB", "SEROG", "THROG1", "TYROH"],
        ["HISCG", "HISND1", "HISCD2", "HISCE1", "HISNE2", "TRPNE1"],
        ["TYRCE1", "TYRCE2", "TYRCZ"],
        [
            "ARGCG",
            "GLNCG",
            "GLUCG",
            "ILECG1",
            "LEUCG",
            "LYSCG",
            "METCG",
            "METSD",
            "PHECG",
            "PHECD1",
            "PHECD2",
            "PHECE1",
            "PHECE2",
            "PHECZ",
            "THRCG2",
            "TRPCG",
            "TRPCD1",
            "TRPCD2",
            "TRPCE2",
            "TRPCE3",
            "TRPCZ2",
            "TRPCZ3",
            "TRPCH2",
            "TYRCG",
            "TYRCD1",
            "TYRCD2",
        ],
        ["ILECG2", "ILECD1", "ILECD", "LEUCD1", "LEUCD2", "METCE", "VALCG1", "VALCG2"],
        ["CYSSG"],
    ]

    atom_indice = {at: idx for idx, ats in enumerate(atom_types) for at in ats}

    def __init__(self):
        data_path = Path(__file__).parent.resolve() / "data"
        self.tobi_a_1 = self._read_potentials(data_path / "TOBI_A1_step1")
        self.tobi_a_2 = self._read_potentials(data_path / "TOBI_A1_step2")

    def _read_potentials(self, data_file_name):
        """Reads TOBIA2 data potentials"""
        return np.loadtxt(data_file_name, dtype=float)


class TOBIA1Adapter(ModelAdapter):
    """Adapts a given Complex to a DockingModel object suitable for this
    TOBIA1 scoring function.
    """

    def _get_docking_model(self, molecule, restraints):
        """Builds a suitable docking model for this scoring function"""
        list_of_coordinates = molecule.atom_coordinates
        parsed_restraints = {}
        for chain in molecule.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    atom.index = TOBIA1Potential.atom_indice[
                        f"{atom.residue_name}{atom.name}"
                    ]
                res_id = (
                    f"{chain.cid}.{residue.name}.{residue.number}{residue.insertion}"
                )
                if restraints and res_id in restraints:
                    parsed_restraints[res_id] = []

        return DockingModel(molecule, list_of_coordinates, parsed_restraints)


class TOBIA1(ScoringFunction):
    """Implements TOBIA1 potential"""

    def __init__(self, weight=1.0):
        super(TOBIA1, self).__init__(weight, anm_support=False)
        self.function = self._default
        self.potential = TOBIA1Potential()
        self.cutoff = DEFAULT_CONTACT_RESTRAINTS_CUTOFF

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        return self.function(receptor, receptor_coordinates, ligand, ligand_coordinates)

    def _default(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        energy = 0.0

        dist_matrix_atom = scipy.spatial.distance.cdist(
            receptor_coordinates, ligand_coordinates
        )  # shape:len(rec.atoms)*len(lig.atoms)
        rec_tobi_atom_types = np.array(
            [atom.index for atom in receptor.objects.atoms]
        )  # shape:len(rec.atoms)*1,aij in [0,17]
        lig_tobi_atom_types = np.array(
            [atom.index for atom in ligand.objects.atoms]
        )  # shape:len(lig.atoms)*1,aij in [0,17]
        searched_indice = tuple(
            np.array(
                np.meshgrid(rec_tobi_atom_types, lig_tobi_atom_types, indexing="ij")
            ).reshape((2, -1))
        )  # generate (array[idxx],array[idxy]) for fancy indexing
        tobi_pot1_matrix = self.potential.tobi_a_1[searched_indice].reshape(
            dist_matrix_atom.shape
        )
        tobi_pot2_matrix = self.potential.tobi_a_2[searched_indice].reshape(
            dist_matrix_atom.shape
        )
        mask1 = dist_matrix_atom <= 4.0
        mask2 = np.logical_and(dist_matrix_atom <= 6.0, dist_matrix_atom > 4.0)
        energy = np.sum(tobi_pot1_matrix, where=mask1) + np.sum(
            tobi_pot2_matrix, where=mask2
        )
        rec_coms = np.array(
            [res.get_central_coordinate() for res in receptor.objects.residues]
        )
        lig_coms = np.array(
            [lig.get_central_coordinate() for lig in ligand.objects.residues]
        )
        dist_matrix_residue = scipy.spatial.distance.cdist(rec_coms, lig_coms)
        res_cutoff_indice = dist_matrix_residue[
            dist_matrix_residue <= self.cutoff
        ].nonzero()
        interface_receptor = set(res_cutoff_indice[0])
        interface_ligand = set()
        # nonzero returns ([],) if res_cutoff_indice only contains 0,
        # else returns (idxx,idxy) for 2D-matrix
        if len(res_cutoff_indice) == 2:
            interface_ligand = set(res_cutoff_indice[1])
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
DefinedScoringFunction = TOBIA1
DefinedModelAdapter = TOBIA1Adapter
