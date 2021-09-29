"""DFIRE potentials scoring functions

S. Liu, C. Zhang, H. Zhou, and Y. Zhou, A physical reference state unifies the structure-derived
potential of mean force for protein folding and binding. Proteins 56, 93-101 (2004)
"""

import os

from lightdock.structure.model import DockingModel
from lightdock.scoring.functions import ModelAdapter, ScoringFunction
from lightdock.scoring.dfire.cython.cdfire import calculate_dfire
from lightdock.util.logger import LoggingManager
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF


log = LoggingManager.get_logger("dfire")


class DFIREPotential(object):
    """Loads DFIRE potentials information"""

    atoms_in_residues = {
        "ALA": ["N", "CA", "C", "O", "CB"],
        "CYS": ["N", "CA", "C", "O", "CB", "SG"],
        "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
        "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
        "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "GLY": ["N", "CA", "C", "O"],
        "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
        "ILE": [
            "N",
            "CA",
            "C",
            "O",
            "CB",
            "CG1",
            "CG2",
            "CD1",
        ],  # or maybe 'CD' if CHARMM
        "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
        "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
        "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
        "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
        "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
        "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
        "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
        "SER": ["N", "CA", "C", "O", "CB", "OG"],
        "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
        "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
        "TRP": [
            "N",
            "CA",
            "C",
            "O",
            "CB",
            "CG",
            "CD1",
            "CD2",
            "CE2",
            "NE1",
            "CE3",
            "CZ3",
            "CH2",
            "CZ2",
        ],
        "TYR": [
            "N",
            "CA",
            "C",
            "O",
            "CB",
            "CG",
            "CD1",
            "CD2",
            "CE1",
            "CE2",
            "CZ",
            "OH",
        ],
    }

    # Recognized residues in order in params file
    RES_3 = [
        "ALA",
        "CYS",
        "ASP",
        "GLU",
        "PHE",
        "GLY",
        "HIS",
        "ILE",
        "LYS",
        "LEU",
        "MET",
        "ASN",
        "PRO",
        "GLN",
        "ARG",
        "SER",
        "THR",
        "VAL",
        "TRP",
        "TYR",
    ]

    # DFIRE only uses 20 distance bins
    dfire_dist_to_bins = [
        1,
        1,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        14,
        15,
        15,
        16,
        16,
        17,
        17,
        18,
        18,
        19,
        19,
        20,
        20,
        21,
        21,
        22,
        22,
        23,
        23,
        24,
        24,
        25,
        25,
        26,
        26,
        27,
        27,
        28,
        28,
        29,
        29,
        30,
        30,
        31,
    ]

    # DFIRE has 167 atom types
    dfire_atom_types = [
        "CYSN",
        "CYSCA",
        "CYSC",
        "CYSO",
        "CYSCB",
        "CYSSG",
        "METN",
        "METCA",
        "METC",
        "METO",
        "METCB",
        "METCG",
        "METSD",
        "METCE",
        "PHEN",
        "PHECA",
        "PHEC",
        "PHEO",
        "PHECB",
        "PHECG",
        "PHECD1",
        "PHECD2",
        "PHECE1",
        "PHECE2",
        "PHECZ",
        "ILEN",
        "ILECA",
        "ILEC",
        "ILEO",
        "ILECB",
        "ILECG1",
        "ILECG2",
        "ILECD1",
        "LEUN",
        "LEUCA",
        "LEUC",
        "LEUO",
        "LEUCB",
        "LEUCG",
        "LEUCD1",
        "LEUCD2",
        "VALN",
        "VALCA",
        "VALC",
        "VALO",
        "VALCB",
        "VALCG1",
        "VALCG2",
        "TRPN",
        "TRPCA",
        "TRPC",
        "TRPO",
        "TRPCB",
        "TRPCG",
        "TRPCD1",
        "TRPCD2",
        "TRPNE1",
        "TRPCE2",
        "TRPCE3",
        "TRPCZ2",
        "TRPCZ3",
        "TRPCH2",
        "TYRN",
        "TYRCA",
        "TYRC",
        "TYRO",
        "TYRCB",
        "TYRCG",
        "TYRCD1",
        "TYRCD2",
        "TYRCE1",
        "TYRCE2",
        "TYRCZ",
        "TYROH",
        "ALAN",
        "ALACA",
        "ALAC",
        "ALAO",
        "ALACB",
        "GLYN",
        "GLYCA",
        "GLYC",
        "GLYO",
        "THRN",
        "THRCA",
        "THRC",
        "THRO",
        "THRCB",
        "THROG1",
        "THRCG2",
        "SERN",
        "SERCA",
        "SERC",
        "SERO",
        "SERCB",
        "SEROG",
        "GLNN",
        "GLNCA",
        "GLNC",
        "GLNO",
        "GLNCB",
        "GLNCG",
        "GLNCD",
        "GLNOE1",
        "GLNNE2",
        "ASNN",
        "ASNCA",
        "ASNC",
        "ASNO",
        "ASNCB",
        "ASNCG",
        "ASNOD1",
        "ASNND2",
        "GLUN",
        "GLUCA",
        "GLUC",
        "GLUO",
        "GLUCB",
        "GLUCG",
        "GLUCD",
        "GLUOE1",
        "GLUOE2",
        "ASPN",
        "ASPCA",
        "ASPC",
        "ASPO",
        "ASPCB",
        "ASPCG",
        "ASPOD1",
        "ASPOD2",
        "HISN",
        "HISCA",
        "HISC",
        "HISO",
        "HISCB",
        "HISCG",
        "HISND1",
        "HISCD2",
        "HISCE1",
        "HISNE2",
        "ARGN",
        "ARGCA",
        "ARGC",
        "ARGO",
        "ARGCB",
        "ARGCG",
        "ARGCD",
        "ARGNE",
        "ARGCZ",
        "ARGNH1",
        "ARGNH2",
        "LYSN",
        "LYSCA",
        "LYSC",
        "LYSO",
        "LYSCB",
        "LYSCG",
        "LYSCD",
        "LYSCE",
        "LYSNZ",
        "PRON",
        "PROCA",
        "PROC",
        "PROO",
        "PROCB",
        "PROCG",
        "PROCD",
    ]

    def __init__(self):
        data_path = os.path.dirname(os.path.realpath(__file__)) + "/data/"

        self.r3_to_numerical = {}
        for x in range(len(DFIREPotential.RES_3)):
            self.r3_to_numerical[DFIREPotential.RES_3[x]] = x

        self.atomnumber = {}
        for x in range(len(DFIREPotential.RES_3)):
            for y in range(
                len(DFIREPotential.atoms_in_residues[DFIREPotential.RES_3[x]])
            ):
                name = "%s%s" % (
                    DFIREPotential.RES_3[x],
                    DFIREPotential.atoms_in_residues[DFIREPotential.RES_3[x]][y],
                )
                self.atomnumber[name] = y

        self.dfire_energy = self._read_potentials(data_path + "DCparams")

    def _read_potentials(self, data_file_name):
        """Reads DFIRE data potentials"""
        dfire_energy = []
        for x in range(len(DFIREPotential.RES_3)):
            dfire_energy.append([])
            for y in range(
                len(DFIREPotential.atoms_in_residues[DFIREPotential.RES_3[x]])
            ):
                dfire_energy[x].append([])
                for a in range(len(DFIREPotential.RES_3)):
                    dfire_energy[x][y].append([])
                    for b in range(
                        len(DFIREPotential.atoms_in_residues[DFIREPotential.RES_3[a]])
                    ):
                        dfire_energy[x][y][a].append([])
                        for _ in range(20):
                            dfire_energy[x][y][a][b].append(99999.9)

        infile = open(data_file_name).readlines()
        count = 0
        for x in range(167):
            residuea = DFIREPotential.dfire_atom_types[x][:3]
            rnuma = self.r3_to_numerical[residuea]
            anuma = self.atomnumber[DFIREPotential.dfire_atom_types[x]]
            for y in range(167):
                residueb = DFIREPotential.dfire_atom_types[y][:3]
                rnumb = self.r3_to_numerical[residueb]
                anumb = self.atomnumber[DFIREPotential.dfire_atom_types[y]]
                for z in range(20):
                    dfire_energy[rnuma][anuma][rnumb][anumb][z] = float(
                        infile[count].strip()
                    )
                    count += 1

        return dfire_energy


class DFIREObject(object):
    def __init__(self, residue_index, dfire_residue_index, atom_index):
        self.residue_index = residue_index
        self.dfire_residue_index = dfire_residue_index
        self.atom_index = atom_index


class DFIREAdapter(ModelAdapter):
    """Adapts a given Complex to a DockingModel object suitable for this
    DFIRE scoring function.
    """

    def _get_docking_model(self, molecule, restraints):
        """Builds a suitable docking model for this scoring function"""
        r3_to_numerical = {}
        for x in range(len(DFIREPotential.RES_3)):
            r3_to_numerical[DFIREPotential.RES_3[x]] = x

        atomnumber = {}
        for x in range(len(DFIREPotential.RES_3)):
            for y in range(
                len(DFIREPotential.atoms_in_residues[DFIREPotential.RES_3[x]])
            ):
                name = "%s%s" % (
                    DFIREPotential.RES_3[x],
                    DFIREPotential.atoms_in_residues[DFIREPotential.RES_3[x]][y],
                )
                atomnumber[name] = y

        parsed_restraints = {}
        dfire_objects = []
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
                    rec_atom_type = rec_atom.residue_name + rec_atom.name
                    rnuma = r3_to_numerical[rec_atom.residue_name]
                    anuma = atomnumber[rec_atom_type]
                    dfire_objects.append(DFIREObject(residue.index, rnuma, anuma))
                    if in_restraint:
                        parsed_restraints[res_id].append(atom_index)
                    atom_index += 1
        try:
            return DockingModel(
                dfire_objects,
                molecule.copy_coordinates(),
                parsed_restraints,
                n_modes=molecule.n_modes.copy(),
            )
        except AttributeError:
            return DockingModel(
                dfire_objects, molecule.copy_coordinates(), parsed_restraints
            )


class DFIRE(ScoringFunction):
    """Implements DFIRE potential"""

    def __init__(self, weight=1.0):
        super(DFIRE, self).__init__(weight)
        self.potential = DFIREPotential()
        self.cutoff = DEFAULT_CONTACT_RESTRAINTS_CUTOFF

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        energy, interface_receptor, interface_ligand = calculate_dfire(
            receptor,
            receptor_coordinates,
            ligand,
            ligand_coordinates,
            self.potential.dfire_dist_to_bins,
            self.potential.dfire_energy,
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
DefinedScoringFunction = DFIRE
DefinedModelAdapter = DFIREAdapter
