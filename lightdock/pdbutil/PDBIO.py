"""Parses Atomic coordinates entries from PDB files"""

import math
from os import linesep
from lightdock.error.lightdock_errors import PDBParsingError, PDBParsingWarning
from lightdock.structure.atom import Atom, HetAtom
from lightdock.structure.residue import Residue
from lightdock.structure.chain import Chain
from lightdock.util.logger import LoggingManager

log = LoggingManager.get_logger("pdb")


def cstrip(string):
    """Remove unwanted symbols from string."""
    return string.strip(" \t\n\r")


def read_atom_line(line, line_type="", atoms_to_ignore=[], residues_to_ignore=[]):
    """Parses a PDB file line starting with 'ATOM' or 'HETATM'"""
    element = cstrip(line[76:78])
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        if math.isnan(x) or math.isnan(y) or math.isnan(z):
            raise Exception()
    except:
        raise PDBParsingError("Wrong coordinates in '%s'" % line)

    try:
        atom_number = int(line[6:11])
    except ValueError:
        raise PDBParsingError("Wrong atom number in '%s'" % line)

    atom_name = cstrip(line[12:16])
    atom_alternative = cstrip(line[16])
    residue_name = cstrip(line[17:21])
    chain_id = cstrip(line[21])
    residue_insertion = line[26]

    try:
        residue_number = int(line[22:26])
    except ValueError:
        raise PDBParsingError("Wrong residue number in '%s'" % line)

    if residue_name in residues_to_ignore:
        raise PDBParsingWarning(
            "Ignored atom %s.%s.%s %s"
            % (chain_id, residue_name, residue_number, atom_name)
        )

    if ("H" in atoms_to_ignore and atom_name[0] == "H") or atom_name in atoms_to_ignore:
        raise PDBParsingWarning(
            "Ignored atom %s.%s.%s %s"
            % (chain_id, residue_name, residue_number, atom_name)
        )

    try:
        occupancy = float(line[54:60])
    except:
        occupancy = 1.0

    try:
        b_factor = float(line[60:66])
    except:
        b_factor = 0.0

    if not line_type:
        line_type = line[0:6].strip()

    if line_type == "ATOM":
        return Atom(
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
        )
    else:
        return HetAtom(
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
        )


def parse_complex_from_file(
    input_file_name, atoms_to_ignore=[], residues_to_ignore=[], verbose=False
):
    """Reads and parses a given input_file_name PDB file.

    TODO: Check if chain have been already created and insert it into the first one
    """
    lines = open(input_file_name).readlines()
    atoms = []
    residues = []
    chains = []
    num_models = 0
    last_chain_id = "#"
    last_residue_name = "#"
    last_residue_number = "#"
    last_residue_insertion = "#"
    current_chain = None
    current_residue = None
    for line in lines:
        # Only first model is going to be read
        if num_models <= 1:
            line_type = line[0:6].strip()
            if line_type == "MODEL":
                num_models += 1
                if num_models > 1:
                    log.warning(
                        "Multiple models found in %s. Only first model will be used."
                        % input_file_name
                    )
            elif line_type == "ATOM" or line_type == "HETATM":
                try:
                    atom = read_atom_line(
                        line, line_type, atoms_to_ignore, residues_to_ignore
                    )
                    atoms.append(atom)
                except PDBParsingWarning as warning:
                    if verbose:
                        print(warning)
                    continue

                if last_chain_id != atom.chain_id:
                    last_chain_id = atom.chain_id
                    current_chain = Chain(last_chain_id)
                    chains.append(current_chain)
                if (
                    last_residue_name != atom.residue_name
                    or last_residue_number != atom.residue_number
                    or last_residue_insertion != atom.residue_insertion
                ):
                    last_residue_name = atom.residue_name
                    last_residue_number = atom.residue_number
                    last_residue_insertion = atom.residue_insertion
                    current_residue = Residue(
                        atom.residue_name, atom.residue_number, atom.residue_insertion
                    )
                    residues.append(current_residue)
                    current_chain.residues.append(current_residue)
                current_residue.atoms.append(atom)

    # Set backbone and side-chain atoms
    for residue in residues:
        residue.set_backbone_and_sidechain()
        try:
            residue.check()
        except Exception as e:
            log.warning("Possible problem: %s" % str(e))

    return atoms, residues, chains


def _format_atom_name(atom_name):
    """Format ATOM name with correct padding"""
    if len(atom_name) == 4:
        return atom_name
    else:
        return " %s" % atom_name


def write_atom_line(atom, atom_coordinates, output):
    """Writes a PDB file format line to output."""
    if atom.__class__.__name__ == "HetAtom":
        atom_type = "HETATM"
    else:
        atom_type = "ATOM  "
    line = "%6s%5d %-4s%-1s%3s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%12s%s" % (
        atom_type,
        atom.number,
        _format_atom_name(atom.name),
        atom.alternative,
        atom.residue_name,
        atom.chain_id,
        atom.residue_number,
        atom.residue_insertion,
        atom_coordinates[atom.index][0],
        atom_coordinates[atom.index][1],
        atom_coordinates[atom.index][2],
        atom.occupancy,
        atom.b_factor,
        atom.element,
        linesep,
    )
    output.write(line)


def write_pdb_to_file(
    molecule, output_file_name, atom_coordinates=None, structure_id=0
):
    """Writes a Complex structure to a file in PDB format."""
    output_file = open(output_file_name, "a")
    for atom in molecule.atoms:
        if atom_coordinates is not None:
            write_atom_line(atom, atom_coordinates, output_file)
        else:
            write_atom_line(atom, molecule.atom_coordinates[structure_id], output_file)
    output_file.close()


def create_pdb_from_points(
    pdb_file_name, points, atom_name="H", res_name="SWR", chain_id="Z", element="H"
):
    """Creates a PDB file which contains an atom_type atom for each point
    in points list.
    """
    with open(pdb_file_name, "w") as output:
        for index, point in enumerate(points):
            line = "ATOM  %5d %-4s %3s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%12s%s" % (
                index + 1,
                _format_atom_name(atom_name),
                res_name,
                chain_id,
                index + 1,
                "",
                point[0],
                point[1],
                point[2],
                1.,
                1.,
                element,
                linesep,
            )

            output.write(line)
