#!/usr/bin/env python3

"""Calculates the reference points of the simulation"""

import os
import argparse
from lightdock.constants import (
    DEFAULT_LIST_EXTENSION,
    DEFAULT_REFERENCE_POINTS_EXTENSION,
)
from lightdock.error.lightdock_errors import MinimumVolumeEllipsoidError
from lightdock.mathutil.ellipsoid import MinimumVolumeEllipsoid
from lightdock.pdbutil.PDBIO import parse_complex_from_file, create_pdb_from_points
from lightdock.structure.complex import Complex
from lightdock.util.logger import LoggingManager


script_name = "reference_points"
log = LoggingManager.get_logger(script_name)


def parse_command_line():
    parser = argparse.ArgumentParser(prog=script_name)
    parser.add_argument(
        "structure", help="structure to calculate reference points", metavar="structure"
    )
    parser.add_argument(
        "--noxt",
        help="Remove OXT atoms",
        dest="noxt",
        action="store_true",
        default=False,
    )
    return parser.parse_args()


def get_pdb_files(input_file_name):
    """Get a list of the PDB files in the input_file_name"""
    structure_file_names = []
    with open(input_file_name) as input_file:
        for line in input_file:
            structure_file_name = line.rstrip(os.linesep)
            if os.path.exists(structure_file_name):
                structure_file_names.append(structure_file_name)
    return structure_file_names


def get_point_respresentation(point):
    return "%8.5f %8.5f %8.5f" % (point[0], point[1], point[2])


if __name__ == "__main__":
    try:
        # Parse command line
        args = parse_command_line()

        atoms_to_ignore = []
        if args.noxt:
            atoms_to_ignore.append("OXT")

        structures = []
        file_names = []
        file_name, file_extension = os.path.splitext(args.structure)
        if file_extension == DEFAULT_LIST_EXTENSION:
            file_names.extend(get_pdb_files(args.structure))
        else:
            file_names.append(args.structure)
        for file_name in file_names:
            log.info("Reading %s PDB file..." % file_name)
            atoms, residues, chains = parse_complex_from_file(
                file_name, atoms_to_ignore
            )
            structures.append(
                {
                    "atoms": atoms,
                    "residues": residues,
                    "chains": chains,
                    "file_name": file_name,
                }
            )
            log.info("%s atoms, %s residues read." % (len(atoms), len(residues)))

        molecule = Complex.from_structures(structures)
        try:
            ellipsoid = MinimumVolumeEllipsoid(molecule.atom_coordinates[0].coordinates)
        except MinimumVolumeEllipsoidError as e:
            log.error(
                "Impossible to calculate minimum volume ellipsoid. Reason: %s" % str(e)
            )
            raise SystemExit("%s finished with error" % script_name)

        output_file_name = (
            molecule.structure_file_names[0] + DEFAULT_REFERENCE_POINTS_EXTENSION
        )
        with open(output_file_name, "w") as output:
            for point in ellipsoid.poles:
                output.write(get_point_respresentation(point) + os.linesep)
            output.write(get_point_respresentation(ellipsoid.center) + os.linesep)
        log.info("Points written to %s" % output_file_name)

        points = [point for point in ellipsoid.poles]
        points.append(ellipsoid.center)
        pdb_file_name = output_file_name + ".pdb"
        create_pdb_from_points(pdb_file_name, points)
        log.info("Points written to %s in pdb format" % pdb_file_name)
        log.info("Done.")

    except KeyboardInterrupt:
        log.info("Caught interrupt...")
        log.info("bye.")
