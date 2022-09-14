#!/usr/bin/env python3

import os
import argparse
from prody import parsePDB, confProDy
from math import cos, sin, radians, pi
from numpy import arange
from lightdock.pdbutil.PDBIO import create_pdb_from_points
from lightdock.util.logger import LoggingManager


# Disable ProDy output
confProDy(verbosity="info")

log = LoggingManager.get_logger("lgd_create_membrane")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="lgd_create_membrane")
    parser.add_argument(
        "input_pdb_file", help="Input PDB file", metavar="input_pdb_file"
    )
    parser.add_argument(
        "anchor_residue",
        help="A residue in the format ChainID.ResidueName.ResidueNumber "
        "which indicates the point to start to build the fake-membrane",
        metavar="anchor_residue",
    )
    parser.add_argument(
        "--angular_resolution",
        "-angular_resolution",
        "-ar",
        help="Angular resolution, in degrees",
        dest="angular_resolution",
        type=float,
        default=15.0,
    )
    parser.add_argument(
        "--radius_resolution",
        "-radius_resolution",
        "-rr",
        help="Radius resolution, in Angstroms",
        dest="radius_resolution",
        type=float,
        default=4.0,
    )
    parser.add_argument(
        "--radius_offset",
        "-radius_offset",
        "-ro",
        help="Radius offset, in Angstroms",
        dest="radius_offset",
        type=float,
        default=30.0,
    )
    parser.add_argument(
        "--output_pdb_file",
        "-output_pdb_file",
        "-o",
        help="Output PDB file name",
        dest="output_pdb_file",
        default="membrane.pdb",
    )

    args = parser.parse_args()

    if not os.path.exists(args.input_pdb_file):
        log.error(f"File {args.input_pdb_file} does not exist")
        raise SystemExit
    molecule = parsePDB(args.input_pdb_file)

    try:
        chain, res_name, res_num = args.anchor_residue.split(".")
    except:
        log.error("Wrong residue format")
        raise SystemExit

    ca_residue = molecule.select(
        f"chain {chain} and resname {res_name} and resnum {res_num} name CA"
    )

    # Get anchor residue CA atom provided by user
    if not ca_residue:
        log.error(f"Residue {args.anchor_residue} was not found")
        raise SystemExit
    anchor_coord = ca_residue.getCoords()[0]

    # Get all atoms close in Z axis to this anchor residue CA
    z = anchor_coord[2]
    z_min = z - 0.5
    z_max = z + 0.5

    # Calculate the atom slice using ProDy
    atoms_slice = molecule.select(f"{z_min} <= z <= {z_max}")

    # Keep only coordinates
    slice_coordinates = [atom.getCoords() for atom in atoms_slice]

    # Calculate the center and radius of the slice
    min_x = min_y = 9999999.0
    max_x = max_y = -9999999.0
    mean_x = mean_y = 0.0
    for a in slice_coordinates:
        mean_x += a[0]
        mean_y += a[1]
        if a[0] < min_x:
            min_x = a[0]
        if a[0] > max_x:
            max_x = a[0]
        if a[1] < min_y:
            min_y = a[1]
        if a[1] > max_y:
            max_y = a[1]
    mean_x /= float(len(slice_coordinates))
    mean_y /= float(len(slice_coordinates))

    # Radius
    r_x = max(abs(min_x), abs(max_x)) - mean_x
    r_y = max(abs(min_y), abs(max_y)) - mean_y
    min_r = min(r_x, r_y)

    # Create ring and include a central bead
    center_coord = [mean_x, mean_y, z]
    points = [center_coord]

    # Ring composed of Lightdock membrane CG beads
    density = 0.0
    for r in arange(min_r, min_r + args.radius_offset, args.radius_resolution):
        step = radians(max(0.5, args.angular_resolution - density))
        for theta in arange(0.0, 2 * pi, step):
            c = cos(theta)
            s = sin(theta)
            px = r * c + mean_x
            py = r * s + mean_y

            points.append([px, py, z])

        density += 1.0

    # Save genearted beads as a new PDB file
    create_pdb_from_points(
        args.output_pdb_file, points, atom_name="BJ", res_name="MMB", element="P"
    )

    log.info(f"Membrane PDB file written to [{args.output_pdb_file}]")
    log.info(f"- Angular resolution: {args.angular_resolution} degrees")
    log.info(f"- Radius resolution: {args.radius_resolution} A")
    log.info(f"- Radius offset: {args.radius_offset} A")
