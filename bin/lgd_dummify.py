#!/usr/bin/env python3

"""Converts any HETATM molecule to dummy beads"""

import os
import argparse
from prody import parsePDB, writePDB, confProDy
from lightdock.util.logger import LoggingManager

# Disable ProDy output
confProDy(verbosity="info")

log = LoggingManager.get_logger("lgd_dummify")


def parse_command_line():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(prog="lgd_dummify")

    parser.add_argument(
        "input_pdb", help="Input PDB file name", metavar="input_pdb"
    )
    parser.add_argument(
        "output_pdb", help="Output PDB file name", metavar="output_pdb"
    )

    return parser.parse_args()


if __name__ == "__main__":

    # Parse command line
    args = parse_command_line()

    molecule = parsePDB(args.input_pdb)

    hetero_selection = "hetero and not water and not lipid and not ion and not heme"
    hetero = molecule.select(hetero_selection)
    not_hetero = molecule.select(f"not ({hetero_selection})")

    if hetero:
        carbon = hetero.select(f"carbon")
        if carbon:
            log.info(f"Number of dummy atom candidates: {len(carbon)}")
            for atom in carbon:
                atom.setName("DU")
                atom.setResname("MMY")
                atom.setElement("P")

            writePDB(args.output_pdb, not_hetero+carbon)

            log.info(f"Original structure with new dummy atoms written to {args.output_pdb}")
        else:
            log.warning("No dummy atom candidates found from carbon selection, stopping")
            raise SystemExit
    else:
        log.warning("No dummy atom candidates found, stopping")
        raise SystemExit
