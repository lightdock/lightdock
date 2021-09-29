#!/usr/bin/env python3

"""Calculates the scoring function value for a pair of receptor and ligand PDB structures"""

import argparse
import importlib
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.util.logger import LoggingManager


log = LoggingManager.get_logger("calculate_scoring")


def parse_command_line():
    parser = argparse.ArgumentParser(prog="calculate_scoring")
    parser.add_argument("scoring_function", help="scoring function")
    parser.add_argument("receptor", help="PDB receptor")
    parser.add_argument("ligand", help="PDB ligand")
    script_args = parser.parse_args()
    return script_args


if __name__ == "__main__":
    args = parse_command_line()

    try:
        scoring_function_module = "lightdock.scoring.%s.driver" % args.scoring_function
        module = importlib.import_module(scoring_function_module)
    except ImportError:
        raise SystemExit("Scoring function not found or not available")

    atoms, residues, chains = parse_complex_from_file(args.receptor)
    receptor = Complex(chains, atoms, structure_file_name=args.receptor)
    atoms, residues, chains = parse_complex_from_file(args.ligand)
    ligand = Complex(chains, atoms, structure_file_name=args.ligand)

    CurrentScoringFunction = getattr(module, "DefinedScoringFunction")
    CurrentModelAdapter = getattr(module, "DefinedModelAdapter")

    adapter = CurrentModelAdapter(receptor, ligand)
    scoring_function = CurrentScoringFunction()

    energy = scoring_function(
        adapter.receptor_model,
        adapter.receptor_model.coordinates[0],
        adapter.ligand_model,
        adapter.ligand_model.coordinates[0],
    )
    print(energy)
