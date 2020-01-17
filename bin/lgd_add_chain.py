#!/usr/bin/env python3

"""Renames or adds chains to PDB files"""

import argparse
import shutil
import os


def valid_file(file_name):
    """Checks if it is a valid file"""
    if not os.path.exists(file_name):
        raise argparse.ArgumentTypeError("The file does not exist")
    return file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="add_chain")
    # PDB structure
    parser.add_argument("pdb_file", help="PDB file name",
                        type=valid_file, metavar="pdb_file")
    # Old chains
    parser.add_argument("old_chains", help="Chains separated by comma",
                        metavar="old_chains")
    # New chains
    parser.add_argument("new_chains", help="Chains separated by comma",
                        metavar="new_chains")

    args = parser.parse_args()

    old_chains = [x.upper() for x in list(set(args.old_chains.split(',')))]
    new_chains = [x.upper() for x in list(set(args.new_chains.split(',')))]

    if len(old_chains) != len(new_chains):
        raise SystemExit("Old and new chains length is different")

    print("Old chains are: %s" % old_chains)
    print("New chains are: %s" % new_chains)

    backup_file = args.pdb_file + '.orig'
    print("Saving a backup to %s" % backup_file)
    shutil.copyfile(args.pdb_file, backup_file)

    translator = {}
    for old_chain, new_chain in zip(old_chains, new_chains):
        translator[old_chain] = new_chain

    with open(backup_file) as pdb_file:
        with open(args.pdb_file, 'w') as output:
            for line in pdb_file:
                try:
                    if line.startswith('ATOM  ') or line.startswith('HETATM'):
                        chain = line[21].strip()
                        if chain in old_chains:
                            line = line[:21] + translator[chain] + line[22:]
                        output.write(line)
                except IndexError:
                    pass
    print("Done.")
