#!/usr/bin/env python3

"""Creates a PDB with atom points representing the position for each of the glowworms of a swarm"""

import argparse
import os
from lightdock.error.lightdock_errors import LightDockError
from lightdock.pdbutil.PDBIO import create_pdb_from_points
from lightdock.util.logger import LoggingManager
from lightdock.constants import DEFAULT_LIST_EXTENSION, DEFAULT_LIGHTDOCK_PREFIX


log = LoggingManager.get_logger('generate_glowworm_positions')


def valid_file(file_name):
    """Checks if it is a valid file"""
    if not os.path.exists(file_name):
        raise argparse.ArgumentTypeError("The file does not exist")
    return file_name


def get_lightdock_structures(input_file):
    """Get a list of the PDB files in the input_file"""
    input_file_name, input_file_extension = os.path.splitext(input_file)
    file_names = []
    if input_file_extension == DEFAULT_LIST_EXTENSION:
        with open(input_file) as input_lines:
            for line in input_lines:
                file_name = line.rstrip(os.linesep)
                lightdock_structure = os.path.join(os.path.dirname(file_name),
                                                   DEFAULT_LIGHTDOCK_PREFIX % os.path.basename(file_name))
                if os.path.exists(lightdock_structure):
                    file_names.append(lightdock_structure)
    else:
        file_name = input_file
        lightdock_structure = os.path.join(os.path.dirname(file_name),
                                           DEFAULT_LIGHTDOCK_PREFIX % os.path.basename(file_name))
        if os.path.exists(lightdock_structure):
            file_names.append(lightdock_structure)
        else:
            raise LightDockError('Structure file %s not found' % lightdock_structure)
    return file_names


def parse_output_file(lightdock_output):
    glowworm_translations = []

    data_file = open(lightdock_output)
    lines = data_file.readlines()
    data_file.close()

    counter = 0
    for line in lines:
        if line[0] == '(':
            counter += 1
            last = line.index(')')
            coord = line[1:last].split(',')
            glowworm_translations.append([float(coord[0]), float(coord[1]), float(coord[2])])
    log.info("Read %s coordinate lines" % counter)
    return glowworm_translations


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="generate_glowworm_positions")
    # Lightdock output file
    parser.add_argument("lightdock_output", help="lightdock output file",
                        type=valid_file, metavar="lightdock_output")
    args = parser.parse_args()

    # Output file
    translations = parse_output_file(args.lightdock_output)

    # Destination path is the same as the lightdock output
    destination_path = os.path.dirname(args.lightdock_output)
    pdb_file_name = os.path.splitext(args.lightdock_output)[0] + '.pdb'

    create_pdb_from_points(os.path.join(destination_path, pdb_file_name), translations)
    log.info("PDB %s file created." % os.path.join(destination_path, pdb_file_name))
