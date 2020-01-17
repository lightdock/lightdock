#!/usr/bin/env python3

"""Prepares a new simulation (init folder) starting from a previous one.

It uses the clustering data if exists the filer cluster.repr in the swarm folder.

This is a very experimental feature and may change in the future without any warning.

"""

import os
import argparse
import numpy as np
import math
from lightdock.constants import DEFAULT_SWARM_FOLDER, CLUSTER_REPRESENTATIVES_FILE, \
    GSO_OUTPUT_FILE, DEFAULT_POSITIONS_FOLDER, DEFAULT_BILD_STARTING_PREFIX, CLUSTERS_CENTERS_FILE, \
    DEFAULT_STARTING_PREFIX
from lightdock.prep.poses import create_file_from_poses
from lightdock.util.logger import LoggingManager
from lightdock.util.analysis import read_lightdock_output
from lightdock.prep.geometry import create_bild_file
from lightdock.pdbutil.PDBIO import create_pdb_from_points


log = LoggingManager.get_logger('lgd_prepare_new_simulation')


def parse_command_line():
    parser = argparse.ArgumentParser(prog='lgd_prepare_new_simulation')
    parser.add_argument("swarm", help="swarm to consider", type=int, metavar="swarm")
    parser.add_argument("steps", help="steps to consider", type=int, metavar="steps")
    parser.add_argument("destination", help="destination folder", metavar="destination")
    parser.add_argument("-nm", "--nm", help="Keep normal modes in starting positions if exist",
                        dest="nm", action="store_true")
    return parser.parse_args()


def read_cluster_representative_data(file_name):
    glowworms = []
    with open(file_name) as input_file:
        for line in input_file:
            line = line.rstrip(os.linesep)
            fields = line.split(':')
            glowworms.append(int(fields[-2]))
    return glowworms


def get_new_positions(pose, n=100, sigma=0.25, tolerance=0.0001, threshold=0.95):
    """Returns random translation and quaternion normalized vectors"""
    positions = []
    positions_generated = 0
    translation = np.array(pose[:3])
    rotation = np.array(pose[3:7])
    while positions_generated < n:
        new_t = [np.random.normal(x, sigma) for x in translation]
        new_q = np.random.randn(4)
        mag2 = sum(n * n for n in new_q)
        if abs(mag2 - 1.0) > tolerance:
            mag = math.sqrt(mag2)
            new_q = np.array(tuple(n / mag for n in new_q))
        dot = np.dot(rotation, new_q)
        if dot > threshold:
            if args.nm:
                positions.append(new_t + list(new_q) + pose[7:])
            else:
                positions.append(new_t + list(new_q))
            positions_generated += 1

    return positions


def create_positions(path, cid, pose):
    positions = get_new_positions(pose)
    positions_file = os.path.join(path, "%s_%d.bild" % (DEFAULT_BILD_STARTING_PREFIX, cid))
    create_bild_file(positions_file, positions)
    pos_file_name = os.path.join(path, '%s_%s.dat' % (DEFAULT_STARTING_PREFIX, cid))
    create_file_from_poses(pos_file_name, positions)


def create_centers_file(path, poses):
    centers_file = os.path.join(path, CLUSTERS_CENTERS_FILE)
    create_pdb_from_points(centers_file, poses)


if __name__ == "__main__":
    try:
        # Parse command line
        args = parse_command_line()

        # Create destination folder
        if os.path.exists(args.destination):
            raise SystemExit("%s already exists. Stop" % args.destination)
        else:
            os.mkdir(args.destination)
            log.info("Folder %s created" % args.destination)

        # Create init folder
        init_folder = os.path.join(args.destination, DEFAULT_POSITIONS_FOLDER)
        os.mkdir(init_folder)

        total_swarms = 0

        swarm_id = args.swarm
        log.info("swarm %d" % swarm_id)
        result_file_name = os.path.join(DEFAULT_SWARM_FOLDER + str(swarm_id), GSO_OUTPUT_FILE % args.steps)
        results = read_lightdock_output(result_file_name)

        cluster_repr_file = os.path.join(DEFAULT_SWARM_FOLDER + str(swarm_id), CLUSTER_REPRESENTATIVES_FILE)
        selected_glowworms = read_cluster_representative_data(cluster_repr_file)

        representative_results = []
        for glowworm_id in selected_glowworms:
            representative_results.append(results[glowworm_id])
            # Create centers and initial positions
            create_positions(init_folder, total_swarms, results[glowworm_id].pose)
            log.info("%8d - Initial positions for glowworm %d created" % (total_swarms, glowworm_id))
            total_swarms += 1
        create_centers_file(init_folder, [result.pose for result in representative_results])
        log.info('Done.')

    except KeyboardInterrupt:
        log.info("Caught interrupt...")
        log.info("bye.")

    except SystemExit as e:
        log.error(str(e))
