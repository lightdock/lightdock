#!/usr/bin/env python3

"""Calculates the ranking files depending of different metrics"""

import os
import argparse
from lightdock.constants import (
    DEFAULT_SWARM_FOLDER,
    GSO_OUTPUT_FILE,
    EVALUATION_FILE,
    SCORING_FILE,
    LIGHTDOCK_PDB_FILE,
    CLUSTER_REPRESENTATIVES_FILE,
)
from lightdock.util.logger import LoggingManager
from lightdock.util.analysis import (
    read_rmsd_and_contacts_data,
    read_lightdock_output,
    write_ranking_to_file,
    read_cluster_representatives_file,
)


log = LoggingManager.get_logger("lgd_rank")


def parse_command_line():
    parser = argparse.ArgumentParser(prog="lgd_rank")
    parser.add_argument(
        "num_swarms",
        help="number of swarms to consider",
        type=int,
        metavar="num_swarms",
    )
    parser.add_argument("steps", help="steps to consider", type=int, metavar="steps")
    parser.add_argument(
        "-c",
        "--clashes_cutoff",
        help="clashes cutoff",
        dest="clashes_cutoff",
        type=float,
    )
    parser.add_argument(
        "-f",
        "--file_name",
        help="lightdock output file to consider",
        dest="result_file",
    )
    parser.add_argument(
        "--ignore_clusters",
        help="Ignore cluster information",
        dest="ignore_clusters",
        action="store_true",
    )
    return parser.parse_args()


if __name__ == "__main__":
    try:
        # Parse command line
        args = parse_command_line()

        solutions = []
        contacts = []
        rmsds = []
        if os.path.isfile(EVALUATION_FILE):
            contacts, rmsds = read_rmsd_and_contacts_data(EVALUATION_FILE)

        num_swarms_found = 0
        for swarm_id in range(args.num_swarms):
            if args.result_file:
                result_file_name = os.path.join(
                    DEFAULT_SWARM_FOLDER + str(swarm_id), args.result_file
                )
            else:
                result_file_name = os.path.join(
                    DEFAULT_SWARM_FOLDER + str(swarm_id), (GSO_OUTPUT_FILE % args.steps)
                )

            cluster_representatives_file = os.path.join(
                DEFAULT_SWARM_FOLDER + str(swarm_id), CLUSTER_REPRESENTATIVES_FILE
            )
            clusters = []
            if (
                os.path.isfile(cluster_representatives_file)
                and not args.ignore_clusters
            ):
                clusters = read_cluster_representatives_file(
                    cluster_representatives_file
                )

            scoring_file_name = os.path.join(
                DEFAULT_SWARM_FOLDER + str(swarm_id), SCORING_FILE
            )
            try:
                results = read_lightdock_output(result_file_name)
                num_swarms_found += 1
                for result in results:
                    result.id_swarm = swarm_id
                    result.pdb_file = LIGHTDOCK_PDB_FILE % result.id_glowworm
                    try:
                        result.rmsd = rmsds[result.id_swarm][result.id_glowworm]
                        result.contacts = contacts[result.id_swarm][result.id_glowworm]
                    except Exception:
                        pass
                    if len(clusters):
                        # Clusters read
                        if result.id_glowworm in clusters:
                            solutions.append(result)
                    else:
                        # Default without clustering
                        solutions.append(result)
            except IOError:
                log.warning("Results %s not found, ignoring." % result_file_name)

        write_ranking_to_file(solutions, args.clashes_cutoff)
        write_ranking_to_file(solutions, args.clashes_cutoff, order_by="luciferin")
        write_ranking_to_file(solutions, args.clashes_cutoff, order_by="rmsd")
        write_ranking_to_file(solutions, args.clashes_cutoff, order_by="scoring")

        log.info("Number of swarms: %d" % args.num_swarms)
        log.info("Number of steps: %d" % args.steps)
        if args.clashes_cutoff:
            log.info("Clashes cutoff: %5.3f" % args.clashes_cutoff)
        if args.result_file:
            log.info("Output files: %s" % args.result_file)
        else:
            log.info("Output files: %s" % (GSO_OUTPUT_FILE % args.steps))
        log.info("Number of swarms processed: %d" % num_swarms_found)
        log.info("Done.")

    except KeyboardInterrupt:
        log.info("Caught interrupt...")
        log.info("bye.")
