#!/usr/bin/env python3

"""Copy clustered structures to new folder for analysis"""

import os
import argparse
import shutil
import re
from lightdock.util.logger import LoggingManager
from lightdock.util.analysis import read_ranking_file


clustered_folder = "clustered"

log = LoggingManager.get_logger("lgd_copy_structures")


def get_structures(ranking, base_path="."):
    structures = []
    for rank in ranking:
        swarm_id = rank.id_swarm
        glowworm_id = rank.id_glowworm
        score = rank.scoring
        structures.append(
            [
                os.path.join(
                    base_path,
                    "swarm_{}".format(swarm_id),
                    "lightdock_{}.pdb".format(glowworm_id),
                ),
                score,
            ]
        )
    return structures


def parse_command_line():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(prog="lgd_copy_structures")

    parser.add_argument(
        "ranking_file", help="Path of ranking to be used", metavar="ranking_file"
    )

    return parser.parse_args()


if __name__ == "__main__":

    # Parse command line
    args = parse_command_line()

    # Get ranking
    ranking = read_ranking_file(args.ranking_file)

    # Get all the PDB structures in a given directory
    base_path = os.path.abspath(os.path.dirname(args.ranking_file))
    structures = get_structures(ranking, base_path)

    if os.path.exists(clustered_folder):
        raise SystemExit("Folder {} already exists".format(clustered_folder))
    else:
        os.makedirs(clustered_folder)

    for pdb_file in structures:
        pdb = pdb_file[0]
        swarm_id = int(re.findall(r"swarm_\d+", pdb)[0].split("_")[-1])
        glowworm_id = int(re.findall(r"lightdock_\d+", pdb)[0].split("_")[-1])
        shutil.copyfile(
            pdb,
            os.path.join(
                clustered_folder, "swarm_{}_{}.pdb".format(swarm_id, glowworm_id)
            ),
        )

    clustered_ranking = os.path.join(clustered_folder, "rank_clustered.list")
    with open(clustered_ranking, "w") as handle:
        for rank in ranking:
            handle.write(
                "swarm_{}_{}.pdb   {:5.3f}  ".format(
                    rank.id_swarm, rank.id_glowworm, rank.scoring
                )
                + os.linesep
            )
