#!/usr/bin/env python3

"""Filter LightDock final swarm results depending on the percentage of restraints satisfied"""


import os
import argparse
import shutil
import re
from prody.measure.contacts import Contacts
from prody import parsePDB, confProDy
from lightdock.util.logger import LoggingManager
from lightdock.util.analysis import read_ranking_file


# Disable ProDy output
confProDy(verbosity="info")
filtered_folder = "filtered"

log = LoggingManager.get_logger("lgd_filter_restraints")


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


def get_restraints(restraints_file, include_passive=False, only_passive=False):
    restraints_receptor = set()
    restraints_ligand = set()
    with open(restraints_file) as handle:
        for line in handle:
            line = line.rstrip(os.linesep)
            if line:
                if only_passive:
                    if line.startswith("R") and line.endswith("P"):
                        restraints_receptor.add(line)
                    if line.startswith("L") and line.endswith("P"):
                        restraints_ligand.add(line)
                else:
                    if line.startswith("R"):
                        if line.endswith("P"):
                            if include_passive:
                                restraints_receptor.add(line)
                        else:
                            restraints_receptor.add(line)
                    else:
                        if line.endswith("P"):
                            if include_passive:
                                restraints_ligand.add(line)
                        else:
                            restraints_ligand.add(line)
    return restraints_receptor, restraints_ligand


def parse_command_line():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(prog="lgd_filter_restraints")

    parser.add_argument(
        "ranking_file", help="Path of ranking to be used", metavar="ranking_file"
    )
    parser.add_argument(
        "restraints_file", help="File including restraints", metavar="restraints_file"
    )
    parser.add_argument(
        "receptor_chains",
        help="Chains on the receptor partner",
        metavar="receptor_chains",
    )
    parser.add_argument(
        "ligand_chains", help="Chains on the receptor partner", metavar="ligand_chains"
    )
    parser.add_argument(
        "--cutoff",
        "-cutoff",
        "-c",
        help="Interaction cutoff",
        dest="cutoff",
        type=float,
        default=5.0,
    )
    parser.add_argument(
        "--fnat",
        "-fnat",
        "-f",
        help="Structures with at least this fraction of native contacts",
        dest="fnat",
        type=float,
    )
    parser.add_argument(
        "--rnuc",
        help="Is receptor molecule a nucleic acid?",
        dest="rnuc",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--lnuc",
        help="Is ligand molecule a nucleic acid?",
        dest="lnuc",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--include_passive",
        "-include_passive",
        help="Include passive restraints into filtering",
        dest="include_passive",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--only_passive",
        "-only_passive",
        help="Filter only by passive restraints",
        dest="only_passive",
        action="store_true",
        default=False,
    )

    return parser.parse_args()


if __name__ == "__main__":

    # Parse command line
    args = parse_command_line()

    # Check that --include_passive and --only_passive are not defined simultaneously
    if args.include_passive and args.only_passive:
        log.error("--include_passive and --only_passive are mutually exclusive")
        raise SystemExit

    log.info("Calculating interface at {:3.1f}A".format(args.cutoff))

    # Get ranking
    ranking = read_ranking_file(args.ranking_file)

    # Get all the PDB structures in a given directory
    base_path = os.path.abspath(os.path.dirname(args.ranking_file))
    structures = get_structures(ranking, base_path)

    restraints_receptor, restraints_ligand = get_restraints(
        args.restraints_file, args.include_passive, args.only_passive
    )

    # Total number of restraints to be satisfied
    total = float(len(restraints_receptor) + len(restraints_ligand))

    if os.path.exists(filtered_folder):
        raise SystemExit("Folder {} already exists".format(filtered_folder))
    else:
        os.makedirs(filtered_folder)

    filter_passed = {}
    percentages = {}
    for pdb_file in structures:
        try:
            contacts_receptor = set()
            contacts_ligand = set()

            pdb = pdb_file[0]
            swarm_id = int(re.findall(r"swarm_\d+", pdb)[-1].split("_")[-1])
            glowworm_id = int(re.findall(r"lightdock_\d+", pdb)[-1].split("_")[-1])
            score = float(pdb_file[-1])

            # Read molecule and split by receptor and ligand
            if score > 0.0:
                rec_chains = [
                    chain.strip() for chain in args.receptor_chains.split(",")
                ]
                rec_chains_rst = " ".join(rec_chains)
                lig_chains = [chain.strip() for chain in args.ligand_chains.split(",")]
                lig_chains_rst = " ".join(lig_chains)

                molecule = parsePDB(pdb)

                if args.rnuc:
                    receptor = molecule.select(f"nucleic and (chain {rec_chains_rst})")
                else:
                    receptor = molecule.select(f"protein and (chain {rec_chains_rst})")
                if args.lnuc:
                    ligand = molecule.select(f"nucleic and (chain {lig_chains_rst})")
                else:
                    ligand = molecule.select(f"protein and (chain {lig_chains_rst})")

                # Contacts on receptor side
                protein_contacts = Contacts(receptor)
                contacts = protein_contacts.select(args.cutoff, ligand)
                if contacts:
                    for contact in contacts:
                        contacts_receptor.add(
                            "R {}.{}.{}".format(
                                contact.getChid(),
                                contact.getResname(),
                                contact.getResnum(),
                            )
                        )

                # Contacts on ligand side
                protein_contacts = Contacts(ligand)
                contacts = protein_contacts.select(args.cutoff, receptor)
                if contacts:
                    for contact in contacts:
                        contacts_ligand.add(
                            "L {}.{}.{}".format(
                                contact.getChid(),
                                contact.getResname(),
                                contact.getResnum(),
                            )
                        )

                # Calculate percentage of satisfied restraints
                perc = (
                    len(contacts_receptor & restraints_receptor)
                    + len(contacts_ligand & restraints_ligand)
                ) / total
                percentages[(swarm_id, glowworm_id)] = perc
                if args.fnat:
                    if perc >= args.fnat:
                        shutil.copyfile(
                            pdb,
                            os.path.join(
                                filtered_folder,
                                "swarm_{}_{}.pdb".format(swarm_id, glowworm_id),
                            ),
                        )
                        try:
                            filter_passed[swarm_id].append(glowworm_id)
                        except:
                            filter_passed[swarm_id] = [glowworm_id]
                print("{:40s}  {:5.3f}".format(pdb, perc))

        except Exception as e:
            log.error(
                "Filtering has failed for structure {}. Please see error:".format(pdb)
            )
            log.error(str(e))

    filtered_ranking = os.path.join(filtered_folder, "rank_filtered.list")
    with open(filtered_ranking, "w") as handle:
        for rank in ranking:
            if (
                rank.id_swarm in filter_passed
                and rank.id_glowworm in filter_passed[rank.id_swarm]
            ):
                handle.write(
                    "swarm_{}_{}.pdb   {:5.3f}  {:5.3f}".format(
                        rank.id_swarm,
                        rank.id_glowworm,
                        rank.scoring,
                        percentages[(rank.id_swarm, rank.id_glowworm)],
                    )
                    + os.linesep
                )
