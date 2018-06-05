#!/usr/bin/env python

"""Cluster LightDock final swarm results using an hierarchical algorithm"""

import os
import argparse
import math
import numpy as np
import scipy.cluster.hierarchy as hier
import Bio.PDB
from lightdock.constants import CLUSTER_ANALYSIS_FILE, DEFAULT_SWARM_FOLDER, DEFAULT_RMSD_EXTENSION, \
    NUMPY_FILE_SAVE_EXTENSION, EVALUATION_FILE, SCORING_FILE, GSO_OUTPUT_FILE, LIGHTDOCK_PDB_FILE, \
    CLUSTER_DEFAULT_NAME, CLUSTER_REPRESENTATIVES_FILE
from lightdock.util.logger import LoggingManager
from lightdock.util.analysis import read_rmsd_and_contacts_data, read_lightdock_output


log = LoggingManager.get_logger('lgd_cluster_hierarchical')


POPULATION_THRESHOLD = 10


def parse_command_line():
    """
    Parses command line arguments
    """
    parser = argparse.ArgumentParser(prog='cluster_poses')

    parser.add_argument("swarm_id", help="swarm to consider for clustering", type=int, metavar="swarm_id")
    parser.add_argument("steps", help="steps to consider", type=int, metavar="steps")
    parser.add_argument("-f", "--file_name", help="lightdock output file to consider", dest="result_file")
    parser.add_argument("-p", "--ponderated", help="Structures selection takes into account cluster population",
                        dest="ponderated", action="store_true")
    parser.add_argument("-r", "--rmsd", help="RMSD matrix in binary numpy format", dest="rmsd_file", action="store_true")
    parser.add_argument("-t", "--threshold", help="Cluster population threshold", dest="population_threshold",
                        type=int, default=POPULATION_THRESHOLD)
    return parser.parse_args()


def calculate_inter_rmsd(swarm_id):
    N = len(solutions)
    distances = np.zeros((N, N))
    indexes = np.triu_indices(N)
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    super_imposer = Bio.PDB.Superimposer()

    ca_atoms = [[] for _ in xrange(N)]
    for i in range(N):
        log.info('Reading structure %d' % i)
        structure_file = os.path.join(DEFAULT_SWARM_FOLDER + str(swarm_id),
                                      LIGHTDOCK_PDB_FILE % i)
        structure = pdb_parser.get_structure("reference", structure_file)
        model = structure[0]
        for chain in model:
            for residue in chain:
                ca_atoms[i].append(residue['CA'])
    log.info('Calculating RMSD Matrix...')
    for i, j in zip(indexes[0], indexes[1]):
        if i != j:
            super_imposer.set_atoms(ca_atoms[i], ca_atoms[j])
            distances[i][j] = super_imposer.rms
            distances[j][i] = distances[i][j]
    numpy_file_name = os.path.join(DEFAULT_SWARM_FOLDER + str(swarm_id),
                                   CLUSTER_DEFAULT_NAME + DEFAULT_RMSD_EXTENSION)
    np.save(numpy_file_name, distances)
    log.info('Done.')
    return distances


def stats(data):
    return round(np.mean(data), 3), round(np.median(data), 3), round(np.std(data), 3)


if __name__ == "__main__":
    try:
        # Parse command line
        args = parse_command_line()

        if not os.path.isfile(EVALUATION_FILE):
            raise SystemExit("%s file not found" % EVALUATION_FILE)

        # Get contacts and rmsds
        contacts, rmsds = read_rmsd_and_contacts_data(EVALUATION_FILE)

        swarm_id = args.swarm_id
        log.info("cluster %d" % swarm_id)
        solutions = []
        if args.result_file:
            result_file_name = os.path.join(DEFAULT_SWARM_FOLDER + str(swarm_id), args.result_file)
        else:
            result_file_name = os.path.join(DEFAULT_SWARM_FOLDER + str(swarm_id), GSO_OUTPUT_FILE % args.steps)
        scoring_file_name = os.path.join(DEFAULT_SWARM_FOLDER + str(swarm_id), SCORING_FILE)
        results = read_lightdock_output(result_file_name)
        for result in results:
            result.id_cluster = swarm_id
            try:
                result.rmsd = rmsds[result.id_cluster][result.id_glowworm]
                result.contacts = contacts[result.id_cluster][result.id_glowworm]
                result.pdb_file = LIGHTDOCK_PDB_FILE % result.id_glowworm
            except:
                pass
            solutions.append(result)

        if args.rmsd_file:
            log.info('Previous RMSD matrix found. Loading...')
            rmsd_matrix_file = os.path.join(DEFAULT_SWARM_FOLDER + str(swarm_id),
                                            CLUSTER_DEFAULT_NAME +
                                            DEFAULT_RMSD_EXTENSION + NUMPY_FILE_SAVE_EXTENSION)
            distances = np.load(rmsd_matrix_file)
        else:
            log.info('Calculating RMSD distances...')
            distances = calculate_inter_rmsd(swarm_id)
        log.info('Done.')

        # Calculate clusters
        clusters = hier.fclusterdata(distances, distances.max(), criterion='maxclust',
                                     metric='euclidean', depth=2, method='complete')

        # Save data
        data_file_name = os.path.join(DEFAULT_SWARM_FOLDER + str(swarm_id), CLUSTER_ANALYSIS_FILE)
        with open(data_file_name, 'w') as output:
            output.write("Clusters found: %d" % max(clusters) + os.linesep)
            output.write(os.linesep)
            output.write(str(clusters))
            output.write(os.linesep)
            output.write(os.linesep)
            solutions_clustered = {}
            for id_solution, solution in enumerate(solutions):
                id_cluster = clusters[id_solution]
                try:
                    solutions_clustered[id_cluster].append(solution)
                except KeyError:
                    solutions_clustered[id_cluster] = [solution]

            solutions_reclustered = {}
            for k,v in zip(solutions_clustered.keys(), solutions_clustered.values()):
                if len(v) >= args.population_threshold:
                    solutions_reclustered[k] = v

            # Analysis
            for id_cluster in sorted(solutions_reclustered, key=lambda k: len(solutions_reclustered[k]), reverse=True):
                cluster_solutions = solutions_reclustered[id_cluster]
                output.write("[%d] %d" % (id_cluster, len(cluster_solutions)) + os.linesep)
                structures = [solution.pdb_file for solution in cluster_solutions]
                output.write("%s%s" % (' '.join(structures), os.linesep))
                output.write("RMSD: %s%s" % (stats([solution.rmsd for solution in cluster_solutions]), os.linesep))
                output.write("Scoring: %s%s" % (stats([solution.scoring for solution in cluster_solutions]), os.linesep))
                output.write("Neighbors: %s%s" % (stats([solution.num_neighbors for solution in cluster_solutions]), os.linesep))
                output.write("Clashes: %s%s" % (stats([solution.contacts for solution in cluster_solutions]), os.linesep))
                output.write(os.linesep)

        cluster_file_name = os.path.join(DEFAULT_SWARM_FOLDER + str(swarm_id), CLUSTER_REPRESENTATIVES_FILE)
        with open(cluster_file_name, 'w') as output:
            solutions_clustered = {}
            for id_solution, solution in enumerate(solutions):
                id_cluster = clusters[id_solution]
                try:
                    solutions_clustered[id_cluster].append(solution)
                except KeyError:
                    solutions_clustered[id_cluster] = [solution]

            solutions_reclustered = {}
            for k,v in zip(solutions_clustered.keys(), solutions_clustered.values()):
                if len(v) >= args.population_threshold:
                    solutions_reclustered[k] = v

            # Analysis
            num_clusters = len(solutions_reclustered.keys())
            population = [len(cluster) for cluster in solutions_reclustered.values()]
            threshold = np.mean(population) - np.std(population)
            log.info("Cluster population threshold is %d" % args.population_threshold)
            for id_cluster in sorted(solutions_reclustered, key=lambda k: len(solutions_reclustered[k]), reverse=True):
                cluster_solutions = solutions_reclustered[id_cluster]
                sorted_by_scoring = sorted(cluster_solutions, key=lambda k: k.scoring, reverse=True)
                if args.ponderated:
                    if len(sorted_by_scoring) >= threshold:
                        n = int(math.ceil((len(sorted_by_scoring) / float(len(solutions)))*num_clusters))
                        for i in range(n):
                            solution = sorted_by_scoring[i]
                            output.write("%d:%d:%8.5f:%d:%s" % (id_cluster,
                                                                len(cluster_solutions),
                                                                solution.scoring,
                                                                solution.id_glowworm,
                                                                solution.pdb_file) + os.linesep)
                else:
                    solution = sorted_by_scoring[0]
                    output.write("%d:%d:%8.5f:%d:%s" % (id_cluster,
                                                        len(cluster_solutions),
                                                        solution.scoring,
                                                        solution.id_glowworm,
                                                        solution.pdb_file) + os.linesep)


        log.info("Cluster: %d" % args.swarm_id)
        log.info("Number of steps: %d" % args.steps)
        log.info("Done.")

    except KeyboardInterrupt:
        log.info("Caught interrupt...")
        log.info("bye.")
