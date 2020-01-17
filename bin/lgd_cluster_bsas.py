#!/usr/bin/env python3

"""Cluster LightDock final swarm results using BSAS algorithm"""

import argparse
import Bio.PDB
from lightdock.util.analysis import read_lightdock_output
from lightdock.util.logger import LoggingManager
from lightdock.constants import CLUSTER_REPRESENTATIVES_FILE


log = LoggingManager.get_logger('lgd_cluster_bsas')


def parse_command_line():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(prog='lgd_cluster_bsas')

    parser.add_argument("gso_output_file", help="LightDock output file", metavar="gso_output_file")

    return parser.parse_args()


def get_ca_atoms(ids_list):
    """Get all Carbon-alpha atoms of the PDB files specified by the ids_list.

    PDB files follow the format lightdock_ID.pdb where ID is in ids_list
    """
    ca_atoms = {}
    try:
        pdb_parser = Bio.PDB.PDBParser(QUIET=True)
        for struct_id in ids_list:
            pdb_file = "lightdock_%d.pdb" % struct_id
            log.info("Reading CA from %s" % pdb_file)
            structure = pdb_parser.get_structure(pdb_file, pdb_file)
            model = structure[0]
            for chain in model:
                for residue in chain:
                    try:
                        res = residue['CA']
                    except KeyError:
                        pass
                    try:
                        res = residue['P']
                    except KeyError:
                        pass
                    try:
                        ca_atoms[struct_id].append(res)
                    except KeyError:
                        ca_atoms[struct_id] = [res]
    except IOError as e:
        log.error('Error found reading a structure: %s' % str(e))
        log.error('Did you generate the LightDock structures corresponding to this output file?')
        raise SystemExit()

    return ca_atoms


def clusterize(sorted_ids):
    """Clusters the structures identified by the IDS inside sorted_ids list"""
    super_imposer = Bio.PDB.Superimposer()

    clusters_found = 0
    clusters = {clusters_found: [sorted_ids[0]]}

    # Read all structures CA's
    ca_atoms = get_ca_atoms(sorted_ids)

    for j in sorted_ids[1:]:
        log.info("Glowworm %d with pdb lightdock_%d.pdb" % (j, j))
        in_cluster = False
        for cluster_id in list(clusters.keys()):
            # For each cluster representative
            representative_id = clusters[cluster_id][0]
            super_imposer.set_atoms(ca_atoms[representative_id], ca_atoms[j])
            rmsd = super_imposer.rms
            log.info('RMSD between %d and %d is %5.3f' % (representative_id, j, rmsd))
            if rmsd <= 4.0:
                clusters[cluster_id].append(j)
                log.info("Glowworm %d goes into cluster %d" % (j, cluster_id))
                in_cluster = True
                break

        if not in_cluster:
            clusters_found += 1
            clusters[clusters_found] = [j]
            log.info("New cluster %d" % clusters_found)
    return clusters


def write_cluster_info(clusters, gso_data):
    """Writes the clustering result"""
    with open(CLUSTER_REPRESENTATIVES_FILE, 'w') as output:
        for id_cluster, ids in clusters.items():
            output.write("%d:%d:%8.5f:%d:%s\n" % (id_cluster, len(ids), gso_data[ids[0]].scoring,
                                                  ids[0], 'lightdock_%d.pdb' % ids[0]))
        log.info("Cluster result written to %s file" % CLUSTER_REPRESENTATIVES_FILE)


if __name__ == '__main__':

    try:
        # Parse command line
        args = parse_command_line()

        # Read LightDock output data
        gso_data = read_lightdock_output(args.gso_output_file)

        # Sort the glowworms data by scoring
        sorted_data = sorted(gso_data, key=lambda k: k.scoring, reverse=True)

        # Get the Glowworm ids sorted by their scoring
        sorted_ids = [g.id_glowworm for g in sorted_data]

        # Calculate the different clusters
        clusters = clusterize(sorted_ids)

        # Write clustering information
        write_cluster_info(clusters, gso_data)

    except Exception as e:
        log.error('Clustering has failed. Please see error:')
        log.error(str(e))
