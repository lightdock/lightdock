#!/usr/bin/env python

"""Cluster LightDock final swarm results using BSAS algorithm"""

import Bio.PDB
import sys
from lightdock.util.analysis import read_lightdock_output


def get_ca_atoms(ids_list):
    ca_atoms = {}
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    for struct_id in ids_list:
        pdb_file = "lightdock_%d.pdb" % struct_id
        print "Reading CA from %s" % pdb_file
        structure = pdb_parser.get_structure(pdb_file, pdb_file)
        model = structure[0]
        for chain in model:
            for residue in chain:
                try:
                    ca_atoms[struct_id].append(residue['CA'])
                except:
                    ca_atoms[struct_id] = [residue['CA']]
    return ca_atoms


def clusterize(sorted_ids):
    N = len(sorted_ids)
    super_imposer = Bio.PDB.Superimposer()

    clusters_found = 0
    clusters = {clusters_found : [sorted_ids[0]]}
    
    # Read all structures CA's
    ca_atoms = get_ca_atoms(sorted_ids)

    for j in sorted_ids[1:]:
        print "Glowworm %d with pdb lightdock_%d.pdb" % (j, j)
        in_cluster = False
        for cluster_id in clusters.keys():
            # For each cluster representative
            representative_id = clusters[cluster_id][0]
            super_imposer.set_atoms(ca_atoms[representative_id], ca_atoms[j])
            rmsd = super_imposer.rms
            print 'RMSD between %d and %d is %5.3f' % (representative_id, j, rmsd)
            if rmsd <= 4.0:
                clusters[cluster_id].append(j)
                print "Glowworm %d goes into cluster %d" % (j, cluster_id)
                in_cluster = True
                break
        
        if not in_cluster:
            clusters_found += 1
            clusters[clusters_found] = [j]
            print "New cluster %d" % clusters_found
    print clusters
    return clusters
    

def write_cluster_info(clusters, gso_data):
    with open('cluster.repr', 'w') as output:
        for id_cluster, ids in clusters.iteritems():
            output.write("%d:%d:%8.5f:%d:%s\n" % (id_cluster, len(ids), gso_data[ids[0]].scoring,
                                                  ids[0], 'lightdock_%d.pdb' % ids[0]))
        

if __name__ == '__main__':
 
    gso_data = read_lightdock_output(sys.argv[1])
    sorted_data = sorted(gso_data, key=lambda k: k.scoring, reverse=True)

    sorted_ids = [g.id_glowworm for g in sorted_data] 
    clusters = clusterize(sorted_ids)

    write_cluster_info(clusters, gso_data)
