#!/bin/bash

################
# You may change these variables according to your needs
COMPLEX="4IZ7"
SWARMS=20
GLOWWORMS=50
STEPS=10
CORES=4
################

# Setup
lightdock_setup ${COMPLEX}_A_noh.pdb ${COMPLEX}_B_noh.pdb ${SWARMS} ${GLOWWORMS} --noxt --noh -anm -rst restraints.list

# Simulation
lightdock setup.json ${STEPS} -s fastdfire -c ${CORES}

# Generate predictions in PDB format
s=`ls -d ./swarm_* | wc -l`
swarms=$((s-1))
for i in $(seq 0 $swarms)
  do
    cd swarm_${i}; lgd_generate_conformations.py ../${COMPLEX}_A_noh.pdb ../${COMPLEX}_B_noh.pdb  gso_${STEPS}.out ${GLOWWORMS}; cd ..;
  done

# Cluster per swarm
for i in $(seq 0 $swarms)
  do
    cd swarm_${i}; lgd_cluster_bsas.py gso_${STEPS}.out; cd ..;
  done

# Generate ranking of predictions
lgd_rank.py ${s} ${STEPS}

# Filter according to 40% of satisfied restraints
lgd_filter_restraints.py -cutoff 5.0 -fnat 0.4 rank_by_scoring.list restraints.list A B

echo "Done."
