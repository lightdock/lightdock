#!/usr/bin/env python3

"""Filter LightDock final swarm results depending on the compatibility with the membrane"""


import sys
import os
import argparse
import shutil
import re
from prody.measure.contacts import Contacts
from prody import parsePDB, confProDy
from lightdock.util.logger import LoggingManager
from lightdock.util.analysis import read_ranking_file
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


# Disable ProDy output
confProDy(verbosity='info')
filtered_folder = 'filtered'

log = LoggingManager.get_logger('lgd_filter_membrane')


def get_structures(ranking, base_path='.'):
    structures = []
    for rank in ranking:
        swarm_id = rank.id_cluster
        glowworm_id = rank.id_glowworm
        structures.append(os.path.join(base_path, 
                                       'swarm_{}'.format(swarm_id), 
                                       'lightdock_{}.pdb'.format(glowworm_id)))
    return structures


def get_restraints(restraints_file):
    restraints_receptor = set()
    restraints_ligand = set()
    with open(restraints_file) as handle:
        for line in handle:
            line = line.rstrip(os.linesep)
            if line:
                if line.startswith('R'):
                    restraints_receptor.add(line.split()[-1])
                if line.startswith('L'):
                    restraints_ligand.add(line.split()[-1])
    return restraints_receptor, restraints_ligand


def calculate_membrane_height(parsed_receptor_file, restraints):
    atoms, residues, chains = parse_complex_from_file(parsed_receptor_file)
    receptor = Complex(chains, atoms)
    z_coord = []
    for restraint in restraints:
        chain_id, residue_name, residue_number = restraint.split(".")
        residue = receptor.get_residue(chain_id, residue_name, residue_number)
        ca = residue.get_calpha()
        z_coord.append(ca.z)
    return min(z_coord)


def parse_command_line():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(prog='lgd_filter_restraints')

    parser.add_argument("ranking_file", help="Path of ranking to be used", metavar="ranking_file")
    parser.add_argument("restraints_file", help="File including restraints", metavar="restraints_file")
    parser.add_argument("parsed_receptor_file", help="Receptor PDB parsed by LightDock", metavar="parsed_receptor_file")
    parser.add_argument("receptor_chains", help="Chains on the receptor partner", metavar="receptor_chains")
    parser.add_argument("ligand_chains", help="Chains on the receptor partner", metavar="ligand_chains")
    parser.add_argument("--cutoff", "-cutoff", "-c", help="Interaction cutoff",
                            dest="cutoff", type=float, default=1.0)

    return parser.parse_args()


if __name__ == '__main__':
 
    # Parse command line
    args = parse_command_line()

    log.info("Cutoff for membrane is {:3.1f}A".format(args.cutoff))

    # Get ranking
    ranking = read_ranking_file(args.ranking_file)

    # Get all the PDB structures in a given directory
    base_path = os.path.abspath(os.path.dirname(args.ranking_file))
    structures = get_structures(ranking, base_path)

    restraints_receptor, restraints_ligand = get_restraints(args.restraints_file)

    membrane_height_z = calculate_membrane_height(args.parsed_receptor_file, restraints_receptor)

    if os.path.exists(filtered_folder):
        raise SystemExit("Folder {} already exists".format(filtered_folder))
    else:
        os.makedirs(filtered_folder)

    filter_passed = {}
    percentages = {}
    for pdb_file in structures:
        try:
            swarm_id = int(re.findall(r'swarm_\d+', pdb_file)[0].split('_')[-1])
            glowworm_id = int(re.findall(r'lightdock_\d+', pdb_file)[0].split('_')[-1])

            # Read molecule and split by receptor and ligand
            molecule = parsePDB(pdb_file)
            ca_ligand = molecule.select('protein and chain {} and calpha'.format(args.ligand_chains))

            # Contacts on ligand side
            out = 0
            for ca in ca_ligand:
                coord = ca.getCoords()
                if coord[-1] >= membrane_height_z:
                    out += 1
            perc = out / float(len(ca_ligand))
            
            if perc >= args.cutoff:
                percentages[(swarm_id, glowworm_id)] = perc
                shutil.copyfile(pdb_file, os.path.join(filtered_folder, 'swarm_{}_{}.pdb'.format(swarm_id, glowworm_id)))
                try:
                    filter_passed[swarm_id].append(glowworm_id)
                except:
                    filter_passed[swarm_id] = [glowworm_id]
            print("{:40s}  {:5.3f}".format(pdb_file, perc))

        except Exception as e:
            log.error('Filtering has failed for structure {}. Please see error:'.format(pdb_file))
            log.error(str(e))


    filtered_ranking = os.path.join(filtered_folder, 'rank_filtered.list')
    with open(filtered_ranking, 'w') as handle:
        for rank in ranking:
            if rank.id_cluster in filter_passed and rank.id_glowworm in filter_passed[rank.id_cluster]:
                handle.write('swarm_{}_{}.pdb   {:5.3f}  {:5.3f}'.format(rank.id_cluster, 
                    rank.id_glowworm, rank.scoring, percentages[(rank.id_cluster, rank.id_glowworm)]) + os.linesep)
