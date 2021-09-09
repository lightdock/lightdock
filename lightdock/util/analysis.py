import numpy as np
import operator
import os
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.constants import (
    RANKING_BY_LUCIFERIN_FILE,
    RANKING_BY_SCORING_FILE,
    RANKING_BY_RMSD_FILE,
    RANKING_FILE,
)
from lightdock.util.logger import LoggingManager


log = LoggingManager.get_logger("analysis")


class DockingResult(object):
    """Represents a LightDock docking result line"""

    def __init__(
        self,
        id_swarm=0,
        id_glowworm=0,
        receptor_id=0,
        ligand_id=0,
        luciferin=0.0,
        num_neighbors=0,
        vision_range=0.0,
        pose=None,
        rmsd=-1.0,
        pdb_file="",
        contacts=0,
        scoring=0.0,
    ):
        self.id_swarm = id_swarm
        self.id_glowworm = id_glowworm
        self.receptor_id = receptor_id
        self.ligand_id = ligand_id
        self.luciferin = luciferin
        self.num_neighbors = num_neighbors
        self.vision_range = vision_range
        self.pose = pose
        self.translation = np.array(pose[:3])
        self.rotation = Quaternion(pose[3], pose[4], pose[5], pose[6]).normalize()
        self.coord = DockingResult.pose_repr(pose)
        self.rmsd = rmsd
        self.pdb_file = pdb_file
        self.contacts = contacts
        self.scoring = scoring

    def __str__(self):
        return "%5d %6d %60s %6d %6d %11.5f %5d %7.3f %8.3f %16s %6d %8.3f" % (
            self.id_swarm,
            self.id_glowworm,
            self.coord,
            self.receptor_id,
            self.ligand_id,
            self.luciferin,
            self.num_neighbors,
            self.vision_range,
            self.rmsd,
            self.pdb_file,
            self.contacts,
            self.scoring,
        )

    def distance_trans(self, other):
        return np.linalg.norm(self.translation - other.translation)

    def distance_rot(self, other):
        return self.rotation.distance(other.rotation)

    @staticmethod
    def pose_repr(coord):
        fields = [("%5.3f" % c) for c in coord]
        return "(%s)" % (", ".join(fields))


def parse_coordinates(line):
    """Parses glowworm's coordinates found in line"""
    first = line.index("(")
    last = line.index(")")
    raw = line[first + 1 : last]
    coord = [float(c) for c in raw.split(",")]
    return coord, first, last


def read_lightdock_output(file_name, initial=None, final=None):
    """Reads a LightDock output file and sorts it by energy"""
    with open(file_name) as fin:
        raw_lines = [line for line in fin if line[0] != "#"]
        results = []
        for id_line, line in enumerate(raw_lines):
            try:
                coord, _, last = parse_coordinates(line)
            except ValueError:
                continue
            rest = line[last + 1 :].split()
            try:
                # Conformer solution
                result = DockingResult(
                    id_glowworm=id_line,
                    receptor_id=int(rest[0]),
                    ligand_id=int(rest[1]),
                    luciferin=float(rest[2]),
                    num_neighbors=int(rest[3]),
                    vision_range=float(rest[4]),
                    pose=coord,
                    scoring=float(rest[5]),
                )
            except ValueError:
                # Default solution
                result = DockingResult(
                    id_glowworm=id_line,
                    receptor_id=0,
                    ligand_id=0,
                    luciferin=float(rest[0]),
                    num_neighbors=int(rest[1]),
                    vision_range=float(rest[2]),
                    pose=coord,
                    scoring=float(rest[3]),
                )
            if initial and final:
                if (id_line + 2) > final:
                    break
                if (id_line + 1) >= initial:
                    results.append(result)
            else:
                results.append(result)

        return results


def read_ranking_file(ranking_file):
    """Reads a LightDock ranking file"""
    with open(ranking_file) as fin:
        raw_lines = fin.readlines()[1:]
        results = []
        for line in raw_lines:
            coord, first, last = parse_coordinates(line)
            raw_fields = line[:first].split()
            id_swarm = int(raw_fields[0])
            id_glowworm = int(raw_fields[1])
            raw_fields = line[last + 1 :].split()
            receptor_id = int(raw_fields[0])
            ligand_id = int(raw_fields[1])
            luciferin = float(raw_fields[2])
            num_neighbors = int(raw_fields[3])
            vision_range = float(raw_fields[4])
            rmsd = float(raw_fields[5])
            pdb_file = raw_fields[6]
            clashes = int(raw_fields[7])
            scoring = float(raw_fields[8])
            result = DockingResult(
                id_swarm=id_swarm,
                id_glowworm=id_glowworm,
                receptor_id=receptor_id,
                ligand_id=ligand_id,
                luciferin=luciferin,
                num_neighbors=num_neighbors,
                vision_range=vision_range,
                pose=coord,
                rmsd=rmsd,
                contacts=clashes,
                pdb_file=pdb_file,
                scoring=scoring,
            )
            results.append(result)
        return results


def write_ranking_to_file(solutions, clashes_cutoff=None, order_by=None):
    """Writes the calculated ranking to a file"""
    if order_by == "luciferin":
        output_file = RANKING_BY_LUCIFERIN_FILE
        solutions.sort(key=operator.attrgetter("luciferin"), reverse=True)
    elif order_by == "scoring":
        output_file = RANKING_BY_SCORING_FILE
        solutions.sort(key=operator.attrgetter("scoring"), reverse=True)
    elif order_by == "rmsd":
        output_file = RANKING_BY_RMSD_FILE
        solutions.sort(key=operator.attrgetter("rmsd"), reverse=False)
    else:
        output_file = RANKING_FILE

    output = open(output_file, "w")
    output.write(
        "Swarm  Glowworm   Coordinates                                             "
        "RecID  LigID  Luciferin  Neigh   VR     RMSD    PDB             Clashes  Scoring\n"
    )
    for solution in solutions:
        if clashes_cutoff:
            if solution.contacts <= clashes_cutoff:
                output.write("%s\n" % (str(solution)))
        else:
            output.write("%s\n" % (str(solution)))
    output.close()


def read_rmsd_and_contacts_data(file_name):
    """Reads a contacts file with columns identified by swarms, glowworm, num_contacts and rmsd"""
    contacts = {}
    rmsds = {}
    if os.path.isfile(file_name):
        with open(file_name) as fin:
            lines = [line.rstrip() for line in fin]
            # Ignore header
            for id_line, line in enumerate(lines[1:]):
                try:
                    fields = line.split()
                    swarm_id = int(fields[0])
                    structure_id = int(fields[1])
                    num_contacts = int(fields[2])
                    rmsd = float(fields[3])
                    if swarm_id in contacts:
                        contacts[swarm_id][structure_id] = num_contacts
                    else:
                        contacts[swarm_id] = {structure_id: num_contacts}
                    if swarm_id in rmsds:
                        rmsds[swarm_id][structure_id] = rmsd
                    else:
                        rmsds[swarm_id] = {structure_id: rmsd}
                except:
                    log.warning("Ignoring line %d in file %s" % (id_line, file_name))
    return contacts, rmsds


def read_cluster_representatives_file(cluster_file_name):
    """Reads a LightDock cluster representatives file"""
    with open(cluster_file_name) as fin:
        raw_lines = fin.readlines()
        glowworm_ids = []
        for line in raw_lines:
            line = line.rstrip(os.linesep)
            fields = line.split(":")
            glowworm_id = int(fields[3])
            glowworm_ids.append(glowworm_id)
        return glowworm_ids
