
import numpy as np
from os import linesep
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.error.lightdock_errors import GSOCoordinatesError
from lightdock.util.logger import LoggingManager


log = LoggingManager.get_logger("analysis")


class LightDockPose:
    """Represents a LightDock pose in translational, rotational and ANM spaces"""
    def __init__(self,
        translation=None,
        rotation=None,
        rec_extents=None,
        lig_extents=None,
    ):
        # Check for NaN
        if np.isnan(np.array(translation, dtype=float)).any() or \
            np.isnan(np.array(rotation, dtype=float)).any() or \
            np.isnan(np.array(rec_extents, dtype=float)).any() or\
            np.isnan(np.array(lig_extents, dtype=float)).any():
            raise GSOCoordinatesError('NaN found in pose')

        self.translation = np.array(translation, dtype=float)
        self.rotation = Quaternion(rotation[0], rotation[1], rotation[2], rotation[3]).normalize()
        self.rec_extents = np.array(rec_extents, dtype=float)
        self.lig_extents = np.array(lig_extents, dtype=float)

    def distance_trans(self, other):
        """Translation distance between two predictions"""
        return np.linalg.norm(self.translation - other.translation)

    def distance_rot(self, other):
        """Quaternion distance between two predictions"""
        return self.rotation.distance(other.rotation)

    def pose_repr(self, precision='.9f'):
        """Pose representation"""

        pose = f'{self.translation[0]:{precision}}, {self.translation[1]:{precision}}, {self.translation[2]:{precision}}'
        pose += f', {self.rotation.w:{precision}}, {self.rotation.x:{precision}}, {self.rotation.y:{precision}}, {self.rotation.z:{precision}}'

        if self.rec_extents is not None and self.rec_extents.any():
            pose += ', ' + ', '.join(np.char.mod(f'%{precision}', self.rec_extents))

        if self.lig_extents is not None and self.lig_extents.any():
            pose += ', ' + ', '.join(np.char.mod(f'%{precision}', self.lig_extents))

        return f'({pose})'


def read_initial_positions_file(file_name, num_anm_rec=0, num_anm_lig=0):
    """Reads and parses the LightDock poses stored in an initial positions file.
    Returns a list of LightDockPose sorted by the same order as in the provided file."""
    poses = []
    try:
        with open(file_name) as ih:
            line_count = 0
            for line in ih:
                line_count += 1
                line = line.rstrip(linesep)
                if line and line[0] != "#":
                    coord = line.split()
                    # Check for length
                    if len(coord) != (7 + num_anm_rec + num_anm_lig):
                        raise GSOCoordinatesError(f'Malformed line {line_count} in {file_name}')

                    try:
                        translation = [float(coord[0]), float(coord[1]), float(coord[2])]

                        rotation = [float(coord[3]), float(coord[4]), float(coord[5]), float(coord[6])]

                        rec_extents = []
                        lig_extents = []
                        if len(coord) > 7:
                            rec_extents = [float(x) for x in coord[7 : 7 + num_anm_rec]]
                            lig_extents = [float(x) for x in coord[-num_anm_lig:]]

                        pose = LightDockPose(translation, rotation, rec_extents, lig_extents)
                        poses.append(pose)

                    except ValueError:
                        raise GSOCoordinatesError(f'Malformed pose in line {line_count}')

            log.info(f"Read {line_count} poses")

    except FileNotFoundError:
        raise GSOCoordinatesError(f'Cannot find or open {file_name}')

    return poses


class LightDockPrediction:
    """Represents a LightDock docking prediction"""

    def __init__(
        self,
        swarm=0,
        glowworm=0,
        pose=None,
        receptor_id=0,
        ligand_id=0,
        luciferin=0.0,
        num_neighbors=0,
        vision_range=0.0,
        scoring=0.0,
    ):
        self.swarm = swarm
        self.glowworm = glowworm
        self.receptor_id = receptor_id
        self.ligand_id = ligand_id
        self.luciferin = luciferin
        self.num_neighbors = num_neighbors
        self.vision_range = vision_range
        self.pose = pose
        self.scoring = scoring


def read_predictions_file(file_name, num_anm_rec=0, num_anm_lig=0):
    """Reads and parses the LightDock predictions.
    Returns a list of LightDockPrediction sorted by the same order as in the provided file."""
    predictions = []
    num_predictions = 0
    try:
        with open(file_name) as ih:
            line_count = 0
            for line in ih:
                line_count += 1
                line = line.rstrip(linesep)
                if line and line[0] != "#":
                    # First, parse coordinates and create a new LightDockPose
                    first = line.index("(")
                    last = line.index(")")
                    coord = [float(c) for c in line[first + 1 : last].split(",")]
                    if len(coord) != (7 + num_anm_rec + num_anm_lig):
                        raise GSOCoordinatesError(f'Malformed line {line_count} in {file_name}')

                    pose = LightDockPose(translation=coord[:3],
                        rotation=coord[3:7],
                        rec_extents=coord[7 : 7 + num_anm_rec],
                        lig_extents=coord[-num_anm_lig:])

                    # Parse the rest of fields
                    try:
                        fields = line[last+1:].split()
                        receptor_id = int(fields[0])
                        ligand_id = int(fields[1])
                        luciferin = float(fields[2])
                        num_neighbors = int(fields[3])
                        vision_range = float(fields[4])
                        scoring = float(fields[5])

                        prediction = LightDockPrediction(pose=pose, receptor_id=receptor_id, ligand_id=ligand_id,
                            luciferin=luciferin, num_neighbors=num_neighbors, vision_range=vision_range, scoring=scoring)
                        predictions.append(prediction)
                        num_predictions += 1

                    except ValueError:
                        raise GSOCoordinatesError(f'Malformed prediction in line {line_count}')

            log.info(f"Read {num_predictions} predictions")

    except FileNotFoundError:
        raise GSOCoordinatesError(f'Cannot find or open {file_name}')

    return predictions
