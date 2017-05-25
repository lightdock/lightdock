"""Module to parse and prepare FTDock poses to be used in lightdock as
starting points/poses.
"""

import os
from math import cos, sin, sqrt
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.mathutil.constants import DEG2RAD
from lightdock.mathutil.cython.cutil import distance2


def classify_ftdock_poses(poses, cluster_centers, radius):
    """Classifies the different FTDock poses depending on the proximity (radius)
    to each of the cluster_centers points.
    """
    dist = radius**2
    clusters = {}
    for pose in poses:
        for id_cluster, center in enumerate(cluster_centers):
            if distance2(center[0], center[1], center[2],
                         pose.translation[0], pose.translation[1], pose.translation[2]) <= dist:
                try:
                    clusters[id_cluster].append(pose)
                except:
                    clusters[id_cluster] = [pose]
    return clusters


class FTDockCoordinatesParser(object):
    """Parses a FTDock result file"""
    @staticmethod
    def get_list_of_poses(ftdock_file, ligand_center):
        """Reads a FTDock file and parses the found poses"""
        span = 1.0
        poses = []
        
        lines = open(ftdock_file).readlines()
        
        for line in lines:
            if line[:6] == "G_DATA":
                line = line.rstrip(os.linesep)
                fields = line.split()
                translation = [float(x)*span for x in fields[5:8]]
                rotation = [float(x) for x in fields[8:11]]
                poses.append(FTDockPose(rotation[0],
                                                rotation[1],
                                                rotation[2],
                                                translation[0],
                                                translation[1],
                                                translation[2]))
            elif line[:21] == 'Global grid cell span':
                span = float(line.split()[-1])
        
        # Post-process of the translations
        for pose in poses:
            offset = pose.q.rotate(ligand_center)
            pose.translation = [x-y for x,y in zip(pose.translation, offset)]
        
        return poses


class FTDockPose(object):
    """Represents a FTDock position and rotation line"""
    def __init__(self, psi, theta, phi, tx, ty, tz):
        self.psi = psi
        self.theta = theta
        self.phi = phi
        self.translation = [tx, ty, tz]
        self.r = self.calculate_rotation_matrix()
        self.q = self.calculate_quaternion()

    def calculate_rotation_matrix(self):
        """Calculates a rotation matrix from the Euler angles"""
        r = [[0.0 for _ in range(3)] for _ in range(3)]
        r[0][0] = -sin(self.psi*DEG2RAD)*sin(DEG2RAD*self.phi) + cos(DEG2RAD*self.psi)*cos(DEG2RAD*self.theta)*cos(DEG2RAD*self.phi)
        r[1][0] = sin(DEG2RAD*self.psi)*cos(DEG2RAD*self.phi) + cos(DEG2RAD*self.psi)*cos(DEG2RAD*self.theta)*sin(DEG2RAD*self.phi)
        r[2][0] = -cos(DEG2RAD*self.psi)*sin(DEG2RAD*self.theta)
        r[0][1] = -cos(DEG2RAD*self.psi)*sin(DEG2RAD*self.phi) - sin(DEG2RAD*self.psi)*cos(DEG2RAD*self.theta)*cos(DEG2RAD*self.phi)
        r[1][1] = cos(DEG2RAD*self.psi)*cos(DEG2RAD*self.phi) - sin(DEG2RAD*self.psi)*cos(DEG2RAD*self.theta)*sin(DEG2RAD*self.phi)
        r[2][1] = sin(DEG2RAD*self.psi)*sin(DEG2RAD*self.theta)
        r[0][2] = sin(DEG2RAD*self.theta)*cos(DEG2RAD*self.phi)
        r[1][2] = sin(DEG2RAD*self.theta)*sin(DEG2RAD*self.phi)
        r[2][2] = cos(DEG2RAD*self.theta)
        return r

    def calculate_quaternion(self):
        """Calculates a quaternion given a rotation matrix. Two quaternions can express the
        same rotation (rotation space is double-covered).
        
        See: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
        """
        trace = self.r[0][0] + self.r[1][1] + self.r[2][2]
        if trace > 0.0:
            s = 0.5 / sqrt(trace + 1.0)
            w = 0.25 / s
            x = (self.r[2][1] - self.r[1][2]) * s
            y = (self.r[0][2] - self.r[2][0]) * s
            z = (self.r[1][0] - self.r[0][1]) * s
            return Quaternion(w,x,y,z)
        else:
            if self.r[0][0] > self.r[1][1] and self.r[0][0] > self.r[2][2]:
                s = 2.0 * sqrt(1.0 + self.r[0][0] - self.r[1][1] - self.r[2][2])
                w = (self.r[2][1] - self.r[1][2]) / s
                x = 0.25 * s
                y = (self.r[0][1] + self.r[1][0]) / s
                z = (self.r[0][2] + self.r[2][0]) / s
                return Quaternion(w,x,y,z)
            elif self.r[1][1] > self.r[2][2]:
                s = 2.0 * sqrt(1.0 + self.r[1][1] - self.r[0][0] - self.r[2][2])
                w = (self.r[0][2] - self.r[2][0]) / s
                x = (self.r[0][1] + self.r[1][0]) / s
                y = 0.25 * s
                z = (self.r[1][2] + self.r[2][1]) / s
                return Quaternion(w,x,y,z)
            else:
                s = 2.0 * sqrt(1.0 + self.r[2][2] - self.r[0][0] - self.r[1][1])
                w = (self.r[1][0] - self.r[0][1]) / s
                x = (self.r[0][2] + self.r[2][0]) / s
                y = (self.r[1][2] + self.r[2][1]) / s
                z = 0.25 * s
                return Quaternion(w,x,y,z)
        
    def __repr__(self):
        return "Psi:%5.3f Theta:%5.3f Phi:%5.3f %s" % (self.psi,
                                                       self.theta,
                                                       self.phi,
                                                       self.translation)
