"""Module to deal with Bild geometric primitives"""

import os
from lightdock.mathutil.cython.quaternion import Quaternion


def sphere(center, radius):
    """Sphere primitive"""
    return ".sphere %f %f %f %f" % (center[0], center[1], center[2], radius)


def axis(pose, length=2):
    """Axis primitive"""
    q = Quaternion(pose[3], pose[4], pose[5], pose[6])

    bild = ".color cornflower blue" + os.linesep
    bild += sphere(pose, 0.3) + os.linesep

    bild += ".color 1 0 0" + os.linesep
    c = [length, 0, 0]
    end = q.rotate(c)
    bild += ".arrow %f %f %f %f %f %f%s" % (pose[0], pose[1], pose[2], pose[0] + end[0], pose[1] + end[1], pose[2] + end[2], os.linesep)

    bild += ".color 1 1 0" + os.linesep
    c = [0, length, 0]
    end = q.rotate(c)
    bild += ".arrow %f %f %f %f %f %f%s" % (pose[0], pose[1], pose[2], pose[0] + end[0], pose[1] + end[1], pose[2] + end[2], os.linesep)

    bild += ".color 0 0 1" + os.linesep
    c = [0, 0, length]
    end = q.rotate(c)
    bild += ".arrow %f %f %f %f %f %f%s" % (pose[0], pose[1], pose[2], pose[0] + end[0], pose[1] + end[1], pose[2] + end[2], os.linesep)

    return bild


def create_bild_file(file_name, poses):
    """Creates a Bild geometry find with the given poses and radius"""
    output = open(file_name, 'w')
    output.write(".color cornflower blue" + os.linesep)
    output.write(".transparency 0.7" + os.linesep)
    output.write(".transparency 0.0" + os.linesep)
    for pose in poses:
        output.write(axis(pose) + os.linesep)
    output.close()
