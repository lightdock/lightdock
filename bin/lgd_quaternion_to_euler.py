#!/usr/bin/env python3

"""Converts a Quaternion to a Rotation Matrix"""

import argparse
from lightdock.util.logger import LoggingManager
from lightdock.util.analysis import read_ranking_file
from lightdock.mathutil.cython.quaternion import Quaternion


log = LoggingManager.get_logger('quaternion2rot')


def parse_command_line():
    parser = argparse.ArgumentParser(prog='quaternion2rot')
    parser.add_argument("ranking_file", help="lightdock ranking file name")
    parser.add_argument("rotation_file", help="pyDock rotation file name")
    args = parser.parse_args()
    return args


def quaternion_to_matrix(q):
    x2 = q.x * q.x
    y2 = q.y * q.y
    z2 = q.z * q.z
    r11 = 1 - 2 * y2 - 2 * z2
    r12 = 2 * (q.x * q.y) + 2 * (q.w * q.z)
    r13 = 2 * (q.x * q.z) - 2 * (q.w * q.y)
    r21 = 2 * (q.x * q.y) - 2 * (q.w * q.z)
    r22 = 1 - 2 * x2 - 2 * z2
    r23 = 2 * (q.y * q.z) + 2 * (q.w * q.x)
    r31 = 2 * (q.x * q.z) + 2 * (q.w * q.y)
    r32 = 2 * (q.y * q.z) - 2 * (q.w * q.x)
    r33 = 1 - 2 * x2 - 2 * y2
    return r11, r12, r13, r21, r22, r23, r31, r32, r33


def create_rot_file(file_name, results):
    t = open(file_name, 'w')
    for i, result in enumerate(results):
        q = Quaternion(result.pose[3], result.pose[4], result.pose[5], result.pose[6])
        r11, r12, r13, r21, r22, r23, r31, r32, r33 = quaternion_to_matrix(q)
        t1 = result.pose[0]
        t2 = result.pose[1]
        t3 = result.pose[2]
        t.write('%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %6d\n' %
        (r11, r12, r13, r21, r22, r23, r31, r32, r33, t1, t2, t3, i + 1))
    t.close()

if __name__ == "__main__":
    try:
        # Parse command line
        args = parse_command_line()
        results = read_ranking_file(args.ranking_file)
        log.info('%d results read from %s' % (len(results), args.ranking_file))
        create_rot_file(args.rotation_file, results)
        log.info("Done.")

    except (KeyboardInterrupt):
        log.info("Caught interrupt...")
        log.info("bye.")
