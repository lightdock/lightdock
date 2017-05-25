"""Module to prepare initial poses for docking"""

import os
from lightdock.pdbutil.PDBIO import create_pdb_from_points
from lightdock.prep.starting_points import calculate_surface_points
from lightdock.prep.ftdock import FTDockCoordinatesParser, classify_ftdock_poses
from lightdock.mathutil.lrandom import MTGenerator, NormalGenerator
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.constants import CLUSTERS_CENTERS_FILE,\
    DEFAULT_PDB_STARTING_PREFIX, DEFAULT_STARTING_PREFIX, DEFAULT_BILD_STARTING_PREFIX, DEFAULT_EXTENT_MU, \
    DEFAULT_EXTENT_SIGMA
from lightdock.prep.geometry import create_bild_file


def get_random_point_within_sphere(number_generator, radius):
    """Generates a random point within a sphere of given radius"""
    r2 = radius**2
    while True:
        x = (2 * number_generator() - 1) * radius
        y = (2 * number_generator() - 1) * radius
        z = (2 * number_generator() - 1) * radius
        if x**2 + y**2 + z**2 <= r2:
            return x, y, z


def populate_poses(to_generate, center, radius, number_generator, rng_nm=None, rec_nm=0, lig_nm=0):
    """Creates new poses around a given center and a given radius"""
    new_poses = []
    for _ in xrange(to_generate):
        x, y, z = get_random_point_within_sphere(number_generator, radius)
        tx = center[0] + x
        ty = center[1] + y
        tz = center[2] + z
        q = Quaternion.random(number_generator)
        op_vector = [tx, ty, tz, q.w, q.x, q.y, q.z]
        if rng_nm:
            if rec_nm > 0:
                op_vector.extend([rng_nm() for _ in xrange(rec_nm)])
            if lig_nm > 0:
                op_vector.extend([rng_nm() for _ in xrange(lig_nm)])
        new_poses.append(op_vector)
    return new_poses


def create_file_from_poses(pos_file_name, poses):
    """Writes to file the initial poses"""
    positions_file = open(pos_file_name, 'w')
    for pose in poses:
        position = ''
        for coord in pose:
            position += str(coord) + ' '
        positions_file.write(position + os.linesep)
    positions_file.close()


def calculate_initial_poses(receptor, ligand, num_clusters, num_glowworms,
                            seed, dest_folder, ftdock_file='', nm_mode=False, nm_seed=0, rec_nm=0, lig_nm=0):
    """Calculates the starting points for each of the glowworms using the center of clusters
    and FTDock poses.
    """
    # Random number generator for poses
    rng = MTGenerator(seed)

    # Random number generator for NM
    if nm_mode:
        rng_nm = NormalGenerator(nm_seed, mu=DEFAULT_EXTENT_MU, sigma=DEFAULT_EXTENT_SIGMA)
    else:
        rng_nm = None
    
    # Calculate cluster centers
    cluster_centers, receptor_diameter, ligand_diameter = calculate_surface_points(receptor, 
                                                                                   ligand, 
                                                                                   num_clusters,
                                                                                   distance_step=1.0)
    pdb_file_name = os.path.join(dest_folder, CLUSTERS_CENTERS_FILE)
    create_pdb_from_points(pdb_file_name, cluster_centers)

    ligand_center = ligand.center_of_coordinates()
    radius = 10.    # ligand_diameter / 2.
    positions_files = []

    # Populate the clusters using the FTDock poses
    if ftdock_file:
        if nm_mode:
            raise NotImplementedError('Using FTDock poses with NM is not supported')

        poses = FTDockCoordinatesParser.get_list_of_poses(ftdock_file, ligand_center)
        clusters = classify_ftdock_poses(poses, cluster_centers, radius)

        for cluster_id, ftdock_poses in clusters.iteritems():
            # Translate FTDock poses into lightdock poses
            poses = []
            for pose in ftdock_poses:
                poses.append([pose.translation[0],
                              pose.translation[1],
                              pose.translation[2],
                              pose.q.w,
                              pose.q.x,
                              pose.q.y,
                              pose.q.z])

            # Populate new poses if needed
            if len(poses) < num_glowworms:
                needed = num_glowworms - len(poses)
                poses.extend(populate_poses(needed, cluster_centers[cluster_id], radius, rng))

            # Save poses as pdb file
            pdb_file_name = os.path.join(dest_folder, '%s_%s.pdb' % (DEFAULT_PDB_STARTING_PREFIX, cluster_id))
            create_pdb_from_points(pdb_file_name, [[pose[0], pose[1], pose[2]] for pose in poses[:num_glowworms]])

            # Save poses as initial_positions file
            pos_file_name = os.path.join(dest_folder, '%s_%s.dat' % (DEFAULT_STARTING_PREFIX, cluster_id))
            bild_file_name = os.path.join(dest_folder, '%s_%s.bild' % (DEFAULT_BILD_STARTING_PREFIX, cluster_id))
            create_file_from_poses(pos_file_name, poses[:num_glowworms])
            positions_files.append(pos_file_name)
            create_bild_file(bild_file_name, poses)
    else:
        for cluster_id, cluster_center in enumerate(cluster_centers):
            poses = populate_poses(num_glowworms, cluster_center, radius, rng, rng_nm, rec_nm, lig_nm)
            # Save poses as pdb file
            pdb_file_name = os.path.join(dest_folder, '%s_%s.pdb' % (DEFAULT_PDB_STARTING_PREFIX, cluster_id))
            create_pdb_from_points(pdb_file_name, [[pose[0], pose[1], pose[2]] for pose in poses[:num_glowworms]])

            # Save poses as initial_positions file
            pos_file_name = os.path.join(dest_folder, '%s_%s.dat' % (DEFAULT_STARTING_PREFIX, cluster_id))
            bild_file_name = os.path.join(dest_folder, '%s_%s.bild' % (DEFAULT_BILD_STARTING_PREFIX, cluster_id))
            create_file_from_poses(pos_file_name, poses[:num_glowworms])
            positions_files.append(pos_file_name)
            create_bild_file(bild_file_name, poses)

    return positions_files
