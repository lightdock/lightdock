"""Calculate the position of a set of points around a protein."""

import warnings
import itertools
from math import sqrt, cos, sin, pi, ceil
from pathlib import Path
import numpy as np
import freesasa
from scipy.cluster.vq import kmeans2
from scipy.spatial import distance, KDTree, ConvexHull
from prody import parsePDB, writePDB, calcDistance, calcCenter
from lightdock.constants import (
    STARTING_POINTS_SEED,
    DEFAULT_SWARM_DISTANCE,
    DEFAULT_SPHERES_PER_CENTROID,
    DEFAULT_CONTACT_RESTRAINTS_CUTOFF,
    DEFAULT_SWARMS_PER_RESTRAINT,
    DEFAULT_SWARM_RADIUS,
    SWARM_DISTANCE_TO_SURFACE_CUTOFF,
)
from lightdock.error.lightdock_errors import SetupError
from lightdock.util.logger import LoggingManager


log = LoggingManager.get_logger("lightdock3_setup")

freesasa.setVerbosity(freesasa.silent)


def points_in_hull(p, hull, tolerance=1e-12):
    """Calculates for a set of p points if they are inside of the convex hull"""
    return np.all(hull.equations[:,:-1] @ p.T + np.repeat(hull.equations[:,-1][None,:], len(p), axis=0).T <= tolerance, 0)

def equidistant_points(p1, p2, parts):
    """Calculates several points at equidistant distance between points p1 and p2"""
    return np.array(list(zip(np.linspace(p1[0], p2[0], parts+1),
                             np.linspace(p1[1], p2[1], parts+1),
                             np.linspace(p1[2], p2[2], parts+1)))
                    )


def points_on_sphere(number_of_points):
    """Creates a list of points using a spiral method.

    Based on method of 'Minimal Discrete Energy on the Sphere' (E. B. Saff, E.A.
    Rakhmanov and Y.M. Zhou), Mathematical Research Letters, Vol. 1 (1994), pp. 647-662.

    Spiral method: Spiral from top of sphere to bottom of sphere, with points
    places at distances the same as the distance between coils of the spiral.
    """
    points = []
    increment = pi * (3.0 - sqrt(5.0))
    offset = 2.0 / number_of_points
    for point in range(number_of_points):
        y = point * offset - 1.0 + (offset / 2.0)
        r = sqrt(1 - y * y)
        phi = point * increment
        points.append([cos(phi) * r, y, sin(phi) * r])
    return points


def calculate_surface_points(
    receptor,
    ligand,
    num_points,
    rec_translation,
    surface_density,
    receptor_restraints,
    blocking_restraints,
    seed=STARTING_POINTS_SEED,
    has_membrane=False,
    num_sphere_points=DEFAULT_SPHERES_PER_CENTROID,
    swarms_at_fixed_distance=DEFAULT_SWARM_DISTANCE,
    swarms_per_restraint=DEFAULT_SWARMS_PER_RESTRAINT,
    dense_sampling=False,
    verbose=True,
    probe_tolerance=0.5,
):
    """Calculates the position of num_points on the surface of the given protein.

    This new implementation differs from 0.9.0 series.

    Steps of the algorithm:
    1. Calculates several swarm centers using the surface atoms of the molecule and
    clustering using K-means.
    2. Places a sphere of points for each swarm center calculated in the previous step.
    3. Filters out overlapping points on those spheres.
    4. Filters out swarms too close to the molecule.
    5. If receptor restraints enabled, filters swarms which vision line goes through the
       receptor molecule.
    6. If not dense_sampling is enabled, clusters the final number of swarms in the
       given input number.

    """
    if num_points < 0:
        raise SetupError("Invalid number of points to generate over the surface")

    receptor_atom_coordinates = receptor.representative(has_membrane)

    # Calculate receptor and ligand max diameters
    distances_matrix_rec = distance.pdist(receptor_atom_coordinates)
    receptor_max_diameter = np.max(distances_matrix_rec)
    distances_matrix_lig = distance.pdist(ligand.representative())
    ligand_max_diameter = np.max(distances_matrix_lig)

    log.info(f"  * Ligand Max Diameter: {ligand_max_diameter:.2f} Å")
    if swarms_at_fixed_distance > 0.:
        # Fixed swarm distance to receptor's surface on user input
        surface_distance = swarms_at_fixed_distance
    else:
        if ligand_max_diameter < DEFAULT_SWARM_RADIUS*2:
            log.warning(f"Ligand radius is below the cutoff, using default swarm radius {DEFAULT_SWARM_RADIUS} as surface distance")
            surface_distance = DEFAULT_SWARM_RADIUS
        else:
            # We will use the ligand size to place the swarms over receptor's surface
            surface_distance = ligand_max_diameter / 4.0
    log.info(f"  * Surface distance: {surface_distance:.2f} Å")

    # Surface
    pdb_file_name = Path(receptor.structure_file_names[receptor.representative_id])
    molecule = parsePDB(pdb_file_name).select("protein or nucleic or (hetero and not water and not resname MMB)")
    if has_membrane:
        pdb_no_membrane = str(
            pdb_file_name.absolute().parent
            / f"{pdb_file_name.stem}_no_membrane{pdb_file_name.suffix}"
        )
        writePDB(pdb_no_membrane, molecule)

    if receptor_restraints:
        res_selection = " or ".join([f"chain {residue.get_chain()} and resnum {residue.number}" for residue in receptor_restraints])
        coords = molecule.select(f"within 10 of ({res_selection}) and surface or ({res_selection})").getCoords()
    else:
        surface = molecule.select("protein and surface or nucleic and name P or (hetero and not water and not resname MMB)")
        coords = surface.getCoords()

    # Automatic number of points
    if receptor_restraints:
        num_points = swarms_per_restraint * len(receptor_restraints)
    else:
        if num_points == 0:
            # Use SASA to get an estimation of points to calculate
            if has_membrane:
                structure = freesasa.Structure(pdb_no_membrane)
            else:
                structure = freesasa.Structure(str(pdb_file_name))
            result = freesasa.calc(structure)
            total_sasa = result.totalArea()
            # fix if using restraints
            num_points = ceil(total_sasa / surface_density)

    # Surface clusters
    if len(coords) > num_points:
        # Extremely important to set the seed in order to get reproducible results
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            surface_clusters = kmeans2(data=coords, k=num_points, minit="points", seed=seed, iter=100)
            surface_centroids = surface_clusters[0]
    else:
        surface_centroids = coords

    # Create points over the surface of each surface cluster
    sampling = []
    for sc in surface_centroids:
        sphere_points = np.array(points_on_sphere(num_sphere_points))
        surface_points = sphere_points * surface_distance + sc
        sampling.append(surface_points)

    # Filter out not compatible points
    centroids_kd_tree = KDTree(surface_centroids)
    for i_centroid in range(len(sampling)):
        # print('.', end="", flush=True)
        centroid = surface_centroids[i_centroid]
        # Search for this centroid neighbors
        centroid_neighbors = centroids_kd_tree.query_ball_point(centroid, r=20.0)
        # For each neighbor, remove points too close
        for n in centroid_neighbors:
            points_to_remove = []
            if n != i_centroid:
                for i_p, p in enumerate(sampling[i_centroid]):
                    if np.linalg.norm(p - surface_centroids[n]) <= surface_distance:
                        points_to_remove.append(i_p)
                points_to_remove = list(set(points_to_remove))
                sampling[i_centroid] = [
                    sampling[i_centroid][i_p]
                    for i_p in range(len(sampling[i_centroid]))
                    if i_p not in points_to_remove
                ]

    s = []
    for points in sampling:
        s.extend(points)

    if verbose:
        log.info(f"Swarms after incompatible filter: {len(s)}")

    # Filter interior points
    surface_swarms = []
    molecule_kd_tree = KDTree(molecule.getCoords())
    for swarm in s:
        if not molecule_kd_tree.query_ball_point(swarm, SWARM_DISTANCE_TO_SURFACE_CUTOFF):
            surface_swarms.append(swarm)
    s = surface_swarms

    if verbose:
        log.info(f"Swarms after interior points filter: {len(s)}")

    # Test for swarms which vision line goes through the receptor
    # We calculate for this the convex hull of surface atoms
    if receptor_restraints:
        # Calculate convex hull of surface atoms
        hull = ConvexHull(surface_centroids)
        swarms_with_visibility = []
        for swarm in s:
            # Find nearest restraint
            rst_centroids = [residue.get_central_atom() for residue in receptor_restraints]
            rst_centroid_coords = [[a.x, a.y, a.z] for a in rst_centroids]
            near_rst_centroid = np.array(rst_centroid_coords[np.argmin(distance.cdist(np.array([swarm]), rst_centroid_coords, 'euclidean'))])

            # Calculate probe positions between swarm and closest restraint
            p1 = np.array(swarm)
            p2 = near_rst_centroid

            num_probes = int(np.ceil(np.linalg.norm(p1 - p2)))
            p = equidistant_points(p1, p2, num_probes)

            # Calculate fraction of probes inside the convex hull described by the surface atoms
            a = points_in_hull(p, hull)

            #num_inside = np.count_nonzero(a)
            try:
                max_consecutive_inside = max([ sum( 1 for _ in group ) for key, group in itertools.groupby( a==True ) if key ])
            except ValueError:
                max_consecutive_inside = 0

            # Accept swarm depending on the ratio of probes found inside
            if max_consecutive_inside <= probe_tolerance * num_probes:
                swarms_with_visibility.append(swarm)

        # Update set of swarms
        s = swarms_with_visibility

        if verbose:
            log.info(f"Swarms after occlusion filter: {len(s)}")


    # Final cluster of points
    if len(s) > num_points and not dense_sampling:
        # Extremely important to set seed in order to get reproducible results
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            s_clusters = kmeans2(data=s, k=num_points, minit="++", seed=seed, iter=500)
            s = s_clusters[0]

    if dense_sampling:
        log.info("Dense sampling is enabled, ignoring user specified number of swarms")

    # Account for translation to origin of coordinates
    for p in s:
        p += rec_translation

    return s, receptor_max_diameter, ligand_max_diameter
