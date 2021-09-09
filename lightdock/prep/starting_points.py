"""Calculate the position of a set of points around a protein."""

from math import sqrt, cos, sin, pi, ceil
from pathlib import Path
import numpy as np
import freesasa
from scipy.cluster.vq import kmeans2
from scipy.spatial import distance, KDTree
from prody import parsePDB, writePDB
from lightdock.constants import STARTING_POINTS_SEED
from lightdock.error.lightdock_errors import SetupError

freesasa.setVerbosity(freesasa.silent)


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
    seed=STARTING_POINTS_SEED,
    has_membrane=False,
    num_sphere_points=100,
):
    """Calculates the position of num_points on the surface of the given protein"""
    if num_points < 0:
        raise SetupError("Invalid number of points to generate over the surface")

    receptor_atom_coordinates = receptor.representative(has_membrane)

    distances_matrix_rec = distance.pdist(receptor_atom_coordinates)
    receptor_max_diameter = np.max(distances_matrix_rec)
    distances_matrix_lig = distance.pdist(ligand.representative())
    ligand_max_diameter = np.max(distances_matrix_lig)
    surface_distance = ligand_max_diameter / 4.0

    # Surface
    pdb_file_name = Path(receptor.structure_file_names[receptor.representative_id])
    molecule = parsePDB(pdb_file_name).select("protein or nucleic")
    if has_membrane:
        pdb_no_membrane = str(
            pdb_file_name.absolute().parent
            / f"{pdb_file_name.stem}_no_membrane{pdb_file_name.suffix}"
        )
        writePDB(pdb_no_membrane, molecule)
    surface = molecule.select("protein and surface or nucleic and name P")
    coords = surface.getCoords()

    # SASA
    if num_points == 0:
        if has_membrane:
            structure = freesasa.Structure(pdb_no_membrane)
        else:
            structure = freesasa.Structure(str(pdb_file_name))
        result = freesasa.calc(structure)
        total_sasa = result.totalArea()
        num_points = ceil(total_sasa / surface_density)

    # Surface clusters
    if len(coords) > num_points:
        # Extremely important to set seed in order to get reproducible results
        np.random.seed(seed)
        surface_clusters = kmeans2(data=coords, k=num_points, minit="points", iter=100)
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

    # Final cluster of points
    if len(s) > num_points:
        # Extremely important to set seed in order to get reproducible results
        np.random.seed(seed)
        s_clusters = kmeans2(data=s, k=num_points, minit="points", iter=100)
        s = s_clusters[0]

    for p in s:
        p += rec_translation

    return s, receptor_max_diameter, ligand_max_diameter
