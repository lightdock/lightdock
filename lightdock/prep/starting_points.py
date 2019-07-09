"""Calculate the position of a set of points around a protein."""

import math
from scipy import spatial
import numpy as np
from lightdock.mathutil.cython.cutil import distance2


def points_on_sphere(number_of_points):
    """Creates a list of points using a spiral method.
    
    Based on method of 'Minimal Discrete Energy on the Sphere' (E. B. Saff, E.A. 
    Rakhmanov and Y.M. Zhou), Mathematical Research Letters, Vol. 1 (1994), pp. 647-662.
    
    Spiral method: Spiral from top of sphere to bottom of sphere, with points 
    places at distances the same as the distance between coils of the spiral.
    """
    points = []
    increment = math.pi * (3. - math.sqrt(5.))
    offset = 2./number_of_points
    for point in range(number_of_points):
        y = point * offset - 1.0 + (offset / 2.0)
        r = math.sqrt(1 - y*y)
        phi = point * increment
        points.append([math.cos(phi)*r, y, math.sin(phi)*r])
    return points


def calculate_surface_points(receptor, ligand, num_points, distance_step=0.5, is_membrane=False):
    """Calculates the position of num_points on the surface of the given protein.
    
    Uses a ray-tracing approach starting from the disposition of the points on the sphere 
    surface.
    """
    if num_points <= 0: 
        return []
    
    sphere_points = points_on_sphere(num_points)
    center_of_mass = receptor.center_of_mass()
    
    receptor_atom_coordinates = receptor.representative(is_membrane)

    distances_matrix_rec = spatial.distance.pdist(receptor_atom_coordinates)
    receptor_max_diameter = np.max(distances_matrix_rec)
    distances_matrix_lig = spatial.distance.pdist(ligand.representative())
    ligand_max_diameter = np.max(distances_matrix_lig)

    # Free memory if it's possible
    del distances_matrix_rec
    del distances_matrix_lig

    max_distance = receptor_max_diameter/2.0 + ligand_max_diameter/2.0
    for point in sphere_points:
        point[0] = point[0]*max_distance + center_of_mass[0]
        point[1] = point[1]*max_distance + center_of_mass[1]
        point[2] = point[2]*max_distance + center_of_mass[2]

    # Ray-tracing
    points = list(sphere_points)
    rays = []
    sub = np.subtract([center_of_mass for _ in range(len(sphere_points))], sphere_points)
    norms = [np.linalg.norm(v) for v in sub]
    for point, norm in zip(sub, norms):
        rays.append([point[0]/norm, point[1]/norm, point[2]/norm])
    
    surface_distance = (ligand_max_diameter/4.0)**2
    
    surface_points = []
    marked = []
    while len(surface_points) < num_points:
        for i, point in enumerate(points):
            if i not in marked:
                for atom in receptor.atoms:
                    if atom.residue_name != 'MMB':
                        if distance2(receptor_atom_coordinates[atom.index][0],
                                     receptor_atom_coordinates[atom.index][1],
                                     receptor_atom_coordinates[atom.index][2],
                                     point[0], point[1], point[2]) <= surface_distance:
                            marked.append(i)
                            surface_points.append(point)
                            break
                point[0] += rays[i][0] * distance_step
                point[1] += rays[i][1] * distance_step
                point[2] += rays[i][2] * distance_step
    
    return surface_points, receptor_max_diameter, ligand_max_diameter
