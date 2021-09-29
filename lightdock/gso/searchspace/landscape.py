"""Objective function landscape representations"""

import numpy as np
from scipy.optimize import fmin_powell

from lightdock.constants import (
    DEFAULT_STEP_SIZE,
    DEFAULT_TRANSLATION_STEP,
    DEFAULT_ROTATION_STEP,
)
from lightdock.mathutil.cython.quaternion import Quaternion


class LandscapePosition(object):
    """Represents glowworm's current position in the objective function space.

    Distance operation is defined for a Cartesian space of N dimensions.
    Different spaces should implement different approaches for move() and
    distance() functions.
    """

    def __init__(self, objective_function, coordinates, step=DEFAULT_STEP_SIZE):
        self.objective_function = objective_function
        self.coordinates = coordinates
        self.step = step

    def evaluate_objective_function(self):
        """Evaluates the objective function at the given coordinates"""
        return self.objective_function(self.coordinates)

    def __eq__(self, other):
        """Compares for equality two landscape positions"""
        return (
            self.coordinates == other.coordinates
            and self.objective_function == other.objective_function
        )

    def __ne__(self, other):
        """Compares for unequality two landscape positions"""
        return not self.__eq__(other)

    def clone(self):
        """Creates a copy of this landscape position"""
        return LandscapePosition(
            self.objective_function, self.coordinates.clone(), self.step
        )

    def __add__(self, other):
        """Adds two landscape positions"""
        return LandscapePosition(
            self.objective_function, self.coordinates + other.coordinates
        )

    def __iadd__(self, other):
        """Adds other to the current landscape position"""
        self.coordinates += other.coordinates
        return self

    def __sub__(self, other):
        """Subtracts two landscape positions"""
        return LandscapePosition(
            self.objective_function, self.coordinates - other.coordinates
        )

    def __isub__(self, other):
        """Subtracts other to the current landscape position"""
        self.coordinates -= other.coordinates
        return self

    def __mul__(self, scalar):
        """Multiplies this landscape position by a scalar"""
        return LandscapePosition(self.objective_function, self.coordinates * scalar)

    def norm(self):
        """Calculates the norm of the coordinates of this landscape position"""
        return self.coordinates.norm()

    def distance(self, other):
        """Calculates the distance between this landscape position and other"""
        delta_x = other - self
        return delta_x.coordinates.norm()

    def distance2(self, other):
        """Calculates the distance^2 between this landscape position and other"""
        delta_x = other - self
        return delta_x.coordinates.sum_of_squares()

    def move(self, other):
        """Move from this landscape position to another given a fixed step"""
        if self != other:
            delta_x = other - self
            delta_x *= self.step / delta_x.norm()
            self += delta_x
        return self

    def update_conformers(self, other, rnd_generator=None, current_scoring=0):
        """Compatibility with GSO test function tests"""
        pass

    def __repr__(self):
        return str(self.coordinates)


class DockingLandscapePosition(LandscapePosition):
    """Represents a current complex in the energy landscape.

    Receptor is fixed and ligand position and orientation depends on the current glowworm
    coordinates (optimization vector).
    """

    def __init__(
        self,
        scoring_function,
        coordinates,
        receptor,
        ligand,
        receptor_id=0,
        ligand_id=0,
        step_translation=DEFAULT_TRANSLATION_STEP,
        step_rotation=DEFAULT_ROTATION_STEP,
        step_nmodes=0,
        num_rec_nmodes=0,
        num_lig_nmodes=0,
    ):
        self.objective_function = scoring_function
        self.translation = np.array(coordinates[:3])
        self.rotation = Quaternion(
            coordinates[3], coordinates[4], coordinates[5], coordinates[6]
        )
        self.receptor = receptor
        self.ligand = ligand
        self.receptor_id = receptor_id
        self.ligand_id = ligand_id
        self.step_translation = step_translation
        self.step_rotation = step_rotation
        self.step_nmodes = step_nmodes
        self.num_rec_nmodes = num_rec_nmodes
        self.num_lig_nmodes = num_lig_nmodes
        # Copy ANM information if required
        self.rec_extent = (
            np.array(coordinates[7 : 7 + self.num_rec_nmodes])
            if self.num_rec_nmodes > 0
            else np.array([])
        )
        self.lig_extent = (
            np.array(coordinates[-self.num_lig_nmodes :])
            if self.num_lig_nmodes > 0
            else np.array([])
        )
        # This part is important, each position needs to retain its own pose coordinates
        self.receptor_pose = self.receptor.coordinates[self.receptor_id].clone()
        self.ligand_pose = self.ligand.coordinates[self.ligand_id].clone()
        self.ligand_reference_points = self.ligand.reference_points.clone()

    def clone(self):
        """Creates a copy of this landscape position"""
        coordinates = [
            self.translation[0],
            self.translation[1],
            self.translation[2],
            self.rotation.w,
            self.rotation.x,
            self.rotation.y,
            self.rotation.z,
        ]
        coordinates.extend(self.rec_extent)
        coordinates.extend(self.lig_extent)
        return DockingLandscapePosition(
            self.objective_function,
            coordinates,
            self.receptor,
            self.ligand,
            self.receptor_id,
            self.ligand_id,
            self.step_translation,
            self.step_rotation,
            self.step_nmodes,
            self.num_rec_nmodes,
            self.num_lig_nmodes,
        )

    def evaluate_objective_function(
        self, receptor_structure_id=None, ligand_structure_id=None
    ):
        """Evaluates the objective function at the given coordinates"""
        # Copy of the coordinates
        if receptor_structure_id:
            rec_id = receptor_structure_id
        else:
            rec_id = self.receptor_id
        self.receptor_pose = self.receptor.coordinates[rec_id].clone()
        if ligand_structure_id:
            lig_id = ligand_structure_id
        else:
            lig_id = self.ligand_id
        self.ligand_pose = self.ligand.coordinates[lig_id].clone()
        self.ligand_reference_points = self.ligand.reference_points.clone()

        # Use normal modes if provided:
        if self.num_rec_nmodes > 0:
            for i in range(self.num_rec_nmodes):
                # Only atoms as True in the mask are moved
                self.receptor_pose.coordinates[self.receptor.nm_mask, :] += (
                    self.receptor.n_modes[i] * self.rec_extent[i]
                )
        if self.num_lig_nmodes > 0:
            for i in range(self.num_lig_nmodes):
                # Only atoms as True in the mask are moved
                self.ligand_pose.coordinates[self.ligand.nm_mask, :] += (
                    self.ligand.n_modes[i] * self.lig_extent[i]
                )

        # We rotate first, ligand it's at initial position
        self.ligand_pose.rotate(self.rotation)
        self.ligand_reference_points.rotate(self.rotation)
        # Then translate
        self.ligand_pose.translate(self.translation)
        self.ligand_reference_points.translate(self.translation)
        return self.objective_function(
            self.receptor, self.receptor_pose, self.ligand, self.ligand_pose
        )

    def __eq__(self, other):
        """Compares for equality"""
        return (
            (self.translation == other.translation).all()
            and self.rotation == other.rotation
            and (self.rec_extent == other.rec_extent).all()
            and (self.lig_extent == other.lig_extent).all()
        )

    def distance(self, other):
        """Calculates the distance between this landscape position and other using reference points."""
        return np.sqrt(self.distance2(other))

    def distance2(self, other):
        """Calculates the distance^2 between this landscape position and other.

        ligand_pose has been already calculated in the update_luciferin
        stage of the algorithm.
        """
        rmsd2 = np.sum(
            (self.ligand_reference_points - other.ligand_reference_points) ** 2
        ) / len(self.ligand_reference_points)
        return rmsd2

    def move(self, other):
        """Move from this landscape position to another given a fixed step for translation
        and rotation movements.
        """
        if self != other:
            # Translation (Euclidian distance)
            delta_x = other.translation - self.translation
            n = np.linalg.norm(delta_x)
            # Only move if required
            if not np.allclose([0.0], [n]):
                delta_x *= self.step_translation / n
                self.translation += delta_x
            # Rotation (Quaternion SLERP)
            self.rotation = self.rotation.slerp(other.rotation, self.step_rotation)
            # NModes
            if self.num_rec_nmodes > 0:
                delta_x = other.rec_extent - self.rec_extent
                n = np.linalg.norm(delta_x)
                # Only move if required
                if not np.allclose([0.0], [n]):
                    delta_x *= self.step_nmodes / n
                    self.rec_extent += delta_x
            if self.num_lig_nmodes > 0:
                delta_x = other.lig_extent - self.lig_extent
                n = np.linalg.norm(delta_x)
                # Only move if required
                if not np.allclose([0.0], [n]):
                    delta_x *= self.step_nmodes / n
                    self.lig_extent += delta_x
        return self

    def update_conformers(self, other, rnd_generator, current_scoring):
        """Updates the structures for receptor and ligand"""
        if self != other:
            _ = rnd_generator.randint(upper_limit=(len(self.receptor) - 1))
            _ = rnd_generator.randint(upper_limit=(len(self.ligand) - 1))
            # Experimental, disabled
            # scoring = self.evaluate_objective_function(random_receptor_id, random_ligand_id)
            # if scoring > current_scoring:
            #    self.receptor_id = random_receptor_id
            #    self.ligand_id = random_ligand_id

    @staticmethod
    def _calculate_scoring(optimization_vector, self):
        """Calculates the energetic scoring at this current position.

        Required for local minimization"""
        self.update_landscape_position(optimization_vector)
        scoring = -1.0 * self.evaluate_objective_function()
        return scoring

    def update_landscape_position(self, optimized_vector):
        """Updates the current pose"""
        self.translation = optimized_vector[:3]
        self.rotation = Quaternion(
            optimized_vector[3],
            optimized_vector[4],
            optimized_vector[5],
            optimized_vector[6],
        )
        self.rec_extent = (
            optimized_vector[7 : 7 + self.num_rec_nmodes]
            if self.num_rec_nmodes > 0
            else np.array([])
        )
        self.lig_extent = (
            optimized_vector[-self.num_lig_nmodes :]
            if self.num_lig_nmodes > 0
            else np.array([])
        )

    def minimize(self):
        """Returns the new scoring after minimizing this landscape position using a local non-grandient
        minimization method.
        """
        optimization_vector = []
        optimization_vector.extend(self.translation)
        q = self.rotation
        optimization_vector.extend([q.w, q.x, q.y, q.z])
        optimization_vector.extend(self.rec_extent)
        optimization_vector.extend(self.lig_extent)
        optimization_vector = np.array(optimization_vector)

        # Minimize using Powell algorythm
        result = fmin_powell(
            DockingLandscapePosition._calculate_scoring,
            optimization_vector,
            args=(self,),
            maxiter=5,
            full_output=1,
            xtol=0.5,
            ftol=0.0001,
            disp=False,
        )
        # Update the landscape position vector
        optimized_vector = result[0]
        self.update_landscape_position(optimized_vector)
        # Return energy
        energy = -1.0 * result[1]
        return energy

    def __repr__(self):
        """String representation of this landscape position"""
        optimization_vector = list(self.translation) + [
            self.rotation.w,
            self.rotation.x,
            self.rotation.y,
            self.rotation.z,
        ]
        optimization_vector.extend(self.rec_extent)
        optimization_vector.extend(self.lig_extent)
        return "(%s) %4d %4d" % (
            ", ".join(["%10.7f" % v for v in optimization_vector]),
            self.receptor_id,
            self.ligand_id,
        )
