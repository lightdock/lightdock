import numpy as np
from lightdock.constants import (
    DEFAULT_LIGHTDOCK_PREFIX,
    DEFAULT_ELLIPSOID_DATA_EXTENSION,
    NUMPY_FILE_SAVE_EXTENSION,
)
from lightdock.gso.searchspace.ofunction import ObjectiveFunction


class ScoringFunction(ObjectiveFunction):
    """Scoring Functions interface"""

    def __init__(self, weight=1.0, anm_support=True):
        self.weight = float(weight)
        self.anm_support = anm_support

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        """Calculates the value of the scoring function.

        The GSO algorithm depends on a positive value for calculating the luciferin.
        If more negative means better in the scoring function, the sign must be changed.
        """
        raise NotImplementedError()

    @staticmethod
    def restraints_satisfied(restraints, interface):
        """Calculates the percentage of satisfied restraints"""
        if not restraints:
            return 0.0

        residues = list(restraints.keys())
        total = len(residues)
        satisfied = 0
        for residue in residues:
            intersection = set(restraints[residue]) & interface
            if len(intersection) > 0:
                satisfied += 1
        return float(satisfied) / total


class ModelAdapter(object):
    """Adapts a given Complex object as a DockingModel suitable for this
    ScoringFunction object.
    """

    def __init__(
        self, receptor, ligand, receptor_restraints=None, ligand_restraints=None
    ):
        self.receptor_model = self._get_docking_model(receptor, receptor_restraints)
        self.ligand_model = self._get_docking_model(ligand, ligand_restraints)

    def _get_docking_model(self, protein, restraints):
        """Complex -> DockingModel interface"""
        raise NotImplementedError()

    @staticmethod
    def load_reference_points(molecule):
        """Load reference points if exist"""
        reference_points = None
        try:
            ellipsoid_data_file = "%s%s%s" % (
                DEFAULT_LIGHTDOCK_PREFIX % molecule.structure_file_names[0],
                DEFAULT_ELLIPSOID_DATA_EXTENSION,
                NUMPY_FILE_SAVE_EXTENSION,
            )
            reference_points = np.load(ellipsoid_data_file)
        except (IOError, ValueError):
            pass
        return reference_points


# Two variables are needed to dynamically load the scoring functions from command line
DefinedScoringFunction = None
DefinedModelAdapter = None
