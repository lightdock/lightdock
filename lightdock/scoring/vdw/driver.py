"""Truncated Van der Waals energy"""

import numpy as np

from lightdock.scoring.functions import ScoringFunction, ModelAdapter
from lightdock.scoring.vdw.data import amber
from lightdock.scoring.vdw.data import vdw
from lightdock.structure.model import DockingModel
import lightdock.scoring.vdw.energy.c.cvdw as cvdw
from lightdock.util.logger import LoggingManager


log = LoggingManager.get_logger('vdw')


class VdWModel(DockingModel):
    """Prepares the necessary structure"""
    def __init__(self, objects, coordinates, vdw_energy, vdw_radii, reference_points=None, n_modes=None):
        super(VdWModel, self).__init__(objects, coordinates, reference_points)
        self.vdw_energy = vdw_energy
        self.vdw_radii = vdw_radii
        self.n_modes = n_modes

    def clone(self):
        """Creates a copy of the current model"""
        return VdWModel(self.objects, self.coordinates.copy(), self.vdw_energy, self.vdw_radii,
                        reference_points=self.reference_points.copy())


class VdWAdapter(ModelAdapter):
    """Adapts a given Complex to a DockingModel object suitable for this scoring function."""
    def _get_docking_model(self, molecule, restraints):
        atoms = molecule.atoms
        # Assign properties to atoms
        for atom in atoms:
            res_name = atom.residue_name
            atom_name = atom.name
            if res_name == "HIS":
                res_name = 'HID'
            if atom_name in amber.translate:
                atom_name = amber.translate[atom.name]
            atom_id = "%s-%s" % (res_name, atom_name)
            atom.amber_type = amber.amber_types[atom_id]
            atom.charge = amber.charges[atom_id]
            atom.mass = amber.masses[atom.amber_type]
            atom.vdw_energy = vdw.vdw_energy[atom.amber_type]
            atom.vdw_radius = vdw.vdw_radii[atom.amber_type]

        # Prepare common model information
        vdw_energies = np.array([atom.vdw_energy for atom in atoms])
        vdw_radii = np.array([atom.vdw_radius for atom in atoms])
        coordinates = molecule.copy_coordinates()

        reference_points = ModelAdapter.load_reference_points(molecule)

        try:
            return VdWModel(atoms, coordinates, vdw_energies, vdw_radii,
                            reference_points=reference_points, n_modes=molecule.n_modes.copy())
        except AttributeError:
            return VdWModel(atoms, coordinates, vdw_energies, vdw_radii, reference_points=reference_points)


class VdW(ScoringFunction):
    def __init__(self, weight=1.0):
        super(VdW, self).__init__(weight)

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        """Computes the truncated VdW energy using receptor and ligand which are
        instances of the DockingModel class"""
        vdw_energy = cvdw.calculate_vdw(receptor_coordinates, ligand_coordinates,
                                        receptor.vdw_energy, ligand.vdw_energy, receptor.vdw_radii, ligand.vdw_radii)
        return vdw_energy * -1.0 * self.weight


# Needed to dynamically load the scoring functions from command line
DefinedScoringFunction = VdW
DefinedModelAdapter = VdWAdapter
