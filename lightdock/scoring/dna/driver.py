"""Implementation of the pyDockDNA scoring function.

C-implementation of the pyDockDNA scoring function (no desolvation) and
custom VdW weight.
"""

import numpy as np
from lightdock.scoring.functions import ScoringFunction, ModelAdapter
from lightdock.structure.model import DockingModel
import lightdock.scoring.dna.energy.c.cdna as cdna
import lightdock.scoring.dna.energy.parameters as parameters
from lightdock.util.logger import LoggingManager
import lightdock.scoring.dna.data.amber as amber
import lightdock.scoring.dna.data.vdw as vdw

log = LoggingManager.get_logger('cdna')


class CPyDockDNAModel(DockingModel):
    """Prepares the structure necessary for the scoring function"""
    def __init__(self, objects, coordinates, restraints, charges, vdw_energy, vdw_radii,
                 reference_points=None, n_modes=None):
        super(CPyDockDNAModel, self).__init__(objects, coordinates, restraints, reference_points)
        self.charges = charges
        self.vdw_energy = vdw_energy
        self.vdw_radii = vdw_radii
        self.n_modes = n_modes

    def clone(self):
        """Creates a copy of the current model"""
        return CPyDockDNAModel(self.objects, self.coordinates.copy(), self.restraints, 
                            self.charges, self.vdw_energy, self.vdw_radii,
                            reference_points=self.reference_points.copy())


class CPyDockDNAAdapter(ModelAdapter):
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
        elec_charges = np.array([atom.charge for atom in atoms])
        vdw_energies = np.array([atom.vdw_energy for atom in atoms])
        vdw_radii = np.array([atom.vdw_radius for atom in atoms])
        coordinates = molecule.copy_coordinates()

        reference_points = ModelAdapter.load_reference_points(molecule)

        try:
            return CPyDockDNAModel(atoms, coordinates, restraints, elec_charges, vdw_energies, vdw_radii,
                                    reference_points=reference_points, n_modes=molecule.n_modes.copy())
        except AttributeError:
            return CPyDockDNAModel(atoms, coordinates, restraints, elec_charges, vdw_energies, vdw_radii,
                                    reference_points=reference_points)


class CPyDockDNA(ScoringFunction):
    def __init__(self, weight=1.0):
        super(CPyDockDNA, self).__init__(weight)
        try:
            with open(parameters.vdw_input_file) as vdw_file:
                self.scoring_vdw_weight = float(vdw_file.readline())
        except (IOError, ValueError) as e:
            log.warning('Error (%s), using default VDW cutoff' % str(e))
            self.scoring_vdw_weight = parameters.scoring_vdw_weight
        log.info('PyDockDNA VDW cutoff is: %3.2f' % self.scoring_vdw_weight)

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        """Computes the pyDockDNA scoring energy using receptor and ligand which are
        instances of DockingModel.
        """
        elec, vdw = cdna.calculate_energy(receptor_coordinates, ligand_coordinates,
                                          receptor.charges, ligand.charges,
                                          receptor.vdw_energy, ligand.vdw_energy,
                                          receptor.vdw_radii, ligand.vdw_radii)
        return (elec + parameters.scoring_vdw_weight * vdw)*-1.


# Needed to dynamically load the scoring functions from command line
DefinedScoringFunction = CPyDockDNA
DefinedModelAdapter = CPyDockDNAAdapter
