"""SwarmDock energy scoring function.

Reference: SwarmDock and the use of Normal Models in Protein-Protein Docking
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2996808/

"""

import numpy as np
from lightdock.scoring.functions import ScoringFunction, ModelAdapter
from lightdock.structure.model import DockingModel
import lightdock.scoring.sd.energy.c.sd as sd
from lightdock.util.logger import LoggingManager
from lightdock.scoring.sd.data.amber import amber_types, masses, charges
import lightdock.scoring.sd.data.vdw as vdw
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF


log = LoggingManager.get_logger('sd')


class SDModel(DockingModel):
    """Prepares the structure necessary for the C-implementation"""
    def __init__(self, objects, coordinates, restraints, 
                 elec_charges, vdw_energy, vdw_radii, reference_points=None, n_modes=None):
        super(SDModel, self).__init__(objects, coordinates, restraints, reference_points)
        self.charges = elec_charges
        self.vdw_energy = vdw_energy
        self.vdw_radii = vdw_radii
        self.n_modes = n_modes

    def clone(self):
        """Creates a copy of the current model"""
        return SDModel(self.objects, self.coordinates.copy(), self.restraints, self.charges.copy(), self.vdw_energy.copy(),
                       self.vdw_radii.copy(), reference_points=self.reference_points.copy())


class SDAdapter(ModelAdapter):
    """Adapts a given Complex to a DockingModel object suitable for this scoring function."""
    def _get_docking_model(self, molecule, restraints):
        atoms = molecule.atoms
        parsed_restraints = {}
        # Assign properties to atoms
        for atom_index, atom in enumerate(atoms):
            res_name = atom.residue_name
            if res_name == "HIS":
                res_name = 'HID'
            res_id = "%s.%s.%s" % (atom.chain_id, atom.residue_name, str(atom.residue_number))
            if restraints and res_id in restraints:
                try:
                    parsed_restraints[res_id].append(atom_index)
                except:
                    parsed_restraints[res_id] = [atom_index]
            atom_id = "%s-%s" % (res_name, atom.name)
            atom.amber_type = amber_types[atom_id]
            atom.charge = charges[atom_id]
            atom.mass = masses[atom.amber_type]
            atom.vdw_energy = vdw.vdw_energy[atom.amber_type]
            atom.vdw_radius = vdw.vdw_radii[atom.amber_type]

        # Prepare common model information
        elec_charges = np.array([atom.charge for atom in atoms])
        vdw_energies = np.array([atom.vdw_energy for atom in atoms])
        vdw_radii = np.array([atom.vdw_radius for atom in atoms])
        coordinates = molecule.copy_coordinates()
        reference_points = ModelAdapter.load_reference_points(molecule)
        try:
            return SDModel(atoms, coordinates, parsed_restraints, elec_charges, vdw_energies, vdw_radii,
                           reference_points=reference_points, n_modes=molecule.n_modes.copy())
        except AttributeError:
            return SDModel(atoms, coordinates, parsed_restraints, elec_charges, vdw_energies, vdw_radii,
                           reference_points=reference_points)


class SD(ScoringFunction):

    def __init__(self, weight=1.0):
        super(SD, self).__init__(weight)

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        """Computes the SD scoring energy using receptor and ligand which are
        instances of DockingModel
        """
        energy, interface_receptor, interface_ligand = sd.calculate_energy(receptor_coordinates, ligand_coordinates,
                                                                           receptor.charges, ligand.charges,
                                                                           receptor.vdw_energy, ligand.vdw_energy,
                                                                           receptor.vdw_radii, ligand.vdw_radii, 
                                                                           DEFAULT_CONTACT_RESTRAINTS_CUTOFF)

        perc_receptor_restraints = ScoringFunction.restraints_satisfied(receptor.restraints, set(interface_receptor))
        perc_ligand_restraints = ScoringFunction.restraints_satisfied(ligand.restraints, set(interface_ligand))
        return (energy + perc_receptor_restraints * energy + perc_ligand_restraints * energy) * self.weight


# Needed to dynamically load the scoring functions from command line
DefinedScoringFunction = SD
DefinedModelAdapter = SDAdapter
