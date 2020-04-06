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
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF
from lightdock.error.lightdock_errors import NotSupportedInScoringError


log = LoggingManager.get_logger('dna')


class DNAModel(DockingModel):
    """Prepares the structure necessary for the scoring function"""
    def __init__(self, objects, coordinates, restraints, charges, vdw_energy, vdw_radii,
                 reference_points=None, n_modes=None):
        super(DNAModel, self).__init__(objects, coordinates, restraints, reference_points)
        self.charges = charges
        self.vdw_energy = vdw_energy
        self.vdw_radii = vdw_radii
        self.n_modes = n_modes

    def clone(self):
        """Creates a copy of the current model"""
        return DNAModel(self.objects, self.coordinates.copy(), self.restraints, 
                            self.charges, self.vdw_energy, self.vdw_radii,
                            reference_points=self.reference_points.copy())


class DNAAdapter(ModelAdapter):
    """Adapts a given Complex to a DockingModel object suitable for this scoring function."""

    to_translate = {'HIS':'HID', 'THY':'DT', 'ADE':'DA', 'CYT':'DC', 'GUA':'DG'}

    def _get_docking_model(self, molecule, restraints):
        atoms = molecule.atoms
        parsed_restraints = {}
        # Assign properties to atoms
        for atom_index, atom in enumerate(atoms):
            res_id = "%s.%s.%s" % (atom.chain_id, atom.residue_name, str(atom.residue_number))
            if restraints and res_id in restraints:
                try:
                    parsed_restraints[res_id].append(atom_index)
                except:
                    parsed_restraints[res_id] = [atom_index]
            try:
                res_name = atom.residue_name
                atom_name = atom.name
                if res_name in DNAAdapter.to_translate:
                    res_name = DNAAdapter.to_translate[res_name]
                if atom_name in amber.translate:
                    atom_name = amber.translate[atom.name]

                atom_id = "%s-%s" % (res_name, atom_name)
                try:
                    atom.amber_type = amber.amber_types[atom_id]
                except Exception as e:
                    # Maybe H N-terminal?
                    if atom_name in ['H1', 'H2', 'H3']:
                        atom_name = 'H'
                        atom_id = "%s-%s" % (res_name, atom_name)
                        atom.amber_type = amber.amber_types[atom_id]
                    else:
                        raise e
            except KeyError:
                raise NotSupportedInScoringError('Residue {} or atom {} not supported. '.format(res_id, atom_name) + 
                    'DNA scoring only supports AMBER94 types.')

            try:
                atom.charge = amber.charges[atom_id]
            except:
                # Go for N-terminal
                atom.charge = amber.nt_charges[atom_id]
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
            return DNAModel(atoms, coordinates, parsed_restraints, elec_charges, vdw_energies, vdw_radii,
                                    reference_points=reference_points, n_modes=molecule.n_modes.copy())
        except AttributeError:
            return DNAModel(atoms, coordinates, parsed_restraints, elec_charges, vdw_energies, vdw_radii,
                                    reference_points=reference_points)


class DNA(ScoringFunction):
    def __init__(self, weight=1.0):
        super(DNA, self).__init__(weight)
        try:
            with open(parameters.vdw_input_file) as vdw_file:
                self.scoring_vdw_weight = float(vdw_file.readline())
        except (IOError, ValueError) as e:
            log.warning('Error (%s), using default VdW cutoff' % str(e))
            self.scoring_vdw_weight = parameters.scoring_vdw_weight
        log.info('DNA VdW cutoff is: %3.2f' % self.scoring_vdw_weight)

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        """Computes the pyDockDNA scoring energy using receptor and ligand which are
        instances of DockingModel.
        """
        elec, vdw, interface_receptor, interface_ligand = cdna.calculate_energy(receptor_coordinates, ligand_coordinates,
                                                                                receptor.charges, ligand.charges,
                                                                                receptor.vdw_energy, ligand.vdw_energy,
                                                                                receptor.vdw_radii, ligand.vdw_radii,
                                                                                DEFAULT_CONTACT_RESTRAINTS_CUTOFF)
        energy = (elec + parameters.scoring_vdw_weight * vdw)*-1.
        perc_receptor_restraints = ScoringFunction.restraints_satisfied(receptor.restraints, set(interface_receptor))
        perc_ligand_restraints = ScoringFunction.restraints_satisfied(ligand.restraints, set(interface_ligand))
        return (energy + perc_receptor_restraints * energy + perc_ligand_restraints * energy) * self.weight


# Needed to dynamically load the scoring functions from command line
DefinedScoringFunction = DNA
DefinedModelAdapter = DNAAdapter
