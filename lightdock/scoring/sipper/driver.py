"""SIPPER scoring function

C-implementation of the SIPPER (https://dx.doi.org/10.1021/ci100353e) scoring function.

"""

import numpy as np
import os
from lightdock.constants import DEFAULT_LIGHTDOCK_PREFIX
from lightdock.scoring.functions import ScoringFunction, ModelAdapter
from lightdock.structure.model import DockingModel
from lightdock.scoring.sipper.data.energy import sipper_energy, res_to_index
import lightdock.scoring.sipper.c.sipper as csipper
from lightdock.util.logger import LoggingManager


log = LoggingManager.get_logger("sipper")


class SIPPERModel(DockingModel):
    """Prepares the structure necessary for the C-implementation of the SIPPER
    scoring function
    """

    def __init__(
        self,
        objects,
        coordinates,
        restraints,
        energy,
        indexes,
        atoms_per_residue,
        oda=None,
        reference_points=None,
        n_modes=None,
    ):
        super(SIPPERModel, self).__init__(
            objects, coordinates, restraints, reference_points
        )
        self.energy = energy
        self.indexes = indexes
        self.atoms_per_residue = atoms_per_residue
        self.oda = oda
        self.n_modes = n_modes

    def clone(self):
        """Creates a copy of the current model"""
        return SIPPERModel(
            self.objects,
            self.coordinates.copy(),
            self.restraints,
            self.energy,
            self.indexes,
            self.atoms_per_residue,
            self.oda,
            reference_points=self.reference_points.copy(),
        )


class SIPPERAdapter(ModelAdapter):
    """Adapts a given Complex to a DockingModel object suitable for this
    SIPPER scoring function.
    """

    def _get_docking_model(self, molecule, restraints):
        atoms = molecule.atoms
        energy = sipper_energy
        parsed_restraints = {}
        indexes = np.array(
            [res_to_index[residue.name] for residue in molecule.residues]
        )
        coordinates = molecule.copy_coordinates()
        atoms_per_residue = np.array(
            [len(residue.atoms) for residue in molecule.residues]
        )

        for atom_index, atom in enumerate(atoms):
            res_id = f"{atom.chain_id}.{atom.residue_name}.{atom.residue_number}{atom.residue_insertion}"
            if restraints and res_id in restraints:
                try:
                    parsed_restraints[res_id].append(atom_index)
                except:
                    parsed_restraints[res_id] = [atom_index]

        oda_values = SIPPERAdapter._read_oda_values(molecule)
        if oda_values:
            oda_values = np.array(oda_values)
        else:
            # Dummy value, for performance
            oda_values = np.array([0.0 for _ in range(len(molecule.residues))])

        reference_points = ModelAdapter.load_reference_points(molecule)

        try:
            return SIPPERModel(
                atoms,
                coordinates,
                parsed_restraints,
                energy,
                indexes,
                atoms_per_residue,
                oda_values,
                reference_points=reference_points,
                n_modes=molecule.n_modes.copy(),
            )
        except AttributeError:
            return SIPPERModel(
                atoms,
                coordinates,
                parsed_restraints,
                energy,
                indexes,
                atoms_per_residue,
                oda_values,
                reference_points=reference_points,
            )

    @staticmethod
    def _read_oda_values(molecule):
        """Reads a ODA file with four fields per line. First field is the id of the residue,
        fourth field is the ODA value for that given residue.

        Returns None if no ODA file is found
        """
        oda_values = None
        oda_file_name = DEFAULT_LIGHTDOCK_PREFIX % (
            molecule.structure_file_names[0] + ".oda"
        )
        if os.path.exists(oda_file_name):
            log.info("ODA %s file found" % oda_file_name)
            with open(oda_file_name) as oda_input:
                oda_values = []
                lines = oda_input.readlines()[2:]
                for line in lines:
                    try:
                        fields = line.split()
                        oda_values.append(float(fields[3]))
                    except ValueError:
                        pass
        return oda_values


class SIPPER(ScoringFunction):
    def __init__(self, weight=1.0):
        super(SIPPER, self).__init__(weight)
        self.energy = sipper_energy

    def __call__(self, receptor, receptor_coordinates, ligand, ligand_coordinates):
        """Computes the pyDock scoring energy using receptor and ligand which are
        instances of DockingModel
        """
        energy, interface_receptor, interface_ligand = csipper.calculate_sipper(
            receptor_coordinates,
            ligand_coordinates,
            self.energy,
            receptor.indexes,
            ligand.indexes,
            receptor.atoms_per_residue,
            ligand.atoms_per_residue,
            len(receptor.atoms_per_residue),
            len(ligand.atoms_per_residue),
            receptor.oda,
            ligand.oda,
        )
        perc_receptor_restraints = ScoringFunction.restraints_satisfied(
            receptor.restraints, interface_receptor
        )
        perc_ligand_restraints = ScoringFunction.restraints_satisfied(
            ligand.restraints, interface_ligand
        )
        return (
            energy + perc_receptor_restraints * energy + perc_ligand_restraints * energy
        ) * self.weight


# Needed to dynamically load the scoring functions from command line
DefinedScoringFunction = SIPPER
DefinedModelAdapter = SIPPERAdapter
