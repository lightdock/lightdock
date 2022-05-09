"""Module to calculate normal modes of a given protein.

It uses the awesome Prody library
"""

import numpy as np
from prody import parsePDB, ANM, extendModel, confProDy
from lightdock.error.lightdock_errors import NormalModesCalculationError
from lightdock.util.logger import LoggingManager
from lightdock.constants import (
    DEFAULT_NUM_SWARMS,
    DEFAULT_NUM_GLOWWORMS,
    DEFAULT_EXTENT_MU,
    DEFAULT_EXTENT_SIGMA,
)

# Disable ProDy output
confProDy(verbosity="none")

log = LoggingManager.get_logger("ANM")


def calculate_nmodes(pdb_file_name, n_modes, rmsd, seed, molecule):
    """Calculates Normal modes for a given molecule"""
    prody_molecule = parsePDB(pdb_file_name)
    # Try first for proteins
    backbone_atoms = prody_molecule.select("name CA")
    if not backbone_atoms:
        # If previous step has failed, maybe we're dealing with DNA
        backbone_atoms = prody_molecule.select("name C4'")
    if not backbone_atoms:
        raise NormalModesCalculationError(
            "Error selecting backbone atoms (protein or nucleic)"
        )
    molecule_anm = ANM("molecule backbone")
    molecule_anm.buildHessian(backbone_atoms)
    molecule_anm.calcModes(n_modes=n_modes)

    num_atoms_prody = prody_molecule.select("protein or nucleic").numAtoms()

    if num_atoms_prody != (molecule.protein_num_atoms + molecule.nucleic_num_atoms):
        log.error(
            "Number of atoms in ProDy (%d) vs LightDock (%d)"
            % (num_atoms_prody, molecule.protein_num_atoms + molecule.nucleic_num_atoms)
        )
        raise NormalModesCalculationError("Number of atoms is different")

    # Check for sanity in atoms from both structures just in case
    prody_molecule_prot_nuc = prody_molecule.select("protein or nucleic")
    for lightdock_atom, prody_atom in zip(
        molecule.get_atoms(protein=True, nucleic=True, dummy=False),
        prody_molecule_prot_nuc.iterAtoms(),
    ):
        if lightdock_atom.name != prody_atom.getName():
            raise NormalModesCalculationError(
                "Atoms differ: %s - %s" % (str(lightdock_atom), str(prody_atom))
            )

    molecule_anm_ext, molecule_all = extendModel(
        molecule_anm, backbone_atoms, prody_molecule, norm=True
    )

    calculated_n_modes = len(molecule_anm_ext)
    variances = molecule_anm_ext.getVariances()
    magnitudes = np.array([abs(mode) for mode in molecule_anm_ext])

    if np.any(variances == 0):
        raise NormalModesCalculationError("One or more modes has zero variance")

    np.random.seed(seed)
    randn = np.random.normal(DEFAULT_EXTENT_MU, DEFAULT_EXTENT_SIGMA,
        size=(DEFAULT_NUM_SWARMS * DEFAULT_NUM_GLOWWORMS, n_modes))
    coef = ((randn ** 2 * variances).sum(1) ** 0.5).mean()
    scale = num_atoms_prody ** 0.5 * rmsd / coef
    scale = scale / magnitudes * variances ** 0.5

    modes = []
    try:
        for i in range(calculated_n_modes):
            nm = molecule_anm_ext.getEigvecs()[:, i].reshape((num_atoms_prody, 3)) * scale[i]
            modes.append(nm)
    except (ValueError, IndexError) as err:
        log.info(
            "Number of atoms of the ANM model: %s" % str(molecule_anm_ext.numAtoms())
        )
        log.info(
            "Number of nodes in the model: %s"
            % str((molecule_anm_ext.getEigvecs()).shape)
        )
        print(err)
        raise NormalModesCalculationError(
            "Number of atoms and ANM model differ. Please, check there are no missing "
            "nucleotides nor residues."
        )
    if calculated_n_modes < n_modes:
        log.warning(
            "Number of non-trivial calculated modes is %d (asked for %d)"
            % (calculated_n_modes, n_modes)
        )
        # Padding
        for i in range(n_modes - calculated_n_modes):
            modes.append(np.zeros((num_atoms_prody, 3)))

    return np.array(modes)


def write_nmodes(n_modes, output_file_name):
    """Writes the previous calculated n_modes to a binary numpy file"""
    np.save(output_file_name, n_modes)


def read_nmodes(file_name):
    """Reads normal modes from a numpy binary file"""
    return np.load(file_name)
