"""calculate_energy pyDock C implementation"""

from lightdock.structure.space import SpacePoints


def calculate_energy(receptor_coordinates: SpacePoints, ligand_coordinates: SpacePoints, rec_charges, lig_charges, rec_vdw, lig_vdw, rec_vdw_radii, lig_vdw_radii, rec_hydrogens, lig_hydrogens, rec_asa, lig_asa, rec_des_energy, lig_des_energy, interface_cutoff: float):
    """
    calculate_energy pyDock C implementation.

    Returns
    -------
    tuple[energy,vdw,solv_rec,solv_lig,interface_receptor,interface_ligand]
    """
    ...
