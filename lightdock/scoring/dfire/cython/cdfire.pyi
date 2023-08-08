from lightdock.structure.space import SpacePoints
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF


def calculate_dfire(receptor, receptor_coordinates: SpacePoints, ligand, ligand_coordinates: SpacePoints,
                    dfire_dist_to_bins: list[int], dfire_energy: list[list[list[list[list[float]]]]], interface_cutoff: float = DEFAULT_CONTACT_RESTRAINTS_CUTOFF) -> tuple[float, set[int], set[int]]:
    """
    calculate_dfire Cython implementation.

    Returns
    -------
    tuple[energy,interface_receptor,interface_receptor]
    """
