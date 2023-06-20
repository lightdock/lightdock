"""calculate_pisa Cython implementation."""

from lightdock.structure.atom import Atom
from lightdock.structure.model import DockingModel
from lightdock.structure.space import SpacePoints


def calculate_pisa(receptor: DockingModel[Atom], receptor_coordinates: SpacePoints, ligand: DockingModel[Atom], ligand_coordinates: SpacePoints, pisa_energy: list[list[list[float]]], interface_cutoff: float) -> tuple[float, set[int], set[int]]:
    """
    calculate_pisa Cython implementation.

    Returns
    -------
    tuple[energy,interface_receptor,interface_ligand]
    """
    ...
