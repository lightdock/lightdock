"""calculate_dfire C implementation"""

from lightdock.structure.model import DockingModel
from lightdock.structure.space import SpacePoints


def calculate_dfire(receptor, ligand, dfire_energy, receptor_coordinates: SpacePoints, ligand_coordinates: SpacePoints, interface_cutoff: float):
    """
    calculate_dfire C implementation
    
    Returns
    -------
    tuple[energy,interface_receptor,interface_ligand]
    """
    ...
