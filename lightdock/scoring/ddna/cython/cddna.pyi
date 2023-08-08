from lightdock.structure.space import SpacePoints


def calculate_ddna(receptor, receptor_coordinates: SpacePoints, ligand, ligand_coordinates: SpacePoints, ddna_potentials: list[float], ddna_map: list[int], interface_cutoff: float = 3.9) -> tuple[float, set[int], set[int]]:
    """
    calculate_ddna Cython implementation.

    Returns
    -------
    tuple[energy,interface_receptor,interface_ligand]
    """
