"""calculate_dfire2 Cython implementation"""


def calculate_dfire2(res_index, atom_index, coordinates, potentials, mol_length: int, interface_cutoff: float) -> tuple[float, set[int], set[int]]:
    """
    calculate_dfire2 DFIRE2 Cython implementation.
    
    Returns
    -------
    tuple[energy,interface_receptor,interface_ligand]
    """
