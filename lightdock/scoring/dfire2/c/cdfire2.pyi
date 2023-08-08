"""calculate_dfire2 C implementation"""


def calculate_dfire2(res_index, atom_index, coordinates, dfire2_energy, mol_length: int, interface_cutoff: float) -> tuple[float, set[int], set[int]]:
    """
    calculate_dfire2 C implementation.
    
    Returns
    -------
    tuple[energy,interface_receptor,interface_ligand]
    """
