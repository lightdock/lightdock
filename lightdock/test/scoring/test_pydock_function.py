"""Tests for CPyDock scoring function module"""

from pathlib import Path
from nose.tools import assert_almost_equal
from lightdock.scoring.cpydock.driver import CPyDock, CPyDockAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestPyDock:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        self.pydock = CPyDock()

    def test_calculate_PyDock_1AY7(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7_rec.pdb.H"
        )
        receptor = Complex(
            chains,
            atoms,
            structure_file_name=(self.golden_data_path / "1AY7_rec.pdb.H"),
        )
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7_lig.pdb.H"
        )
        ligand = Complex(
            chains,
            atoms,
            structure_file_name=(self.golden_data_path / "1AY7_lig.pdb.H"),
        )
        adapter = CPyDockAdapter(receptor, ligand)
        assert_almost_equal(
            -15.923994756,
            self.pydock(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            ),
        )
