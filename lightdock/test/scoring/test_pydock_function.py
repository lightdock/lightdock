"""Tests for CPyDock scoring function module"""

import pytest
from pathlib import Path
from lightdock.scoring.cpydock.driver import CPyDock, CPyDockAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestPyDock:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        self.pydock = CPyDock()

    def test_calculate_PyDock_1AY7(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7_rec.pdb"
        )
        receptor = Complex(
            chains,
            atoms,
            structure_file_name=(self.golden_data_path / "1AY7_rec.pdb"),
        )
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7_lig.pdb"
        )
        ligand = Complex(
            chains,
            atoms,
            structure_file_name=(self.golden_data_path / "1AY7_lig.pdb"),
        )
        adapter = CPyDockAdapter(receptor, ligand)
        assert -15.923994756 == pytest.approx(
            self.pydock(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )
