"""Tests for TOBIBAHAR scoring function module"""

import pytest
from pathlib import Path
from lightdock.scoring.tobibahar.driver import TOBIBAHARPotential, TOBIBAHAR, TOBIBAHARAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestTOBIBAHARPotential:
    def test_create_TOBIBAHARPotential_interface(self):
        potential = TOBIBAHARPotential()

        assert len(potential.tobibahar) == 22

        assert -3.56 == pytest.approx(potential.tobibahar[0][0])
        assert 1.82 == pytest.approx( potential.tobibahar[-1][-1])
        assert 1.6 == pytest.approx(potential.tobibahar[1][20])


class TestTOBIBAHAR:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        self.tobisc = TOBIBAHAR()

    def test_calculate_TOBIBAHAR_1PPE(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPErec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPElig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBIBAHARAdapter(receptor, ligand)
        assert 27.67 == pytest.approx(
            self.tobisc(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_TOBIBAHAR_1EAW(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWrec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWlig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBIBAHARAdapter(receptor, ligand)
        assert -14.64 == pytest.approx(
            self.tobisc(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_TOBIBAHAR_1AY7(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7rec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7lig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBIBAHARAdapter(receptor, ligand)
        assert -65.94 == pytest.approx(
            self.tobisc(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_TOBIBAHAR_1CZY(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1czy_protein.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1czy_peptide.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBIBAHARAdapter(receptor, ligand)
        assert 14.70 == pytest.approx(
            self.tobisc(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )