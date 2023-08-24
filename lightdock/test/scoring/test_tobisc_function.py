"""Tests for TOBISC scoring function module"""

import pytest
from pathlib import Path
from lightdock.scoring.tobisc.driver import TOBISCPotential, TOBISC, TOBISCAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestTOBISCPotential:
    def test_create_TOBISCPotential_interface(self):
        potential = TOBISCPotential()

        assert len(potential.tobi_sc_1) == 22
        assert len(potential.tobi_sc_2) == 22

        assert -0.59 == pytest.approx(potential.tobi_sc_1[0][0])
        assert -0.09 == pytest.approx( potential.tobi_sc_1[-1][-1])
        assert 1.37 == pytest.approx(potential.tobi_sc_1[1][20])

        assert -0.58 == pytest.approx(potential.tobi_sc_2[0][0])
        assert -0.24 == pytest.approx(potential.tobi_sc_2[-1][-1])
        assert 0.39 == pytest.approx(potential.tobi_sc_2[3][20])


class TestTOBISC:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        self.tobisc = TOBISC()

    def test_calculate_TOBISC_1PPE(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPErec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPElig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBISCAdapter(receptor, ligand)
        assert 17.58 == pytest.approx(
            self.tobisc(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_TOBISC_1EAW(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWrec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWlig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBISCAdapter(receptor, ligand)
        assert -9.87 == pytest.approx(
            self.tobisc(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_TOBISC_1AY7(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7rec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7lig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBISCAdapter(receptor, ligand)
        assert 2.34 == pytest.approx(
            self.tobisc(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_TOBISC_1CZY(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1czy_protein.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1czy_peptide.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBISCAdapter(receptor, ligand)
        assert 16.28 == pytest.approx(
            self.tobisc(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )