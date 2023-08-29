"""Tests for TOBIA2 scoring function module"""

import pytest
from pathlib import Path
from lightdock.scoring.tobia2.driver import TOBIA2Potential, TOBIA2, TOBIA2Adapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestTOBIA2Potential:
    def test_create_TOBIA2Potential_interface(self):
        potential = TOBIA2Potential()

        assert len(potential.tobi_a_1) == 18
        assert len(potential.tobi_a_2) == 18

        assert -1.02 == pytest.approx(potential.tobi_a_1[0][0])
        assert 10.00 == pytest.approx(potential.tobi_a_1[-1][-1])
        assert 2.73 == pytest.approx(potential.tobi_a_1[1][17])

        assert 0.10 == pytest.approx(potential.tobi_a_2[0][0])
        assert 10.00 == pytest.approx(potential.tobi_a_2[-1][-1])
        assert 1.66 == pytest.approx(potential.tobi_a_2[3][17])


class TestTOBIA2:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        self.tobiA2 = TOBIA2()

    def test_calculate_TOBIA2_1PPE(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPErec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPElig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBIA2Adapter(receptor, ligand)
        assert -358.14 == pytest.approx(
            self.tobiA2(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_TOBIA2_1EAW(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWrec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWlig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBIA2Adapter(receptor, ligand)
        assert -280.54 == pytest.approx(
            self.tobiA2(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_TOBIA2_1AY7(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7rec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7lig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBIA2Adapter(receptor, ligand)
        assert -230.36 == pytest.approx(
            self.tobiA2(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_TOBIA2_1CZY(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1czy_protein.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1czy_peptide.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = TOBIA2Adapter(receptor, ligand)
        assert 43.23 == pytest.approx(
            self.tobiA2(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )