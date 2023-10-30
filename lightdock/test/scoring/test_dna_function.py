"""Tests for CPyDockDNA scoring function module"""

import pytest
from pathlib import Path
from lightdock.scoring.dna.driver import DNA, DNAAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestPyDockDNA:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        self.dna = DNA()

    def test_calculate_DNA_3MFK(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "3mfk_homodimer.pdb"
        )
        receptor = Complex(
            chains,
            atoms,
            structure_file_name=(self.golden_data_path / "3mfk_homodimer.pdb"),
        )
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "3mfk_dna.pdb"
        )
        ligand = Complex(
            chains, atoms, structure_file_name=(self.golden_data_path / "3mfk_dna.pdb")
        )
        adapter = DNAAdapter(receptor, ligand)
        assert -2716.68018700585 == pytest.approx(
            self.dna(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_DNA_3MFK_with_hydrogens(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "3mfk_homodimer_with_H.pdb"
        )
        receptor = Complex(
            chains,
            atoms,
            structure_file_name=(self.golden_data_path / "3mfk_homodimer_with_H.pdb"),
        )
        atoms, _t, chains = parse_complex_from_file(
            self.golden_data_path / "3mfk_dna.pdb"
        )
        ligand = Complex(
            chains, atoms, structure_file_name=(self.golden_data_path / "3mfk_dna.pdb")
        )
        adapter = DNAAdapter(receptor, ligand)
        assert 688.1703668834168 == pytest.approx(
            self.dna(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )
