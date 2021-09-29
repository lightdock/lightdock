"""Tests for CPyDockDNA scoring function module"""

from pathlib import Path
from nose.tools import assert_almost_equal
from lightdock.scoring.dna.driver import DNA, DNAAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestPyDockDNA:
    def __init__(self):
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
        assert_almost_equal(
            -2716.68018700585,
            self.dna(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            ),
        )

    def test_calculate_DNA_3MFK_with_hydrogens(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "3mfk_homodimer.pdb.H"
        )
        receptor = Complex(
            chains,
            atoms,
            structure_file_name=(self.golden_data_path / "3mfk_homodimer.pdb.H"),
        )
        atoms, _t, chains = parse_complex_from_file(
            self.golden_data_path / "3mfk_dna.pdb"
        )
        ligand = Complex(
            chains, atoms, structure_file_name=(self.golden_data_path / "3mfk_dna.pdb")
        )
        adapter = DNAAdapter(receptor, ligand)
        assert_almost_equal(
            688.1703668834168,
            self.dna(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            ),
        )
