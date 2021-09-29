"""Tests for DFIRE scoring function module"""

from pathlib import Path
from nose.tools import assert_almost_equal
from lightdock.scoring.dfire.driver import DFIREPotential, DFIRE, DFIREAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestDFIREPotential:
    def test_create_DFIREPotential_interface(self):
        potential = DFIREPotential()
        # Check if there is energy for the 20 residues
        assert len(potential.dfire_energy) == 20


class TestDFIRE:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        self.dfire = DFIRE()

    def test_calculate_DFIRE_1PPE(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPErec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPElig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = DFIREAdapter(receptor, ligand)
        assert_almost_equal(
            -17.3749982699,
            self.dfire(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            ),
        )

    def test_calculate_DFIRE_1EAW(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWrec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWlig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = DFIREAdapter(receptor, ligand)
        assert_almost_equal(
            -16.2182794457,
            self.dfire(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            ),
        )

    def test_calculate_DFIRE_1AY7(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7rec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7lig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = DFIREAdapter(receptor, ligand)
        assert_almost_equal(
            -20.7486309727,
            self.dfire(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            ),
        )
