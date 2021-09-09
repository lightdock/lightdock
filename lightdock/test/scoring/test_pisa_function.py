"""Tests for PISA scoring function module"""

from pathlib import Path
from nose.tools import assert_almost_equal
from lightdock.scoring.pisa.driver import PISAPotential, PISA, PISAAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestPISAPotential:
    def test_create_PISAPotential_interface(self):
        potential = PISAPotential()
        assert potential is not None


class TestPISA:
    """Original PISA scoring energy goes from negative to positive"""

    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        self.pisa = PISA()

    def test_calculate_PISA_1PPE(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPErec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPElig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = PISAAdapter(receptor, ligand)
        assert_almost_equal(
            -0.4346,
            round(
                self.pisa(
                    adapter.receptor_model,
                    adapter.receptor_model.coordinates[0],
                    adapter.ligand_model,
                    adapter.ligand_model.coordinates[0],
                ),
                4,
            ),
        )

    def test_calculate_PISA_1EAW(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWrec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWlig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = PISAAdapter(receptor, ligand)
        assert_almost_equal(
            -0.2097,
            round(
                self.pisa(
                    adapter.receptor_model,
                    adapter.receptor_model.coordinates[0],
                    adapter.ligand_model,
                    adapter.ligand_model.coordinates[0],
                ),
                4,
            ),
        )

    def test_calculate_PISA_1AY7(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7rec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7lig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = PISAAdapter(receptor, ligand)
        assert_almost_equal(
            -0.2141,
            round(
                self.pisa(
                    adapter.receptor_model,
                    adapter.receptor_model.coordinates[0],
                    adapter.ligand_model,
                    adapter.ligand_model.coordinates[0],
                ),
                4,
            ),
        )
