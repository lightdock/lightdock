"""Tests for MJ3h potentials module"""

from pathlib import Path
from nose.tools import assert_almost_equal
from lightdock.scoring.mj3h.driver import MJPotential, MJ3h, MJ3hAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestMJPotential:
    def test_create_MJPotential_interface(self):
        mj = MJPotential()

        assert len(mj.potentials) == 4

    def test_create_MJ3h(self):
        mj3h = MJ3h()

        assert len(mj3h.potentials) == 20

        assert_almost_equal(-0.84, mj3h.potentials[0][0])
        assert_almost_equal(0.05, mj3h.potentials[15][3])
        assert_almost_equal(0.76, mj3h.potentials[19][19])


class TestMJ3hAdapter:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"

    def test_create_adapter(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPErec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPElig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(receptor, ligand)

        assert adapter.receptor_model.coordinates[0].coordinates.shape == (223, 3)
        assert adapter.ligand_model.coordinates[0].coordinates.shape == (29, 3)
        assert len(adapter.receptor_model.objects) == 223
        assert len(adapter.ligand_model.objects) == 29


class TestMJ3h:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        # Test for creating object in the interface tests for speeding up purposes
        self.mj3h = MJ3h()

    def test_calculate_MJ3h_1PPE(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPErec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPElig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(receptor, ligand)
        assert_almost_equal(
            2.02,
            self.mj3h(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            ),
        )
        # Original potential without min cutoff
        # assert_almost_equal(-17.94/2, self.mj3h(receptor, receptor.residue_coordinates, ligand, ligand.residue_coordinates))

    def test_calculate_MJ3h_1EAW(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWrec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWlig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(receptor, ligand)
        assert_almost_equal(
            -1.55,
            self.mj3h(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            ),
        )
        # Original potential without min cutoff
        # assert_almost_equal(-3.14/2, self.mj3h(receptor, receptor.residue_coordinates, ligand, ligand.residue_coordinates))

    def test_calculate_MJ3h_1AY7(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7rec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7lig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = MJ3hAdapter(receptor, ligand)
        assert_almost_equal(
            -12.92,
            self.mj3h(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            ),
        )
        # Original potential without min cutoff
        # assert_almost_equal(9.06/2, self.mj3h(receptor, receptor.residue_coordinates, ligand, ligand.residue_coordinates))
