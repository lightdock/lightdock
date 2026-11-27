"""Tests for DFIRE scoring function module"""

import pytest
from pathlib import Path
from lightdock.scoring.dfire.driver import DFIREPotential, DFIRE, DFIREAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


class TestDFIREPotential:
    def test_create_DFIREPotential_interface(self):
        potential = DFIREPotential()
        # Check if there is energy for the 20 residues
        assert len(potential.dfire_energy) == 20


class TestDFIRE:
    def setup_class(self):
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
        assert -17.3749982699 == pytest.approx(
            self.dfire(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
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
        assert -16.2182794457 == pytest.approx(
            self.dfire(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
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
        assert -20.7486309727 == pytest.approx(
            self.dfire(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    @pytest.mark.parametrize(
        "dist,expected_contact",
        [
            (0.25, True),
            (2.4, True),
            (2.5, True),
            (3.0, True),
            (3.89, True),
            (3.9, True),
            (3.91, False),
        ],
    )
    def test_interface_cutoff_uses_distance_in_angstroms(self, dist, expected_contact):
        from lightdock.structure.atom import Atom
        from lightdock.structure.residue import Residue
        from lightdock.structure.chain import Chain
        from lightdock.scoring.dfire.cython.cdfire import calculate_dfire

        a1 = Atom(
            atom_number=1,
            atom_name="CA",
            chain_id="A",
            residue_name="ALA",
            residue_number=1,
            x=0.0,
            y=0.0,
            z=0.0,
        )
        a2 = Atom(
            atom_number=1,
            atom_name="CA",
            chain_id="B",
            residue_name="ALA",
            residue_number=1,
            x=dist,
            y=0.0,
            z=0.0,
        )
        receptor = Complex([Chain("A", [Residue("ALA", 1, atoms=[a1])])], [a1])
        ligand = Complex([Chain("B", [Residue("ALA", 1, atoms=[a2])])], [a2])
        adapter = DFIREAdapter(
            receptor,
            ligand,
            receptor_restraints=["A.ALA.1"],
            ligand_restraints=["B.ALA.1"],
        )

        energy, intf_rec, intf_lig = calculate_dfire(
            adapter.receptor_model,
            adapter.receptor_model.coordinates[0],
            adapter.ligand_model,
            adapter.ligand_model.coordinates[0],
            self.dfire.potential.dfire_dist_to_bins,
            self.dfire.potential.dfire_energy,
            interface_cutoff=3.9,
        )
        expected = {0} if expected_contact else set()
        assert intf_rec == expected
        assert intf_lig == expected

        score = self.dfire(
            adapter.receptor_model,
            adapter.receptor_model.coordinates[0],
            adapter.ligand_model,
            adapter.ligand_model.coordinates[0],
        )
        expected_score = energy * 3.0 if expected_contact else energy
        assert score == pytest.approx(expected_score)
