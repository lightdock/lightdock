"""Tests for DFIRE scoring function module"""

import pytest
from pathlib import Path
from lightdock.scoring.dfire.driver import DFIREPotential, DFIRE, DFIREAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.structure.atom import Atom
from lightdock.structure.residue import Residue
from lightdock.structure.chain import Chain
from lightdock.scoring.dfire.cython.cdfire import calculate_dfire
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF


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
        "distance, in_contact",
        [(2.0, True), (3.0, True), (3.89, True), (3.91, False), (5.0, False)],
    )
    def test_interface_cutoff(self, distance, in_contact):
        def single_atom_complex(chain_id):
            atom = Atom(1, "CA", chain_id=chain_id, residue_name="ALA", residue_number=1)
            residue = Residue("ALA", 1, atoms=[atom])
            return Complex([Chain(chain_id, [residue])])

        adapter = DFIREAdapter(single_atom_complex("A"), single_atom_complex("B"))
        receptor, ligand = adapter.receptor_model, adapter.ligand_model
        pose = ligand.coordinates[0].clone()
        pose.translate([distance, 0.0, 0.0])
        scoring = self.dfire
        _, interface_receptor, interface_ligand = calculate_dfire(
            receptor,
            receptor.coordinates[0],
            ligand,
            pose,
            scoring.potential.dfire_dist_to_bins,
            scoring.potential.dfire_energy,
            interface_cutoff=DEFAULT_CONTACT_RESTRAINTS_CUTOFF,
        )

        expected = {0} if in_contact else set()
        assert set(interface_receptor) == expected
        assert set(interface_ligand) == expected
