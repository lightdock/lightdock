"""Tests for C implementation of DFIRE scoring function module"""

import pytest
from pathlib import Path
from lightdock.scoring.fastdfire.driver import DFIRE, DFIREAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.structure.atom import Atom
from lightdock.structure.residue import Residue
from lightdock.structure.chain import Chain
from lightdock.scoring.fastdfire.c.cdfire import calculate_dfire
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF


class TestFastDFIRE:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        self.dfire = DFIRE()

    def test_calculate_FastDFIRE_1PPE(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPErec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPElig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = DFIREAdapter(receptor, ligand)
        assert -17.3745706065 == pytest.approx(
            self.dfire(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_FastDFIRE_1EAW(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWrec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWlig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = DFIREAdapter(receptor, ligand)
        assert -16.2239702546 == pytest.approx(
            self.dfire(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_FastDFIRE_1AY7(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7rec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7lig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = DFIREAdapter(receptor, ligand)
        assert -20.7459619159 == pytest.approx(
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
            ligand,
            scoring.potential.dfire_energy,
            receptor.coordinates[0],
            pose,
            DEFAULT_CONTACT_RESTRAINTS_CUTOFF,
        )

        expected = {0} if in_contact else set()
        assert set(interface_receptor) == expected
        assert set(interface_ligand) == expected
