"""Tests for DFIRE2 scoring function module"""

import pytest
import numpy as np
from pathlib import Path
from lightdock.scoring.dfire2.driver import DFIRE2, DFIRE2Adapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.structure.atom import Atom
from lightdock.structure.residue import Residue
from lightdock.structure.chain import Chain
from lightdock.scoring.dfire2.c.cdfire2 import calculate_dfire2
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF


class TestDFIRE2:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        # FIXME: Segmentation fault when adapter is loaded here

    def test_calculate_DFIRE2_1PPE(self):
        dfire2 = DFIRE2()
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPErec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1PPElig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = DFIRE2Adapter(receptor, ligand)
        assert -398.7303561600074 == pytest.approx(
            dfire2(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_DFIRE2_1EAW(self):
        dfire2 = DFIRE2()
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWrec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1EAWlig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = DFIRE2Adapter(receptor, ligand)
        assert -488.34640492000244 == pytest.approx(
            dfire2(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    def test_calculate_DFIRE2_1AY7(self):
        dfire2 = DFIRE2()
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7rec.pdb"
        )
        receptor = Complex(chains, atoms)
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1AY7lig.pdb"
        )
        ligand = Complex(chains, atoms)
        adapter = DFIRE2Adapter(receptor, ligand)
        assert -283.19129030999665 == pytest.approx(
            dfire2(
                adapter.receptor_model,
                adapter.receptor_model.coordinates[0],
                adapter.ligand_model,
                adapter.ligand_model.coordinates[0],
            )
        )

    @staticmethod
    def _complex(chain_id, num_atoms=1):
        """Creates a complex of a single residue with all its atoms at the origin"""
        atoms = [
            Atom(i + 1, "CA", chain_id=chain_id, residue_name="ALA", residue_number=1)
            for i in range(num_atoms)
        ]
        return Complex([Chain(chain_id, [Residue("ALA", 1, atoms=atoms)])])

    def _interface(self, distance, rec_atoms=1, lig_atoms=1):
        """Interface of two molecules whose atoms are separated by distance"""
        adapter = DFIRE2Adapter(
            self._complex("A", rec_atoms), self._complex("B", lig_atoms)
        )
        receptor, ligand = adapter.receptor_model, adapter.ligand_model
        pose = ligand.coordinates[0].clone()
        pose.translate([distance, 0.0, 0.0])
        dfire2 = DFIRE2()
        # Initializes the cached residue and atom indexes
        dfire2(receptor, receptor.coordinates[0], ligand, pose)
        _, interface_receptor, interface_ligand = calculate_dfire2(
            dfire2.res_index,
            dfire2.atom_index,
            np.append(receptor.coordinates[0].coordinates, pose.coordinates).reshape(
                (-1, 3)
            ),
            dfire2.potential.energy,
            dfire2.molecule_length,
            len(receptor.objects),
            DEFAULT_CONTACT_RESTRAINTS_CUTOFF,
        )
        return interface_receptor, interface_ligand

    @pytest.mark.parametrize(
        "distance, in_contact",
        [(1.0, True), (3.0, True), (3.89, True), (3.91, False), (5.0, False)],
    )
    def test_interface_cutoff(self, distance, in_contact):
        interface_receptor, interface_ligand = self._interface(distance)

        expected = {0} if in_contact else set()
        assert set(interface_receptor) == expected
        assert set(interface_ligand) == expected

    def test_interface_uses_indexes_of_each_molecule(self):
        # Many more contacts than atoms, and atoms of the same molecule overlapping
        interface_receptor, interface_ligand = self._interface(3.0, 20, 30)

        assert len(interface_receptor) == len(interface_ligand) == 20 * 30
        assert set(interface_receptor) == set(range(20))
        assert set(interface_ligand) == set(range(30))
