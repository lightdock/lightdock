"""Tests for DDNA scoring function module"""

import pytest
from pathlib import Path
from lightdock.scoring.ddna.driver import DDNA, DDNAAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.structure.atom import Atom
from lightdock.structure.residue import Residue
from lightdock.structure.chain import Chain
from lightdock.scoring.ddna.cython.cddna import calculate_ddna
from lightdock.constants import DEFAULT_CONTACT_RESTRAINTS_CUTOFF


class TestDDNA:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"
        self.dna = DDNA()

    def test_calculate_DNA_1AZP(self):
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1azp_prot.pdb"
        )
        receptor = Complex(
            chains, atoms, structure_file_name=(self.golden_data_path / "1azp_prot.pdb")
        )
        atoms, _, chains = parse_complex_from_file(
            self.golden_data_path / "1azp_dna.pdb"
        )
        ligand = Complex(
            chains, atoms, structure_file_name=(self.golden_data_path / "1azp_dna.pdb")
        )
        adapter = DDNAAdapter(receptor, ligand)
        assert 6.915295143021656 == pytest.approx(
            self.dna(
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

    @pytest.mark.parametrize(
        "distance, in_contact",
        [(1.0, True), (3.0, True), (3.89, True), (3.91, False), (5.0, False)],
    )
    def test_interface_cutoff(self, distance, in_contact):
        adapter = DDNAAdapter(self._complex("A"), self._complex("B"))
        receptor, ligand = adapter.receptor_model, adapter.ligand_model
        pose = ligand.coordinates[0].clone()
        pose.translate([distance, 0.0, 0.0])
        _, interface_receptor, interface_ligand = calculate_ddna(
            receptor,
            receptor.coordinates[0],
            ligand,
            pose,
            self.dna.potential.potentials,
            self.dna.potential.map,
            interface_cutoff=DEFAULT_CONTACT_RESTRAINTS_CUTOFF,
        )

        expected = {0} if in_contact else set()
        assert set(interface_receptor) == expected
        assert set(interface_ligand) == expected
