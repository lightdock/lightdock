"""Tests for ANM calculus related module"""

import pytest
import numpy as np
from pathlib import Path
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.structure.nm import calculate_nmodes, write_nmodes, read_nmodes
from lightdock.constants import STARTING_NM_SEED, DEFAULT_ANM_RMSD
from lightdock.error.lightdock_errors import NormalModesCalculationError


class TestNM:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"

    def test_calculate_anm_protein_1(self):
        pdb_file = self.golden_data_path / "nm_prot" / "2UUY_lig.pdb"
        _, _, chains = parse_complex_from_file(pdb_file)
        molecule = Complex(chains)

        nmodes = calculate_nmodes(
            pdb_file,
            n_modes=10,
            rmsd=DEFAULT_ANM_RMSD,
            seed=STARTING_NM_SEED,
            molecule=molecule
        )

        expected_nmodes = read_nmodes(
            self.golden_data_path / "nm_prot" / "lightdock_lig.nm.npy"
        )

        assert np.allclose(expected_nmodes, nmodes)

    def test_calculate_anm_protein_2(self):
        pdb_file = self.golden_data_path / "nm_prot" / "2UUY_rec.pdb"
        _, _, chains = parse_complex_from_file(pdb_file)
        molecule = Complex(chains)

        nmodes = calculate_nmodes(
            pdb_file,
            n_modes=10,
            rmsd=DEFAULT_ANM_RMSD,
            seed=STARTING_NM_SEED,
            molecule=molecule
        )

        expected_nmodes = read_nmodes(
            self.golden_data_path / "nm_prot" / "lightdock_rec.nm.npy"
        )

        assert np.allclose(expected_nmodes, nmodes)

    def test_calculate_anm_dna(self):
        pdb_file = self.golden_data_path / "nm_dna" / "1DIZ_lig.pdb"
        _, _, chains = parse_complex_from_file(pdb_file)
        molecule = Complex(chains)

        nmodes = calculate_nmodes(
            pdb_file,
            n_modes=10,
            rmsd=DEFAULT_ANM_RMSD,
            seed=STARTING_NM_SEED,
            molecule=molecule
        )

        expected_nmodes = read_nmodes(
            self.golden_data_path / "nm_dna" / "lightdock_lig.nm.npy"
        )

        assert np.allclose(expected_nmodes, nmodes)

    def test_calculate_anm_wrong_extension(self):
        pdb_file = self.golden_data_path / "nm_dna" / "1DIZ_lig.pdb.H"
        _, _, chains = parse_complex_from_file(pdb_file)
        molecule = Complex(chains)

        # PATCH: assert_raises is not working properly with ProDy debug mode
        with pytest.raises(NormalModesCalculationError):
            nmodes = calculate_nmodes(
                pdb_file,
                n_modes=10,
                rmsd=DEFAULT_ANM_RMSD,
                seed=STARTING_NM_SEED,
                molecule=molecule
            )

    def test_read_write(self, tmp_path):
        pdb_file = self.golden_data_path / "nm_dna" / "1DIZ_lig.pdb"
        _, _, chains = parse_complex_from_file(pdb_file)
        molecule = Complex(chains)

        nmodes = calculate_nmodes(
            pdb_file,
            n_modes=10,
            rmsd=DEFAULT_ANM_RMSD,
            seed=STARTING_NM_SEED,
            molecule=molecule
        )
        write_nmodes(nmodes, tmp_path / "test_nm")

        expected_nmodes = read_nmodes(
            self.golden_data_path / "nm_dna" / "lightdock_lig.nm.npy"
        )
        other_nmodes = read_nmodes(tmp_path / "test_nm.npy")

        assert np.allclose(expected_nmodes, other_nmodes)
