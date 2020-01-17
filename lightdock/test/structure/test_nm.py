"""Tests for ANM calculus related module"""

import os
import shutil
import filecmp
import numpy as np
from nose.tools import assert_almost_equal
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.structure.nm import calculate_nmodes, write_nmodes, read_nmodes


class TestNM:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = os.path.join(self.path, 'scratch_nm')
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.join(
                                    os.path.normpath(os.path.dirname(os.path.realpath(__file__))),
                                    'golden_data')

    def tearDown(self):
        try:
            shutil.rmtree(self.test_path)
        except:
            pass

    def test_calculate_anm_protein_1(self):
        pdb_file = os.path.join(self.golden_data_path, 'nm_prot', '2UUY_lig.pdb')
        atoms, residues, chains = parse_complex_from_file(pdb_file)
        molecule = Complex(chains)

        nmodes = calculate_nmodes(pdb_file, n_modes=10, molecule=molecule)

        expected_nmodes = read_nmodes(os.path.join(self.golden_data_path, 'nm_prot', 'lightdock_lig.nm.npy'))

        assert np.allclose(expected_nmodes, nmodes)

    def test_calculate_anm_protein_2(self):
        pdb_file = os.path.join(self.golden_data_path, 'nm_prot', '2UUY_rec.pdb')
        atoms, residues, chains = parse_complex_from_file(pdb_file)
        molecule = Complex(chains)

        nmodes = calculate_nmodes(pdb_file, n_modes=10, molecule=molecule)

        expected_nmodes = read_nmodes(os.path.join(self.golden_data_path, 'nm_prot', 'lightdock_rec.nm.npy'))

        assert np.allclose(expected_nmodes, nmodes)

    def test_calculate_anm_dna(self):
        pdb_file = os.path.join(self.golden_data_path, 'nm_dna', '1DIZ_lig.pdb.H')
        atoms, residues, chains = parse_complex_from_file(pdb_file)
        molecule = Complex(chains)

        nmodes = calculate_nmodes(pdb_file, n_modes=10, molecule=molecule)

        expected_nmodes = read_nmodes(os.path.join(self.golden_data_path, 'nm_dna', 'lightdock_lig.nm.npy'))

        assert np.allclose(expected_nmodes, nmodes)

    def test_read_write(self):
        pdb_file = os.path.join(self.golden_data_path, 'nm_dna', '1DIZ_lig.pdb.H')
        atoms, residues, chains = parse_complex_from_file(pdb_file)
        molecule = Complex(chains)

        nmodes = calculate_nmodes(pdb_file, n_modes=10, molecule=molecule)
        write_nmodes(nmodes, os.path.join(self.test_path, 'test_nm'))

        expected_nmodes = read_nmodes(os.path.join(self.golden_data_path, 'nm_dna', 'lightdock_lig.nm.npy'))
        other_nmodes = read_nmodes(os.path.join(self.test_path, 'test_nm.npy'))

        assert np.allclose(expected_nmodes, other_nmodes)
