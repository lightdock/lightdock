"""Test module for rotamer predictor and interface"""

import os
import math
from nose.tools import assert_almost_equals
from lightdock.scoring.dfire.driver import DFIRE, DFIREAdapter
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.rotamer.predictor import get_interface_residues, calculate_chi_angles, calculate_phi_psi


class TestPredictor:
    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'
        self.dfire = DFIRE()

    def tearDown(self):
        pass

    def test_calculate_interface(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPErec.pdb')
        receptor = Complex(chains, atoms)
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        ligand = Complex(chains, atoms)
        adapter = DFIREAdapter(receptor, ligand)
        receptor_residue_indexes, ligand_residue_indexes, dist_matrix, atom_indexes = get_interface_residues(
            adapter.receptor_model,
            adapter.ligand_model,
            max_cutoff=1.)
        receptor_residues = []
        for residue_index in receptor_residue_indexes:
            residue = receptor.residues[residue_index]
            receptor_residues.append("%s.%d" % (residue.name, residue.number))
        ligand_residues = []
        for residue_index in ligand_residue_indexes:
            residue = ligand.residues[residue_index]
            ligand_residues.append("%s.%d" % (residue.name, residue.number))

        # Checked according to UCSF chimera
        expected_receptor = ['SER.61', 'ALA.56', 'TYR.59', 'TYR.94', 'GLY.62', 'ILE.63', 'MET.104', 'ILE.106']
        expected_ligand = ['LEU.7', 'LEU.8', 'GLU.24', 'HIS.25', 'GLY.26', 'TYR.27', 'CYS.28']

        assert len(set(expected_receptor) & set(receptor_residues)) == 8
        assert len(set(expected_ligand) & set(ligand_residues)) == 7

    def test_calculate_chi_angles(self):
        # Command in chimera: angle #0:28@n,ca,cb,sg
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        protein = Complex(chains, atoms)
        for residue in protein.residues:
            if residue.name == 'ARG' and residue.number == 1:
                chi_angles = calculate_chi_angles(residue)
                for key, value in chi_angles.items():
                    if value:
                        if key == 'x1':
                            # angle #0:1@n,ca,cb,cg
                            assert_almost_equals(-173.216, round(math.degrees(value), 3))
                        if key == 'x2':
                            # angle #0:1@ca,cb,cg,cd
                            assert_almost_equals(171.375, round(math.degrees(value), 3))
                        if key == 'x3':
                            # angle #0:1@cb,cg,cd,ne
                            assert_almost_equals(62.852, round(math.degrees(value), 3))
                        if key == 'x4':
                            # angle #0:1@cg,cd,ne,cz
                            assert_almost_equals(100.22, round(math.degrees(value), 2))
            if residue.name == 'VAL' and residue.number == 2:
                # angle #0:2@n,ca,cb,cg1
                chi_angles = calculate_chi_angles(residue)
                for key, value in chi_angles.items():
                    if value:
                        if key == 'x1':
                            # angle #0:1@n,ca,cb,cg
                            assert_almost_equals(-173.629, round(math.degrees(value), 3))
            if residue.name == 'CYS' and residue.number == 28:
                chi_angles = calculate_chi_angles(residue)
                for key, value in chi_angles.items():
                    if value:
                        if key == 'x1':
                            assert_almost_equals(-60.2554, round(math.degrees(value), 4))

    def test_calculate_phi_psi(self):
        atoms, residues, chains = parse_complex_from_file(self.golden_data_path + '1PPElig.pdb')
        protein = Complex(chains, atoms)
        # psi: angle #0:1@n,ca,c #0:2@n
        # phi: angle #0:1@c #0:2@n,ca,c

        phi_angles = [-105.92428251619579, -134.402235889132, -58.32268858533758, -85.62997439535678, -129.666484600813,
                      -77.00076813772478, -142.09891098624075, -82.10672119029674, -163.14606891327375,
                      -109.37900096123484, -138.72905680654182, -59.09699793329797, -60.06774387010816,
                      -74.41030551527874, -99.82766540256617, -92.6495110068149, 54.969041241310705,
                      -104.60151419194615, -67.57074855137641, -123.83574594954692, -85.90313254423194,
                      -87.7781803331676, -66.345484249271, -64.51513795752882, 108.23656098935888, -129.62530277139578,
                      -71.90658189461674, -170.4460918036806]
        psi_angles = [138.38576328505278, 105.22472788100255, 106.42882930892199, 150.65572151747787, 72.08329638522976,
                      130.19890858175336, 115.48238807519739, 132.48041144914038, 163.35191386073618,
                      151.17756189538443, -28.310569696143393, 162.66293554938997, -32.25480696024475,
                      -20.28436719199857, -11.444789534534305, 163.38578466073147, 150.2534549328882,
                      -128.53524744082424, 20.01260634937939, 151.96710290169335, 159.55519588393594,
                      115.07091589216549, 152.8911959270869, -24.04765297807205, -14.890186424782046, 15.86273088398991,
                      152.7552784042674, 146.11762131430552]

        for i in range(1, len(protein.residues)):
            phi, psi = calculate_phi_psi(protein.residues[i], protein.residues[i - 1])
            assert_almost_equals(phi_angles[i - 1], math.degrees(phi))
            assert_almost_equals(psi_angles[i - 1], math.degrees(psi))
