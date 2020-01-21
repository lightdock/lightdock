"""Tests for poses related module"""

import os
import shutil
import filecmp
import numpy as np
from nose.tools import assert_almost_equal
from lightdock.prep.poses import normalize_vector, quaternion_from_vectors, get_quaternion_for_restraint, \
        get_random_point_within_sphere, estimate_membrane, upper_layer, populate_poses, calculate_initial_poses
from lightdock.mathutil.cython.quaternion import Quaternion
from lightdock.structure.residue import Residue
from lightdock.mathutil.lrandom import MTGenerator, NMExtentGenerator
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.prep.simulation import parse_restraints_file, get_restraints


class TestPoses:

    def setUp(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_poses/'
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
        os.mkdir(self.test_path)
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'

    def tearDown(self):
        try:
            shutil.rmtree(self.test_path)
        except:
            pass
    
    def test_normalize_vector1(self):
        v = np.array([1., 0., 0.])
        e = np.array([1., 0., 0.])
        n = normalize_vector(v)
        assert np.allclose(n, e)

    def test_normalize_vector2(self):
        v = np.array([4., 0., 0.])
        e = np.array([1., 0., 0.])
        n = normalize_vector(v)
        assert np.allclose(n, e)

    def test_normalize_vector3(self):
        v = np.array([0., 0., 0.])
        e = np.array([0., 0., 0.])
        n = normalize_vector(v)
        assert np.allclose(n, e)

    def test_normalize_vector4(self):
        v = np.array([2., -2., 0.])
        e = np.array([0.70710678, -0.70710678, 0.])
        n = normalize_vector(v)
        assert np.allclose(n, e)

    def test_quaternion_from_vectors1(self):
        a = np.array([1., 0., 0.])
        b = np.array([1., 0., 0.])
        q = quaternion_from_vectors(a, b)
        e = Quaternion()
        assert e == q

    def test_quaternion_from_vectors2(self):
        a = np.array([1., 0., 0.])
        b = np.array([2., 0., 0.])
        q = quaternion_from_vectors(a, b)
        e = Quaternion()
        assert e == q

    def test_quaternion_from_vectors3(self):
        a = np.array([1., 0., 0.])
        b = np.array([0., 2., 0.])
        q = quaternion_from_vectors(a, b)
        # 90 degrees in Z axis
        e = Quaternion(w=0.70710678, x=0., y=0., z=0.70710678)
        assert e == q

    def test_quaternion_from_vectors4(self):
        a = np.array([1., 0., 0.])
        b = np.array([-1., 0., 0.])
        q = quaternion_from_vectors(a, b)
        # Orthogonal rotation, two possibilities:
        e1 = Quaternion(w=0., x=0., y=1.0, z=0.)
        e2 = Quaternion(w=0., x=0., y=-1.0, z=0.)
        assert e1 == q or e2 == q

    def test_get_quaternion_for_restraint1(self):
        # Origin is 0,0,0
        # Ligand center
        l_center = [-5., 5., -5.]
        # Ligand residue
        l_res = [-7., 7., -7.]
        # Receptor residue
        r_res = [3., 1., 0.]
        # Calculate quaternion for rotation from l_res to point to r_res
        rec_residue = Residue.dummy(r_res[0], r_res[1], r_res[2])
        lig_residue = Residue.dummy(l_res[0], l_res[1], l_res[2])

        q = get_quaternion_for_restraint(rec_residue, lig_residue,
                                         l_center[0], l_center[1], l_center[2], # translation
                                         [0., 0., 0.], # receptor translation 
                                         [5., -5., 5.] # ligand translation
                                         )

        e = Quaternion(w=0.14518697, x=0.19403814, y=-0.58211441, z=-0.77615254)

        assert e == q

    def test_get_quaternion_for_restraint2(self):
        # Origin is 0,0,0
        # Ligand center
        l_center = [-5., 5., -5.]
        # Ligand residue
        l_res = [-7., 6., -7.]
        # Receptor residue
        r_res = [3., 1., 0.]
        # Calculate quaternion for rotation from l_res to point to r_res
        rec_residue = Residue.dummy(r_res[0], r_res[1], r_res[2])
        lig_residue = Residue.dummy(l_res[0], l_res[1], l_res[2])

        q = get_quaternion_for_restraint(rec_residue, lig_residue,
                                         l_center[0], l_center[1], l_center[2], # translation
                                         [0., 0., 0.], # receptor translation 
                                         [5., -5., 5.] # ligand translation
                                         )

        e = Quaternion(0.10977233, -0.44451098, -0.88902195, 0.)

        assert e == q

    def test_get_quaternion_for_restraint2d(self):
        # Origin is 0,0,0
        # Ligand center
        l_center = [5., 5., 0.]
        # Ligand residue
        l_res = [7., 7., 0.]
        # Receptor residue
        r_res = [3., 1., 0.]
        # Calculate quaternion for rotation from l_res to point to r_res
        rec_residue = Residue.dummy(r_res[0], r_res[1], r_res[2])
        lig_residue = Residue.dummy(l_res[0], l_res[1], l_res[2])

        q = get_quaternion_for_restraint(rec_residue, lig_residue,
                                         l_center[0], l_center[1], l_center[2], # translation
                                         [0., 0., 0.], # receptor translation 
                                         [5., 5., 0.] # ligand translation
                                         )

        e = Quaternion(0.16018224, 0., 0., -0.98708746)

        assert e == q

    def test_get_quaternion_for_restraint2d_different_quadrant(self):
        # Origin is 0,0,0
        # Ligand center
        l_center = [5., -5., 0.]
        # Ligand residue
        l_res = [7., -7., 0.]
        # Receptor residue
        r_res = [-3., 1., 0.]
        # Calculate quaternion for rotation from l_res to point to r_res
        rec_residue = Residue.dummy(r_res[0], r_res[1], r_res[2])
        lig_residue = Residue.dummy(l_res[0], l_res[1], l_res[2])

        q = get_quaternion_for_restraint(rec_residue, lig_residue,
                                         l_center[0], l_center[1], l_center[2], # translation
                                         [0., 0., 0.], # receptor translation 
                                         [5., -5., 0.] # ligand translation
                                         )

        e = Quaternion(0.07088902, 0., 0., -0.99748421)

        assert e == q

    def test_get_random_point_within_sphere(self):
        seed = 1984
        rng = MTGenerator(seed)
        to_generate = 10
        radius = 10.0

        points = [get_random_point_within_sphere(rng, radius) for _ in range(to_generate)]

        # Check all generated points are within a sphere of a given radius and
        # centered at the origin
        for p in points:
            assert p[0] <= radius and p[0] >= -radius
            assert p[1] <= radius and p[1] >= -radius
            assert p[2] <= radius and p[2] >= -radius

    def test_estimate_membrane(self):
        seed = 1984
        np.random.seed(seed)
        z_coordinates_two_layers = np.random.rand(10, 3) * 30
        z_coordinates_two_layers[::2][:,2] *= -1

        layers = estimate_membrane(z_coordinates_two_layers[:,2], 15.)

        assert len(layers) == 2

        z_coordinates_one_layer = np.random.rand(10, 3) * 30

        layers = estimate_membrane(z_coordinates_one_layer[:,2])

        assert len(layers) == 1

    def test_upper_layer(self):
        seed = 1984
        np.random.seed(seed)
        z_coordinates_two_layers = np.random.rand(10, 3) * 30
        z_coordinates_two_layers[::2][:,2] *= -1
        expected = np.array([14.865534412778462, 19.42594208380253, 21.84073297851065, 23.828220097290853, 23.927022857098866])
        
        layers = estimate_membrane(z_coordinates_two_layers[:,2], 15.)

        upper = upper_layer(layers)
        
        assert np.allclose(upper, expected)


    def test_populate_poses(self):
        to_generate = 3
        center = [0., 0., 0.]
        radius = 10.
        seed = 1984
        number_generator = MTGenerator(seed)
        # No needed if no restraints are provided
        rec_translation = [0., 0., 0.]
        lig_translation = [-15.0, -15.0, -15.0]
        
        poses = populate_poses(to_generate, center, radius, number_generator, rec_translation, lig_translation)

        expected = [[-2.7294328828938386, -0.11588636361606675, -3.2077982565639607, -0.6725323168859286, 0.5755106199757826, -0.37756439233549577, 0.27190612107759393], 
                    [-2.2842918231021314, -5.118850363927159, 3.2622525507180566, -0.6132005094145724, 0.757439137867041, -0.0367297456106423, -0.22118321244687617], 
                    [2.906625032282104, -3.6251691461446978, -6.096162381118752, -0.8727759595729266, -0.22555077047578467, -0.2572320124803007, 0.3481675833340481]]

        # We generate 3 poses
        assert len(poses) == 3
        # We generate poses without ANM
        assert len(poses[0]) == 7
        # We generate the expected poses
        assert np.allclose(expected, poses) 

    def test_populate_poses_both_restraints(self):
        to_generate = 3
        center = [15., 15., 15.]
        radius = 10.
        seed = 1984
        number_generator = MTGenerator(seed)
        # No needed if no restraints are provided
        rec_translation = [0., 0., 0.]
        lig_translation = [-15.0, -15.0, -15.0]
        # Restraints
        receptor_restraints = [Residue.dummy(1.0, 1.0, 1.0)]
        ligand_restraints = [Residue.dummy(16.0, 16.0, 16.0)]
        ligand_diameter = 10.

        poses = populate_poses(to_generate, center, radius, number_generator, rec_translation, lig_translation,
                               receptor_restraints=receptor_restraints, ligand_restraints=ligand_restraints,
                               ligand_diameter=ligand_diameter)

        expected = [[12.270567117106161, 14.884113636383933, 11.792201743436038, 0.05643308208136652, 0.7572287178886188, -0.11715469622945059, -0.6400740216591683], 
                    [21.986658842426436, 12.715708176897868, 9.881149636072841, 0.17754420770338322, 0.179865123253213, -0.7681474466249519, 0.5882823233717389], 
                    [22.833743558777538, 15.523806353077699, 17.906625032282104, 0.0847943426151146, -0.2599970108635702, -0.5376137514110186, 0.7976107622745888]]

        # We generate 3 poses
        assert len(poses) == 3
        # We generate poses with ANM
        assert len(poses[0]) == 7
        # We generate the expected poses
        assert np.allclose(expected, poses)

    def test_populate_poses_restraints_only_receptor(self):
        to_generate = 3
        center = [15., 15., 15.]
        radius = 10.
        seed = 1984
        number_generator = MTGenerator(seed)
        # No needed if no restraints are provided
        rec_translation = [0., 0., 0.]
        lig_translation = [-15.0, -15.0, -15.0]
        # Restraints
        receptor_restraints = [Residue.dummy(1.0, 1.0, 1.0)]
        ligand_diameter = 10.

        poses = populate_poses(to_generate, center, radius, number_generator, rec_translation, lig_translation,
                               receptor_restraints=receptor_restraints, ligand_diameter=ligand_diameter)

        expected = [[12.270567117106161, 14.884113636383933, 11.792201743436038, -0.6725323168859286, 0.5755106199757826, -0.37756439233549577, 0.27190612107759393], 
                    [12.715708176897868, 9.881149636072841, 18.262252550718056, -0.6132005094145724, 0.757439137867041, -0.0367297456106423, -0.22118321244687617], 
                    [17.906625032282104, 11.374830853855302, 8.903837618881248, -0.8727759595729266, -0.22555077047578467, -0.2572320124803007, 0.3481675833340481]]

        # We generate 3 poses
        assert len(poses) == 3
        # We generate poses with ANM
        assert len(poses[0]) == 7
        # We generate the expected poses
        assert np.allclose(expected, poses)

    def test_populate_poses_restraints_only_ligand(self):
        to_generate = 3
        center = [15., 15., 15.]
        radius = 10.
        seed = 1984
        number_generator = MTGenerator(seed)
        # No needed if no restraints are provided
        rec_translation = [0., 0., 0.]
        lig_translation = [-15.0, -15.0, -15.0]
        # Restraints
        ligand_restraints = [Residue.dummy(16.0, 16.0, 16.0)]
        ligand_diameter = 30.

        poses = populate_poses(to_generate, center, radius, number_generator, rec_translation, lig_translation,
                               ligand_restraints=ligand_restraints, ligand_diameter=ligand_diameter)

        expected = [[12.270567117106161, 14.884113636383933, 11.792201743436038, 0.7092076459798842, 0.5346980891115, -0.08272585379354264, -0.4519722353179574],
                    [22.833743558777538, 15.523806353077699, 17.906625032282104, 0.24053986913227657, -0.25327548418133206, -0.5237151871050496, 0.7769906712863817],
                    [8.903837618881248, 8.747779486586737, 19.195006601282643, 0.7560980480558669, -0.46621474071247004, 0.45925053854810693, 0.00696420216436307]]

        # We generate 3 poses
        assert len(poses) == 3
        # We generate poses with ANM
        assert len(poses[0]) == 7
        # We generate the expected poses
        assert np.allclose(expected, poses)

    def test_calculate_initial_poses(self):

        atoms, residues, chains = parse_complex_from_file(os.path.join(self.golden_data_path, 
                                                                       '3p0g', 'receptor_membrane.pdb'))
        receptor = Complex(chains)
        atoms, residues, chains = parse_complex_from_file(os.path.join(self.golden_data_path, 
                                                                       '3p0g', 'ligand.pdb'))
        ligand = Complex(chains)
        
        rec_translation = receptor.move_to_origin()
        lig_translation = ligand.move_to_origin()

        num_swarms = 10
        num_glowworms = 10
        seed = 324324
        restraints = parse_restraints_file(os.path.join(self.golden_data_path, '3p0g', 'restraints.list'))
        receptor_restraints = get_restraints(receptor, restraints['receptor'])
        ligand_restraints = get_restraints(ligand, restraints['ligand'])
        rec_restraints = receptor_restraints['active'] + receptor_restraints['passive']
        lig_restraints = ligand_restraints['active'] + ligand_restraints['passive']

        positions_files = calculate_initial_poses(receptor, ligand, num_swarms, num_glowworms,
                                                  seed, rec_restraints, lig_restraints,
                                                  rec_translation, lig_translation,
                                                  dest_folder=self.test_path, is_membrane=True)

        assert filecmp.cmp(positions_files[0], os.path.join(self.golden_data_path, '3p0g', 'init', 'initial_positions_0.dat'))
        assert filecmp.cmp(positions_files[1], os.path.join(self.golden_data_path, '3p0g', 'init', 'initial_positions_1.dat'))
