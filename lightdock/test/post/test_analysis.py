"""Tests for analysis module"""

import pytest
import numpy as np
from pathlib import Path
from lightdock.post.analysis import (
    read_initial_positions_file,
    LightDockPose,
    read_predictions_file,
    LightDockPrediction,
)
from lightdock.error.lightdock_errors import GSOCoordinatesError


class TestInitialPositionsFiles:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / 'golden_data' / 'initial_positions_files'

    def test_read_initial_positions_file_good_no_anm(self):

        poses = read_initial_positions_file(self.golden_data_path / 'good_10_glowworms_no_anm.dat')

        assert len(poses) == 10

        first = poses[0]
        assert first.rec_extents.size == 0
        assert first.lig_extents.size == 0
        assert np.allclose([15.345706317, -14.691716422, 12.472000065], first.translation)
        assert -0.105840692 == pytest.approx(first.rotation.w)
        assert -0.484936911 == pytest.approx(first.rotation.x)
        assert 0.599743027 == pytest.approx(first.rotation.y)
        assert -0.627648183 == pytest.approx(first.rotation.z)

        last = poses[-1]
        assert last.rec_extents.size == 0
        assert last.lig_extents.size == 0
        assert np.allclose([17.549675220, -14.619564091, 27.076656519], last.translation)
        assert -0.235809453 == pytest.approx(last.rotation.w)
        assert 0.215334879 == pytest.approx(last.rotation.x)
        assert -0.938641307 == pytest.approx(last.rotation.y)
        assert 0.130296927 == pytest.approx(last.rotation.z)

    def test_read_initial_positions_file_good_no_anm_with_comment(self):

        poses = read_initial_positions_file(self.golden_data_path / 'good_10_glowworms_no_anm_comment.dat')

        assert len(poses) == 10

        first = poses[0]
        assert first.rec_extents.size == 0
        assert first.lig_extents.size == 0
        assert np.allclose([15.345706317, -14.691716422, 12.472000065], first.translation)
        assert -0.105840692 == pytest.approx(first.rotation.w)
        assert -0.484936911 == pytest.approx(first.rotation.x)
        assert 0.599743027 == pytest.approx(first.rotation.y)
        assert -0.627648183 == pytest.approx(first.rotation.z)

        last = poses[-1]
        assert last.rec_extents.size == 0
        assert last.lig_extents.size == 0
        assert np.allclose([17.549675220, -14.619564091, 27.076656519], last.translation)
        assert -0.235809453 == pytest.approx(last.rotation.w)
        assert 0.215334879 == pytest.approx(last.rotation.x)
        assert -0.938641307 == pytest.approx(last.rotation.y)
        assert 0.130296927 == pytest.approx(last.rotation.z)

    def test_read_initial_positions_file_malformed_nan(self):
        with pytest.raises(GSOCoordinatesError, match=r"NaN found in pose"):
            poses = read_initial_positions_file(self.golden_data_path / 'malformed_10_glowworms_no_anm_nan.dat')

    def test_read_initial_positions_file_malformed_missing(self):
        with pytest.raises(GSOCoordinatesError, match=r"Malformed line 2"):
            poses = read_initial_positions_file(self.golden_data_path / 'malformed_10_glowworms_no_anm_less.dat')

    def test_read_initial_positions_file_malformed_string(self):
        with pytest.raises(GSOCoordinatesError, match=r"Malformed pose in line 10"):
            poses = read_initial_positions_file(self.golden_data_path / 'malformed_10_glowworms_no_anm_string.dat')

    def test_read_initial_positions_file_not_found(self):
        with pytest.raises(GSOCoordinatesError, match=r"Cannot find or open"):
            poses = read_initial_positions_file(self.golden_data_path / 'wrong_file_name.dat')

    def test_distance_trans(self):
        poses = read_initial_positions_file(self.golden_data_path / 'good_10_glowworms_no_anm.dat')
        first = poses[0]
        last = poses[-1]

        assert 14.770195497 == pytest.approx(first.distance_trans(last))

    def test_distance_rot(self):
        poses = read_initial_positions_file(self.golden_data_path / 'good_10_glowworms_no_anm.dat')
        first = poses[0]
        last = poses[-1]

        assert 0.4755491286 == pytest.approx(first.distance_rot(last))

    def test_pose_repr(self):
        poses = read_initial_positions_file(self.golden_data_path / 'good_10_glowworms_no_anm.dat')
        first = poses[0]

        assert '(15.345706317, -14.691716422, 12.472000065, -0.105840692, -0.484936911, 0.599743027, -0.627648183)' == first.pose_repr()

    def test_pose_repr_anm(self):
        translation = [15.345706317, -14.691716422, 12.472000065]
        rotation = [-0.105840692, -0.484936911, 0.599743027, -0.627648183]
        rec_extents = [16.0, 17.0, 18.0]
        lig_extents = [19.0, 20.0, 21.0]
        pose = LightDockPose(translation, rotation, rec_extents, lig_extents)

        assert ('(15.345706317, -14.691716422, 12.472000065, -0.105840692, -0.484936911, 0.599743027, -0.627648183,'
                ' 16.000000000, 17.000000000, 18.000000000, 19.000000000, 20.000000000, 21.000000000)') == pose.pose_repr()

    def test_read_initial_positions_file_good_with_anm(self):

        poses = read_initial_positions_file(self.golden_data_path / 'good_10_glowworms_anm.dat', num_anm_rec=10, num_anm_lig=10)

        assert len(poses) == 10

        first = poses[0]
        assert first.rec_extents.size == 10
        assert first.lig_extents.size == 10
        assert np.allclose([4.670706612, 29.108347070, -3.667749226], first.translation)
        assert -0.105840692 == pytest.approx(first.rotation.w)
        assert -0.484936911 == pytest.approx(first.rotation.x)
        assert 0.599743027 == pytest.approx(first.rotation.y)
        assert -0.627648183 == pytest.approx(first.rotation.z)
        assert -1.057746227 == pytest.approx(first.rec_extents[0])
        assert -1.360405829 == pytest.approx(first.rec_extents[-1])
        assert -0.341520034 == pytest.approx(first.lig_extents[0])
        assert 0.002150371 == pytest.approx(first.lig_extents[-1])

        last = poses[-1]
        assert last.rec_extents.size == 10
        assert last.lig_extents.size == 10
        assert np.allclose([6.874675516, 29.180499401, 10.936907228], last.translation)
        assert -0.235809453 == pytest.approx(last.rotation.w)
        assert 0.215334879 == pytest.approx(last.rotation.x)
        assert -0.938641307 == pytest.approx(last.rotation.y)
        assert 0.130296927 == pytest.approx(last.rotation.z)
        assert 0.805922681 == pytest.approx(last.rec_extents[0])
        assert -0.487216972 == pytest.approx(last.rec_extents[-1])
        assert -0.934010084 == pytest.approx(last.lig_extents[0])
        assert 0.129468724 == pytest.approx(last.lig_extents[-1])


class TestPredictions:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / 'golden_data' / 'predictions_files'

    def test_read_predictions_file(self):

        predictions = read_predictions_file(self.golden_data_path / 'gso_100_good.out', 10, 10)

        assert len(predictions) == 10

        first = predictions[0]
        assert 5.7140324 == pytest.approx(first.pose.translation[0])
        assert 0.3004736 == pytest.approx(first.pose.rotation.w)
        assert 10 == first.pose.rec_extents.size
        assert 10 == first.pose.lig_extents.size
        assert -0.3569222 == pytest.approx(first.pose.lig_extents[-1])
        assert 0 == first.receptor_id
        assert 0 == first.ligand_id
        assert 10.64906252 == pytest.approx(first.luciferin)
        assert 11 == first.num_neighbors
        assert 0.08 == pytest.approx(first.vision_range)
        assert 7.02732804 == pytest.approx(first.scoring)

        last = predictions[-1]
        assert 3.4458835 == pytest.approx(last.pose.translation[0])
        assert -0.1630263 == pytest.approx(last.pose.rotation.w)
        assert 10 == last.pose.rec_extents.size
        assert 10 == last.pose.lig_extents.size
        assert 0.2243666 == pytest.approx(last.pose.lig_extents[-1])
        assert 0 == last.receptor_id
        assert 0 == last.ligand_id
        assert 10.29378007 == pytest.approx(last.luciferin)
        assert 3 == last.num_neighbors
        assert 1.96 == pytest.approx(last.vision_range)
        assert 6.84322529 == pytest.approx(last.scoring)

    def test_read_predictions_file_wrong_anm(self):
        with pytest.raises(GSOCoordinatesError, match=r"Malformed line"):
            predictions = read_predictions_file(self.golden_data_path / 'gso_100_good.out')

    def test_read_predictions_file_wrong_luciferin(self):
        with pytest.raises(GSOCoordinatesError, match=r"Malformed prediction"):
            predictions = read_predictions_file(self.golden_data_path / 'gso_100_wrong_luciferin.out', 10, 10)

    def test_read_predictions_file_not_found(self):
        with pytest.raises(GSOCoordinatesError, match=r"Cannot find or open"):
            poses = read_predictions_file(self.golden_data_path / 'wrong_file_name.dat')
