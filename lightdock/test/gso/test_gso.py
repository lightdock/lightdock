"""Tests for GSO class"""

import os
from pathlib import Path
from math import pi
import filecmp
import shutil
from lightdock.gso.parameters import GSOParameters
from lightdock.gso.searchspace.benchmark_ofunctions import J1, J2, J3, J4, J5
from lightdock.gso.algorithm import GSOBuilder
from lightdock.gso.boundaries import Boundary, BoundingBox
from lightdock.mathutil.lrandom import MTGenerator


class TestGSOBuilderInJ1:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_gso"
        self.golden_data_path = self.path / "golden_data"
        self.gso_parameters = GSOParameters()
        self.gso_parameters.initial_vision_range = 3.0
        self.gso_parameters.max_vision_range = 3.0

    def setup(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass
        os.mkdir(self.test_path)

    def teardown(self):
        try:
            shutil.rmtree(self.test_path)
        except OSError:
            pass

    def found_peaks(
        self, peak_coordinates, dimension, glowworms, minimum_matches=3, tolerance=0.05
    ):
        peaks_found = [False for _ in range(len(peak_coordinates))]
        for i_peak, peak in enumerate(peak_coordinates):
            for glowworm in glowworms:
                coordinates = glowworm.landscape_positions[0].coordinates
                found = True
                for j in range(dimension):
                    found = found and (abs(coordinates[j] - peak[j]) <= tolerance)
                if found:
                    peaks_found[i_peak] = True
        return minimum_matches <= peaks_found.count(True)

    def test_GSO_with_J1(self):
        objective_function = J1()
        bounding_box = BoundingBox([Boundary(-3.0, 3.0), Boundary(-3.0, 3.0)])
        number_of_glowworms = 50
        random_number_generator = MTGenerator(324324)
        builder = GSOBuilder()
        gso = builder.create(
            number_of_glowworms,
            random_number_generator,
            self.gso_parameters,
            objective_function,
            bounding_box,
        )

        gso.run(200)

        # Function peak coordinates
        peak_coordinates = [[1.28, 0.0], [0.0, 1.58], [-0.46, -0.63]]

        # Check with auxiliar function the position of the glowworms
        assert self.found_peaks(peak_coordinates, 2, gso.swarm.glowworms)

    def test_GSO_with_J1_initializing_from_file(self):
        objective_function = J1()
        bounding_box = BoundingBox([Boundary(-3.0, 3.0), Boundary(-3.0, 3.0)])
        number_of_glowworms = 50
        random_number_generator = MTGenerator(324324)
        builder = GSOBuilder()
        gso = builder.create_from_file(
            number_of_glowworms,
            random_number_generator,
            self.gso_parameters,
            objective_function,
            bounding_box,
            self.golden_data_path / "initial_positions.txt",
        )

        gso.run(200)

        # Function peak coordinates
        peak_coordinates = [[1.28, 0.0], [0.0, 1.58], [-0.46, -0.63]]

        # Check with auxiliar function the position of the glowworms
        assert self.found_peaks(peak_coordinates, 2, gso.swarm.glowworms, 2)

    def test_GSO_with_J2(self):
        objective_function = J2()
        self.gso_parameters.initial_vision_range = 2.0
        self.gso_parameters.max_vision_range = 2.0
        bounding_box = BoundingBox([Boundary(-1.0, 1.0), Boundary(-1.0, 1.0)])
        number_of_glowworms = 70
        random_number_generator = MTGenerator(324324)
        builder = GSOBuilder()
        gso = builder.create(
            number_of_glowworms,
            random_number_generator,
            self.gso_parameters,
            objective_function,
            bounding_box,
        )

        gso.run(200)

        # Function peak coordinates
        peak_coordinates = [[-0.5, -0.5], [-0.5, 0.5], [0.5, -0.5], [0.5, 0.5]]

        # Check with auxiliar function the position of the glowworms
        assert self.found_peaks(peak_coordinates, 2, gso.swarm.glowworms, 3)

    def test_GSO_with_J3(self):
        objective_function = J3()
        self.gso_parameters.initial_vision_range = 2.0
        self.gso_parameters.max_vision_range = 2.0
        bounding_box = BoundingBox([Boundary(-10.0, 10.0), Boundary(-10.0, 10.0)])
        number_of_glowworms = 70
        random_number_generator = MTGenerator(324324)
        builder = GSOBuilder()
        gso = builder.create(
            number_of_glowworms,
            random_number_generator,
            self.gso_parameters,
            objective_function,
            bounding_box,
        )

        gso.run(50)

        # Save last step
        gso.swarm.save(50, self.test_path, "gso_j3_50.out")

        assert filecmp.cmp(
            self.test_path / "gso_j3_50.out", self.golden_data_path / "gso_j3_50.out"
        )

    def test_GSO_with_J4(self):
        objective_function = J4()
        self.gso_parameters.initial_vision_range = 0.75
        self.gso_parameters.max_vision_range = 0.75
        bounding_box = BoundingBox([Boundary(-2.0, 2.0), Boundary(-2.0, 2.0)])
        number_of_glowworms = 100
        random_number_generator = MTGenerator(324324)
        builder = GSOBuilder()
        gso = builder.create(
            number_of_glowworms,
            random_number_generator,
            self.gso_parameters,
            objective_function,
            bounding_box,
        )

        gso.run(50)

        # Save last step
        gso.swarm.save(50, self.test_path, "gso_j4_50.out")

        assert filecmp.cmp(
            self.test_path / "gso_j4_50.out", self.golden_data_path / "gso_j4_50.out"
        )

    def test_GSO_with_J5(self):
        objective_function = J5()
        self.gso_parameters.initial_vision_range = 3.0
        self.gso_parameters.max_vision_range = 3.0
        bounding_box = BoundingBox(
            [Boundary(-2.0 * pi, 2.0 * pi), Boundary(-2.0 * pi, 2.0 * pi)]
        )
        number_of_glowworms = 100
        random_number_generator = MTGenerator(324324)
        builder = GSOBuilder()
        gso = builder.create(
            number_of_glowworms,
            random_number_generator,
            self.gso_parameters,
            objective_function,
            bounding_box,
        )

        gso.run(70)

        # Save last step
        gso.swarm.save(70, self.test_path, "gso_j5_70.out")

        assert filecmp.cmp(
            self.test_path / "gso_j5_70.out", self.golden_data_path / "gso_j5_70.out"
        )

    def test_GSO_with_report_in_file_and_saving_intermediary_files(self):
        objective_function = J5()
        self.gso_parameters.initial_vision_range = 3.0
        self.gso_parameters.max_vision_range = 3.0
        bounding_box = BoundingBox(
            [Boundary(-2.0 * pi, 2.0 * pi), Boundary(-2.0 * pi, 2.0 * pi)]
        )
        number_of_glowworms = 5
        random_number_generator = MTGenerator(324324)
        builder = GSOBuilder()
        gso = builder.create_from_file(
            number_of_glowworms,
            random_number_generator,
            self.gso_parameters,
            objective_function,
            bounding_box,
            self.golden_data_path / "initial_positions_redux.txt",
        )

        gso.run(
            5,
            saving_path=self.test_path,
            save_intermediary=True,
            save_all_intermediary=True,
        )

        for i in range(5):
            assert (self.test_path / f"gso_{i+1}.out").exists()

        gso.report(self.test_path / "report.out")
        lines = open(self.test_path / "report.out").readlines()
        assert len(lines) == 14

    def test_GSO_with_report(self):
        objective_function = J5()
        self.gso_parameters.initial_vision_range = 3.0
        self.gso_parameters.max_vision_range = 3.0
        bounding_box = BoundingBox(
            [Boundary(-2.0 * pi, 2.0 * pi), Boundary(-2.0 * pi, 2.0 * pi)]
        )
        number_of_glowworms = 5
        random_number_generator = MTGenerator(324324)
        builder = GSOBuilder()
        gso = builder.create(
            number_of_glowworms,
            random_number_generator,
            self.gso_parameters,
            objective_function,
            bounding_box,
        )

        gso.run(5)

        report = gso.report()
        assert len(report.split(os.linesep)) == 14
        assert str(gso) == report
