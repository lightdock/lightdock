"""Tests for RandomNumberGeneratorFromFile class"""

import os
import shutil
from pathlib import Path
from nose.tools import assert_almost_equals, raises
from lightdock.mathutil.lrandom import RandomNumberGeneratorFromFile
from lightdock.error.lightdock_errors import RandomNumberError


class TestRandomNumberGeneratorFromFile:
    def __init__(self):
        self.path = Path(__file__).absolute().parent
        self.test_path = self.path / "scratch_random_file"
        self.golden_data_path = self.path / "golden_data"

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

    def test_get_random_number(self):
        gen = RandomNumberGeneratorFromFile(
            self.golden_data_path / "random_numbers.txt"
        )
        expected_file = open(self.golden_data_path / "random_numbers.txt")
        for line in expected_file:
            if line.startswith("#seed"):
                seed = int(line.rstrip(os.linesep).split("=")[1])
            else:
                expected = float(line)
                assert_almost_equals(expected, gen())
        assert seed == 25

    @raises(RandomNumberError)
    def test_get_list_exhausted(self):
        gen = RandomNumberGeneratorFromFile(
            self.golden_data_path / "random_numbers.txt"
        )
        for i in range(101):
            gen()
        assert i != 101

    @raises(RandomNumberError)
    def test_wrong_seed(self):
        gen = RandomNumberGeneratorFromFile(self.golden_data_path / "wrong_seed.txt")
        assert gen()

    def test_wrong_line(self):
        gen = RandomNumberGeneratorFromFile(self.golden_data_path / "wrong_line.txt")
        count = 0
        for _ in range(99):
            gen()
            count += 1
        assert count == 99
