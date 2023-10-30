"""Tests for RandomNumberGeneratorFromFile class"""

import pytest
from os import linesep
from pathlib import Path
from lightdock.mathutil.lrandom import RandomNumberGeneratorFromFile
from lightdock.error.lightdock_errors import RandomNumberError


class TestRandomNumberGeneratorFromFile:
    def setup_class(self):
        self.path = Path(__file__).absolute().parent
        self.golden_data_path = self.path / "golden_data"

    def test_get_random_number(self):
        gen = RandomNumberGeneratorFromFile(
            self.golden_data_path / "random_numbers.txt"
        )
        expected_file = open(self.golden_data_path / "random_numbers.txt")
        for line in expected_file:
            if line.startswith("#seed"):
                seed = int(line.rstrip(linesep).split("=")[1])
            else:
                expected = float(line)
                assert expected == pytest.approx(gen())
        assert seed == 25

    def test_get_list_exhausted(self):
        with pytest.raises(RandomNumberError):
            gen = RandomNumberGeneratorFromFile(
                self.golden_data_path / "random_numbers.txt"
            )
            for i in range(101):
                gen()
            assert i != 101

    def test_wrong_seed(self):
        with pytest.raises(RandomNumberError):
            gen = RandomNumberGeneratorFromFile(self.golden_data_path / "wrong_seed.txt")
            assert gen()

    def test_wrong_line(self):
        gen = RandomNumberGeneratorFromFile(self.golden_data_path / "wrong_line.txt")
        count = 0
        for _ in range(99):
            gen()
            count += 1
        assert count == 99
