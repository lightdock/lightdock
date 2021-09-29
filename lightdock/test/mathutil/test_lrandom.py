"""Tests for RandomNumberGenerator interface"""

from lightdock.mathutil.lrandom import (
    RandomNumberGenerator,
    NormalGenerator,
    NMExtentGenerator,
)
from nose.tools import raises
import numpy as np


class TestRandomNumberGenerator:
    @raises(NotImplementedError)
    def test_call_interface(self):
        gen = RandomNumberGenerator()
        gen()


class TestNormalGenerator:
    def test_using_generator(self):
        rnd = NormalGenerator(seed=666, mu=0.0, sigma=5.0)

        random_numbers = [rnd() for _ in range(10)]

        expected = [
            4.1209404165876515,
            2.3998300155052408,
            5.867340060460123,
            4.545240343925104,
            -2.858607259497295,
            -0.5474863386290562,
            0.095141324542312,
            -4.718805322646974,
            3.2028657656895856,
            -3.932215858160615,
        ]

        assert np.allclose(expected, random_numbers)


class TestNMExtentGenerator:
    def test_using_generator(self):
        rnd = NMExtentGenerator(seed=666, mu=0.0, sigma=3.0)

        random_numbers = [rnd() for _ in range(10)]

        expected = [
            2.472564249952591,
            1.4398980093031444,
            3.5204040362760733,
            2.7271442063550624,
            1.715164355698377,
            0.32849180317743365,
            0.057084794725387196,
            2.8312831935881846,
            1.9217194594137512,
            2.359329514896369,
        ]

        assert np.allclose(expected, random_numbers)
