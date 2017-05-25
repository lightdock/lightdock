"""Tests for RandomNumberGeneratorFromFile class"""

from nose.tools import assert_almost_equals
from nose.tools import raises
import os
from lightdock.mathutil.lrandom import RandomNumberGeneratorFromFile
from lightdock.error.lightdock_errors import RandomNumberError


class TestRandomNumberGeneratorFromFile:

    def setUp(self):
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'

    def tearDown(self):
        pass
    
    def test_get_random_number(self):
        gen = RandomNumberGeneratorFromFile(self.golden_data_path + 'random_numbers.txt')
        expected_file = open(self.golden_data_path + 'random_numbers.txt')
        for line in expected_file:
            if line.startswith('#seed'):
                seed = int(line.rstrip(os.linesep).split('=')[1])
            else:
                expected = float(line)
                assert_almost_equals(expected, gen())
        assert 25 == seed
    
    @raises(RandomNumberError)
    def test_get_list_exhausted(self):
        gen = RandomNumberGeneratorFromFile(self.golden_data_path + 'random_numbers.txt')
        for i in range(101):
            gen()
        assert i != 101
        
    @raises(RandomNumberError)
    def test_wrong_seed(self):
        gen = RandomNumberGeneratorFromFile(self.golden_data_path + 'wrong_seed.txt')
        assert gen()
        
    def test_wrong_line(self):
        gen = RandomNumberGeneratorFromFile(self.golden_data_path + 'wrong_line.txt')
        count = 0
        for i in range(99):
            gen()
            count += 1
        assert 99 == count
