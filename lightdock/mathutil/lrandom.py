"""Random number generator class for wrapping Python implementation"""

import random
import os
import numpy as np
from lightdock.error.lightdock_errors import RandomNumberError


class RandomNumberGenerator(object):
    """Random number generator interface"""

    def __call__(self):
        raise NotImplementedError()


class MTGenerator(RandomNumberGenerator):
    """Python uses the Mersenne Twister as the core generator.
    It produces 53-bit precision floats and has a period of 2**19937-1
    """

    def __init__(self, seed):
        self.seed = seed
        self.random = random.Random()
        self.random.seed(self.seed, version=1)

    def __call__(self, lower_limit=0.0, upper_limit=1.0):
        return self.random.uniform(lower_limit, upper_limit)

    def randint(self, lower_limit=0, upper_limit=9):
        return int(self() * (upper_limit + 1)) + lower_limit


class RandomNumberGeneratorFromFile(RandomNumberGenerator):
    """Class to interact with a previously generated list of random numbers
    uniformly distributed between 0.0 and 1.0
    """

    def __init__(self, file_name):
        """Reads file_name which contains for each line a random number between
        0.0 and 1.0.
        """
        self._numbers = []
        self._index = 0
        numbers_file = open(file_name)
        for line in numbers_file:
            if line.startswith("#seed"):
                try:
                    self.seed = int(line.rstrip(os.linesep).split("=")[1])
                except:
                    raise RandomNumberError("Invalid seed")
            else:
                try:
                    self._numbers.append(float(line))
                except:
                    pass

    def __call__(self):
        try:
            number = self._numbers[self._index]
            self._index += 1
            return number
        except:
            raise RandomNumberError("Not enough random numbers")


class NormalGenerator(RandomNumberGenerator):
    """Generates random numbers following a gaussian distribution"""

    def __init__(self, seed, mu, sigma):
        self.seed = seed
        self.mu = mu
        self.sigma = sigma
        np.random.seed(self.seed)
        self.random = np.random

    def __call__(self):
        return self.random.normal(self.mu, self.sigma)


class NMExtentGenerator(RandomNumberGenerator):
    """Generates random numbers following a gaussian-uniform distribution"""

    def __init__(self, seed, mu, sigma):
        self.seed = seed
        self.mu = mu
        self.sigma = sigma
        np.random.seed(self.seed)
        self.random = np.random

    def __call__(self):
        n = abs(self.random.normal(self.mu, self.sigma))
        if n < self.mu:
            return self.random.uniform(0.0, self.mu)
        return n
