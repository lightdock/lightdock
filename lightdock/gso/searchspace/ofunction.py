"""The objective function to optimize using the GSO algorithm"""


class ObjectiveFunction(object):
    """Objective functions interface"""
    def __call__(self, coordinates):
        raise NotImplementedError()
