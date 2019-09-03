"""Reads GSO parameters from configuration file"""

from configparser import ConfigParser
import os
from lightdock.error.lightdock_errors import GSOParameteresError


class GSOParameters(object):
    """Represents the set of the variables of the algorithm"""
    def __init__(self, file_name=None):
        self._config = ConfigParser()
        try:
            if file_name:
                self._config.readfp(open(file_name))
            else:
                self._config.readfp(open("%s%s" % (os.environ['LIGHTDOCK_CONF_PATH'], 'glowworm.conf')))
        except Exception as e:
            raise GSOParameteresError(str(e))
        
        try:
            self.rho = float(self._config.get('GSO', 'rho'))
            self.gamma = float(self._config.get('GSO', 'gamma'))
            self.beta = float(self._config.get('GSO', 'beta'))
            self.initial_luciferin = float(self._config.get('GSO', 'initialLuciferin'))
            self.initial_vision_range = float(self._config.get('GSO', 'initialVisionRange'))
            self.max_vision_range = float(self._config.get('GSO', 'maximumVisionRange'))
            self.max_neighbors = int(self._config.get('GSO', 'maximumNeighbors'))
        
        except Exception as e:
            raise GSOParameteresError("Problem parsing GSO parameters file. Details: %s" % str(e))
