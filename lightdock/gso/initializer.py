"""Module to generate initial populations of glowworms agents used by the GSO algorithm"""

from lightdock.gso.swarm import Swarm
from lightdock.gso.coordinates import CoordinatesFileReader, Coordinates
from lightdock.error.lightdock_errors import GSOCoordinatesError
from lightdock.gso.searchspace.landscape import LandscapePosition, DockingLandscapePosition
from lightdock.mathutil.lrandom import MTGenerator


class Initializer(object):
    """Generates a population of glowworms.
    
    The landscape is determined by the ObjectiveFunction object.
    The glowworms will use the given algorithm parameters.
    """
    def __init__(self, objective_functions, number_of_glowworms, gso_parameters):
        self.objective_functions = objective_functions
        self.number_of_glowworms = number_of_glowworms
        self.parameters = gso_parameters
        self.positions = []
    
    def generate_glowworms(self):
        """Creates an initial population of glowworms"""
        self.positions = self.generate_landscape_positions()
        return Swarm(self.positions, self.parameters)
    
    def generate_landscape_positions(self):
        """Generates the initial positions of each glowworm"""
        raise NotImplementedError()


class RandomInitializer(Initializer):
    """Generates a population of glowworms with random positions"""
    def __init__(self, objective_functions, number_of_glowworms, gso_parameters,
                 bounding_box, random_number_generator):
        super(RandomInitializer, self).__init__(objective_functions, number_of_glowworms, gso_parameters)
        self.bounding_box = bounding_box
        self.random_number_generator = random_number_generator
    
    def generate_landscape_positions(self):
        """Generates a list of landscape positions that have been read
        from initial_population_file.
        """
        positions = []
        for index in range(self.number_of_glowworms):
            coordinates = []
            for dimension in range(self.bounding_box.dimension):
                bound = self.bounding_box.get_boundary_of_dimension(dimension)
                coord = self.random_number_generator(bound.lower_limit, bound.upper_limit)
                coordinates.append(coord)
            positions.append(LandscapePosition(self.objective_functions[0], Coordinates(coordinates)))
        return [positions]


class FromFileInitializer(Initializer):
    """Generates a population of glowworms with initial positions that
    are read from a given file.
    """
    def __init__(self, objective_functions, number_of_glowworms, gso_parameters,
                 dimensions, initial_population_file):
        super(FromFileInitializer, self).__init__(objective_functions,
                                                  number_of_glowworms,
                                                  gso_parameters)
        self.dimensions = dimensions
        self.initial_population_file = initial_population_file
    
    def generate_landscape_positions(self):
        """Generates a list of landscape positions that have been read
        from initial_population_file.
        """
        reader = CoordinatesFileReader(self.dimensions)
        coordinates = reader.get_coordinates_from_file(self.initial_population_file)
        
        if not coordinates:
            raise GSOCoordinatesError("No coordinates have been read from %s file" % self.initial_population_file)
        
        if len(coordinates) != self.number_of_glowworms:
            raise GSOCoordinatesError("Number of coordinates read and number of glowworms does not correspond")
        
        positions = []
        for index in range(self.number_of_glowworms):
            positions.append(LandscapePosition(self.objective_functions[0], coordinates[index]))

        return [positions]

    
class LightdockFromFileInitializer(Initializer):
    """This initializer takes into account the complex provided by the adapter"""
    def __init__(self, adapters, scoring_functions, number_of_glowworms, gso_parameters,
                 dimensions, initial_population_file, step_translation, step_rotation,
                 random_number_generator, step_nmodes, anm_rec, anm_lig):
        super(LightdockFromFileInitializer, self).__init__(scoring_functions, number_of_glowworms, gso_parameters)
        self.dimensions = dimensions
        self.initial_population_file = initial_population_file
        self.adapters = adapters
        self.step_translation = step_translation
        self.step_rotation = step_rotation
        self.step_nmodes = step_nmodes
        # Patch to not mess with old simulations
        self.random_number_generator = MTGenerator(random_number_generator.seed)
        self.anm_rec = anm_rec
        self.anm_lig = anm_lig
    
    def generate_landscape_positions(self):
        """Generates a list of landscape positions that have been read
        from initial_population_file.
        """
        reader = CoordinatesFileReader(self.dimensions)
        coordinates = reader.get_coordinates_from_file(self.initial_population_file)
        
        if not coordinates:
            raise GSOCoordinatesError("No coordinates have been read from %s file" % self.initial_population_file)
        
        if len(coordinates) != self.number_of_glowworms:
            raise GSOCoordinatesError("Number of coordinates read and number of glowworms does not correspond")
        
        positions = []
        for i, adapter in enumerate(self.adapters):
            positions.append([])
            for index in range(self.number_of_glowworms):
                receptor_index = self.random_number_generator.randint(0, len(self.adapters[0].receptor_model)-1)
                ligand_index = self.random_number_generator.randint(0, len(self.adapters[0].ligand_model)-1)
                positions[i].append(DockingLandscapePosition(self.objective_functions[i], coordinates[index],
                                                             adapter.receptor_model, adapter.ligand_model,
                                                             receptor_index, ligand_index,
                                                             self.step_translation, self.step_rotation,
                                                             self.step_nmodes, self.anm_rec, self.anm_lig))
        return positions
