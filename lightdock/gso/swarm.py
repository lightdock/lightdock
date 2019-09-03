"""The set of swarms of glowworm agents used in the algorithm"""

from operator import attrgetter
from lightdock.gso.glowworm import Glowworm


class Swarm(object):
    """A swarm of glowworms"""
    def __init__(self, landscape_positions, parameters):
        """Creates a glowworm population using a landscape_positons list and parameters"""
        positions_per_glowworm = [[] for _ in range(len(landscape_positions[0]))]
        for function in landscape_positions:
            for glowworm_id, position in enumerate(function):
                positions_per_glowworm[glowworm_id].append(position)
        self.glowworms = [Glowworm(positions, parameters) for positions in positions_per_glowworm]
        self.docking = (landscape_positions[0][0].__class__.__name__ != 'LandscapePosition')

    def update_luciferin(self):
        """Updates luciferin of each glowworm"""
        for glowworm in self.glowworms:
            glowworm.compute_luciferin()

    def movement_phase(self, rnd_generator):
        """Updates luciferin and probabilities of each glowworm to move if required
        following GSO algorithm.
        """
        selected = []
        positions = {}
        num_glowworms = self.get_size()
        for i in range(num_glowworms):
            glowworm = self.glowworms[i]
            glowworm.search_neighbors(self.glowworms)
            glowworm.compute_probability_moving_toward_neighbor()
            selected.append(glowworm.select_random_neighbor(rnd_generator()))
            positions[i] = [landscape_position.clone() for landscape_position in selected[-1].landscape_positions]

        for i in range(num_glowworms):
            glowworm = self.glowworms[i]
            neighbor = selected[i]
            position = positions[i]
            glowworm.move(neighbor, position)
            glowworm.update_conformers(neighbor, rnd_generator)
            glowworm.update_vision_range()

    def minimize_best(self):
        """Minimizes the glowworm with better energy using a local non-gradient minimization method"""
        best_glowworm = max(self.glowworms, key=attrgetter('scoring'))
        best_glowworm.minimize()

    def get_size(self):
        """Gets the population size of this swarm of glowworms"""
        return len(self.glowworms)

    def save(self, step, destination_path, file_name=''):
        """Saves actual population status to a file"""
        if file_name:
            dest_file_name = "%s/%s" % (destination_path, file_name)
        else:
            dest_file_name = "%s/gso_%d.out" % (destination_path, step)

        dest_file = open(dest_file_name, 'w')
        dest_file.write(str(self))
        dest_file.close()

    def __repr__(self):
        """String representation of the population"""
        if self.docking:
            representation = "#Coordinates  RecID  LigID  Luciferin  Neighbor's number  Vision Range  Scoring\n"
        else:
            representation = "#Coordinates  Luciferin  Neighbor's number  Vision Range  Scoring\n"
        for glowworm in self.glowworms:
            representation += str(glowworm) + "\n"
        return representation
