#!/usr/bin/env python3

"""Calculates few statistics about the result of the simulation"""

import sys
import os

from lightdock.util.logger import LoggingManager
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

log = LoggingManager.get_logger('stats')


def usage():
    """Displays usage parameters and exists with error"""
    log.error("Wrong arguments")
    raise SystemExit("usage: %s number_of_steps number_of_glowworms" % (sys.argv[0]))


def parse_command_line():
    # Arguments parsing
    if len(sys.argv) != 3:
        usage()
    try:
        num_steps = int(sys.argv[1])
    except Exception as e:
        log.error(str(e))
        
    try:
        num_glowworms = int(sys.argv[2])
    except Exception as e:
        log.error(str(e))
    
    return num_steps, num_glowworms


def parse_file(file_name):
    """Parses a given GSO step output file. Return a list of luciferin, neighbors
    and vision range for each glowworm found in the file
    """
    num_glowworms = 0
    glowworms = []
    lines = open(file_name).readlines()[1:]
    for line in lines:
        if line[0] == '(':
            values = line.split(')')[1].split()
            luciferin = float(values[-4])
            neighbors = int(values[-3])
            vision_range = float(values[-2])
            glowworms.append({'luciferin':luciferin, 
                              'neighbors':neighbors,
                              'vision_range':vision_range}) 
            num_glowworms += 1
            
    return glowworms


def plot_stats(num_glowworms, num_steps, values, file_name, x_label='', y_label=''):
    """Plots for each glowworm and step, the given values"""
    fig = Figure(figsize=(20,20))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    ax.grid(True,linestyle='-',color='0.75')
    ax.set_xlabel(x_label,fontsize=12)
    ax.set_ylabel(y_label,fontsize=12)
    
    for i_glowworm in range(num_glowworms):
        ax.plot(range(1,num_steps+1), values[i_glowworm]);
    canvas.print_figure(file_name,gdpi=500)
    print("Generated %s plot." % file_name)


if __name__ == "__main__":
    # Parse arguments
    num_steps, num_glowworms = parse_command_line()
    
    # Output csv files
    energies_file_name = 'energies.csv'
    neighbors_file_name = 'neighbors.csv'
    vision_range_file_name = 'vision_range.csv'
    
    energies_file = open(energies_file_name, 'w')
    neighbors_file = open(neighbors_file_name, 'w')
    vision_range_file = open(vision_range_file_name, 'w')
    
    # Matrices to store values per glowworm/row
    energies = [[] for i in range(num_glowworms)]
    neighbors = [[] for i in range(num_glowworms)]
    vision_range = [[] for i in range(num_glowworms)]
    
    for step in range(1, num_steps+1):
        # Parse each stored step file
        file_name = 'gso_%d.out' % step 
        glowworm_stats = parse_file(file_name)
        
        # Write step at the beginning of the files
        energies_file.write("%4d  " % step)
        neighbors_file.write("%4d  " % step)
        vision_range_file.write("%4d  " % step)
        
        for i_glowworm, stat in enumerate(glowworm_stats):
            # Energy
            energies_file.write("%8.3f  " % stat['luciferin'])
            energies[i_glowworm].append(stat['luciferin'])
            # Number of neighbors
            neighbors_file.write("%d  " % stat['neighbors'])
            neighbors[i_glowworm].append(stat['neighbors'])
            # Vision range
            vision_range_file.write("%5.2f  " % stat['vision_range'])
            vision_range[i_glowworm].append(stat['vision_range'])
        
        energies_file.write(os.linesep)
        neighbors_file.write(os.linesep)
        vision_range_file.write(os.linesep)
        
    energies_file.close()
    neighbors_file.close()
    vision_range_file.close()
    
    # Plot energies
    plot_stats(num_glowworms, num_steps, energies, 'energies.png', 'Step', 'Luciferin')
    
    # Plot neighbors
    plot_stats(num_glowworms, num_steps, neighbors, 'neighbors.png', 'Step', 'Number of neighbors')
    
    # Plot vision range
    plot_stats(num_glowworms, num_steps, vision_range, 'vision_range.png', 'Step', 'Vision range')
