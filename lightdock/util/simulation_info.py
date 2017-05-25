import os
import time
from lightdock.constants import DEFAULT_LIGHTDOCK_INFO
from lightdock.error.lightdock_errors import LightDockError


def show_parameters(log, parser):
    """Informs the user about the parameters of the simulation"""
    log.info("Parameters are:")
    log.info("Receptor: %s" % parser.args.info_receptor_pdb)
    log.info("Ligand: %s" % parser.args.info_ligand_pdb)
    log.info("Number of clusters: %s" % parser.args.clusters)
    log.info("Number of glowworms per cluster: %s" % parser.args.glowworms)
    log.info("Simulation steps: %s" % parser.args.steps)
    log.info("GSO seed: %s" % parser.args.gso_seed)
    log.info("Starting points seed: %s" % parser.args.starting_points_seed)
    log.info("Translation step: %s" % parser.args.translation_step)
    log.info("Rotation step: %s" % parser.args.rotation_step)
    if parser.args.configuration_file:
        log.info("Configuration file: %s" % parser.args.info_configuration_file)
    if parser.args.ftdock_file:
        log.info("FTDock file: %s" % parser.args.info_ftdock_file)


def create_simulation_info_file(parser, path='.', file_name=DEFAULT_LIGHTDOCK_INFO):
    """Creates a simulation file from which recover from in a new simulation"""
    # Create the simulation info file. If it exists, includes a number
    # in the extension to avoid collision
    output_file_name = os.path.join(path, file_name)
    if os.path.isfile(output_file_name):
        original_file_name = output_file_name
        i = 1
        while os.path.isfile(output_file_name) and i < 255:
            output_file_name = "%s.%d" % (original_file_name, i)
            i += 1
        if i == 255:
            raise LightDockError('Too many simulation files')

    # Data to store
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    data = {'start_time': now}
    data.update(vars(parser.args))

    # Store the data in the file sorted alphabetically
    with open(output_file_name, 'w') as output:
        for key in sorted(data.keys()):
            output.write("%s=%s%s" % (str(key), str(data[key]), os.linesep))

    return output_file_name
