"""Execution controller

Depending on the environment, executes a MPI or multiprocessing version.
"""

# Check if Python interpreter version is suitable
import platform
import traceback
py_version = float(platform.python_version()[:3])
if py_version < 2.7:
    raise SystemExit("lightdock [ERROR] required Python version is 2.7.x")

from lightdock.util.logger import LoggingManager
from lightdock.util.parser import CommandLineParser


log = LoggingManager.get_logger('lightdock')


if __name__ == "__main__":

    try:
        parser = CommandLineParser()
        mpi_support = parser.args.mpi
        if mpi_support:
            from docking_mpi import run_simulation as mpi_simulation
            mpi_simulation(parser)
        else:
            from docking_multiprocessing import run_simulation as multiprocessing_simulation
            multiprocessing_simulation(parser)

    except Exception, e:
        log.error("Lightdock has failed, please check traceback:")
        traceback.print_exc()
