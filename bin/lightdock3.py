#!/usr/bin/env python3

"""Execution controller

Depending on the environment, executes a MPI or multiprocessing version.
"""

import traceback
from lightdock.util.logger import LoggingManager
from lightdock.util.parser import CommandLineParser


log = LoggingManager.get_logger("lightdock")


if __name__ == "__main__":

    try:
        parser = CommandLineParser()
        mpi_support = parser.args.mpi
        if mpi_support:
            from lightdock.simulation.docking_mpi import (
                run_simulation as mpi_simulation,
            )

            mpi_simulation(parser)
        else:
            from lightdock.simulation.docking_multiprocessing import (
                run_simulation as multiprocessing_simulation,
            )

            multiprocessing_simulation(parser)

    except Exception:
        log.error("Lightdock has failed, please check traceback:")
        traceback.print_exc()
