"""
A simple demonstration of using ADIOS2 directly, based on the high-level API examples at
https://adios2.readthedocs.io/en/latest/api_high/api_high.html#python-read-step-by-step-example
"""

from mpi4py import MPI

import numpy as np
import adios2

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

with adios2.open("data", "r",
                 config_file="adios2_config.xml",
                 io_in_config_file="solution-io") as fh:

    for fstep in fh:

        step = fstep.current_step()
        print(step)

        vars = fstep.available_variables()
        for name, info in vars.items():
            print(name)
            print(info)

        pp = fstep.read("pp")
        print(pp.shape)
