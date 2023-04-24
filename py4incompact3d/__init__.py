# Copyright 2018 G. Deskos
# Copyright 2021 University of Edinburgh
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.
"""
FILE: __init__.py
"""

from mpi4py import MPI

from .extras import *

from .postprocess import *
from .tools import *
from .deriv import *
    
# Set MPI variables
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
    
#################################################################################

def report_p4i3d_status():
    """Reports a summary of Py4Incompact3D status after initialisation."""

    line_sep = "=" * 72

    print(line_sep)
    print("Py4Incompact3D init status:")
    print("- MPI:")
    print(f"+- running on {size} ranks")
    print(f"- HDF5 enabled: {extras.HAVE_HDF5}")
    print(f"- ADIOS2 enabled: {extras.HAVE_ADIOS2}")
    print(line_sep)

if (rank == 0):
    report_p4i3d_status()

def toggle_hdf5():
    """Toggle HDF5 support on/off."""

    if HDF5_AVAILABLE:
        extras.HAVE_HDF5 = not extras.HAVE_HDF5
    elif (rank == 0):
        print("WARNING: HDF5 is not available!")

    if (rank == 0):
        report_p4i3d_status()

def toggle_adios2():
    """Toggle ADIOS2 support on/off."""

    if ADIOS2_AVAILABLE:
        HAVE_ADIOS2 = not HAVE_ADIOS2
    elif (rank == 0):
        print("WARNING: ADIOS2 is not available!")

    if (rank == 0):
        report_p4i3d_status()
