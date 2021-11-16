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

from .postprocess import *
from .tools import *
from .deriv import *

# Try to import h5py for use with HDF5 files
HAVE_HDF5=False
try:
    import h5py
except ImportError:
    pass
else:
    HAVE_HDF5 = True
    
#################################################################################

def report_p4i3d_status():
    """Reports a summary of Py4Incompact3D status after initialisation."""

    line_sep = "=" * 72

    print(line_sep)
    print("Py4Incompact3D init status:")
    print(f"- HDF5 enabled: {HAVE_HDF5}")
    print(line_sep)

report_p4i3d_status()
