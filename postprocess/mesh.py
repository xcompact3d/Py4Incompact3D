# Copyright 2018 Georgios Deskos

# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.

import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import griddata

class Mesh():
    """
    Mesh is a model object representing

    """

    def __init__(self,instance_dictionary):

        super().__init__()

        self.description = instance_dictionary["description"]

        properties = instance_dictionary["properties"]
        self.Nx  = properties["Nx"]
        self.Ny  = properties["Ny"]
        self.Nz  = properties["Nz"]
        self.Lx  = properties["Lx"]
        self.Ly  = properties["Ly"]
        self.Lz  = properties["Lz"]
        self.BCx = properties["BCx"]
        self.BCy = properties["BCy"]
        self.BCz = properties["BCz"]

        # Once we know the mesh layout we can set the derivative variables
        self.compute_derivvars()

    def compute_derivvars(self):
        """ Compute variables required by derivative functions. """
	if (self.BCx==0):
        	self.dx = self.Lx / float(self.Nx)
        else:
        	self.dx = self.Lx / float(self.Nx-1)
	if (self.BCy==0):	
		self.dy = self.Ly / float(self.Ny) # XXX This will not be correct for stretched grids
        else:
		self.dy = self.Ly / float(self.Ny-1) # XXX This will not be correct for stretched grids
	if (self.BCz==0):	
        	self.dz = self.Lz / float(self.Nz)
        else:
		self.dz=self.Lz/float(self.Nz-1)

        self.alpha = 1.0 / 3.0
        self.a = 14.0 / 9.0
        self.b = 1.0 / 9.0
