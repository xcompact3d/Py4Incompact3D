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

    def get_grid(self):

        x = np.zeros(self.Nx)
        y = np.zeros(self.Ny)
        z = np.zeros(self.Nz)

        for i in range(self.Nx):
            x[i] = i * self.dx
        for i in range(self.Ny):
            y[i] = i * self.dy
        for i in range(self.Nz):
            z[i] = i * self.dz

        return x, y, z

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

    def get_pencil_layout(self, comm_size, comm_rank):
        """ Determine optimal pencil layout for the given mesh and MPI communicator.

        :param comm_size: Size of the MPI communicator.
        :param comm_rank: Rank of the processor within the MPI communicator.

        :type comm_size: int
        :type comm_rank: int
        """

        # XXX Sub-optimal, good enough maybe?
        self.prow = int(self.Ny / self.Nx) * comm_size
        self.pcol = int(comm_size / self.prow)
        assert(self.prow * self.pcol == comm_size)

    def get_pencil(self, axis_size, ndiv, pgrid_distance):
        """ Determine the size of the pencil to be assigned to a rank.

        :param axis_size: Size of the axis.
        :param ndiv: Number of ranks to divide axis by.
        :param pgrid_distance: Distance along the axis in the processor grid.

        :type axis_size: int
        :type ndiv: int
        :type pgrid_distance: int

        By convention the axis is first split equally using integer division, the remainder is then
        distributed to the ranks lower in the processor grid by increasing their size by one.

        The grids are laid out as:

        Xgrid
       y^
        | 3 4 5
        | 0 1 2
        +------>z

        Ygrid
       z^
        | 3 4 5
        | 0 1 2
        +------>x

        Zgrid
        ^
        | 3 4 5
        | 0 1 2
        +------>
        """

        size = axis_size // ndiv

        start = pgrid_distance * size

        rem = axis_size - ndiv * size
        if (pgrid_distance < rem):
            size += 1

        start += min(rem, pgrid_distance)
        end = start + (size - 1)

        return start, end, size
        
    def decompose2d(self, comm_size, comm_rank):
        """ Decompose a mesh using 2D pencil decomposition. 

        :param comm_size: Size of the MPI communicator.
        :param comm_rank: Rank of the processor within the MPI communicator.

        :type comm_size: int
        :type comm_rank: int
        """

        # Setup the decomposition
        self.get_pencil_layout(comm_size, comm_rank)

        # Determine pencil sizes, starts, stops...
        self.xsize = [0] * 3
        self.ysize = [0] * 3
        self.zsize = [0] * 3

        self.xsize[0] = self.Nx
        self.xsize[1] = self.get_pencil_size(self.Ny, p_row)
        self.xsize[2] = self.get_pencil_size(self.Nz, p_col)
