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

        self.pen_start = [[0, 0, 0],
                          [0, 0, 0],
                          [0, 0, 0]]
        self.pen_end = [[self.Nx, 0, 0],
                        [0, self.Ny, 0],
                        [0, 0, self.Nz]]
        
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

    def compute_pencil(self, axis, rank, prow, pcol):
        """ """

        if 0 == axis:
            self.start[0][1] = (self.Ny // prow) * (rank // pcol)
            self.end[0][1] = min((self.Ny // prow) * (1 + rank // pcol), self.Ny)

            self.start[0][2] = (self.Nz // pcol) * (rank - (rank // prow) * pcol)
            self.end[0][2] = min((self.Nz // pcol) * (1 + (rank - (rank // prow) * pcol)), self.Nz)
        elif 1 == axis:
            self.start[1][0] = (self.Nx // pcol) * (rank - (rank // prow) * pcol)
            self.end[1][0] = min((self.Nx // pcol) * (1 + (rank - (rank // prow) * pcol)), self.Nx)

            self.start[1][2] = (self.Nz // prow) * (rank // pcol)
            self.end[1][2] = min((self.Nz // prow) * (1 + rank // pcol), self.Nz)
        else:
            self.start[2][0] = (self.Nx // pcol) * (rank - (rank // prow) * pcol)
            self.end[2][0] = min((self.Nx // pcol) * (1 + (rank - (rank // prow) * pcol)), self.Nx)

            self.start[2][1] = (self.Ny // prow) * (rank // pcol)
            self.end[2][1] = min((self.Ny // prow) * (1 + rank // pcol), self.Ny)
        
    def compute_decomposition(self, comm_rank, comm_size):
        """ Compute the "ideal" decomposition. """

        prow = int(comm_size**0.5)
        pcol = comm_size // prow

        self.compute_pencil(0, comm_rank, prow, pcol)
        self.compute_pencil(1, comm_rank, prow, pcol)
        self.compute_pencil(2, comm_rank, prow, pcol)

        return prow, pcol
        
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
       x^
        | 3 4 5
        | 0 1 2
        +------>y
        """

        size = axis_size // ndiv

        start = pgrid_distance * size

        rem = axis_size - ndiv * size
        if (pgrid_distance < rem):
            size += 1

        start += min(rem, pgrid_distance)
        end = start + (size - 1)

        return start, end, size

    def get_pgrid_distance(self, pgrid_length, comm_rank):
        """ Returns the position of rank along an axis of the processor grid. """

        row = comm_rank // self.pcol
        if pgrid_length == self.prow:
            return row
        else:
            return comm_rank - row * self.pcol
        
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
        self.xstart = [0] * 3
        self.xend = [0] * 3
        
        self.ysize = [0] * 3
        self.ystart = [0] * 3
        self.yend = [0] * 3
        
        self.zsize = [0] * 3
        self.zstart = [0] * 3
        self.zend = [0] * 3

        self.xstart[0], self.xend[0], self.xsize[0] = self.get_pencil(self.Nx, 1, 1)
        self.xstart[1], self.xend[1], self.xsize[1] = self.get_pencil_size(self.Ny, self.prow,
                                                                           self.get_pgrid_distance(self.prow,
                                                                                                   comm_rank))
        self.xstart[2], self.xend[2], self.xsize[2] = self.get_pencil_size(self.Nz, self.pcol,
                                                                           self.get_pgrid_distance(self.pcol,
                                                                                                   comm_rank))

        self.ystart[0], self.yend[0], self.ysize[0] = self.get_pencil(self.Nx, 1, 1)
        self.ystart[1], self.yend[1], self.ysize[1] = self.get_pencil_size(self.Ny, self.prow,
                                                                           self.get_pgrid_distance(self.prow,
                                                                                                   comm_rank))
        self.ystart[2], self.yend[2], self.ysize[2] = self.get_pencil_size(self.Nz, self.pcol,
                                                                           self.get_pgrid_distance(self.pcol,
                                                                                                   comm_rank))

        self.zstart[0], self.zend[0], self.zsize[0] = self.get_pencil(self.Nx, 1, 1)
        self.zstart[1], self.zend[1], self.zsize[1] = self.get_pencil_size(self.Ny, self.prow,
                                                                           self.get_pgrid_distance(self.prow,
                                                                                                   comm_rank))
        self.zstart[2], self.zend[2], self.zsize[2] = self.get_pencil_size(self.Nz, self.pcol,
                                                                           self.get_pgrid_distance(self.pcol,
                                                                                                   comm_rank))

