# Copyright 2018 Georgios Deskos

# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.

import py4incompact3d

from warnings import warn

import math

import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import griddata

if py4incompact3d.HAVE_DECOMP2D:
    import decomp2d

class Mesh():
    """
    Mesh is a model object representing

    """

    def __init__(self, *arg, **kwargs):

        super().__init__()

        # Set default values
        self.Nx = 0; self.Ny = 0; self.Nz = 0
        self.Lx = 0; self.Ly = 0; self.Lz = 0
        self.BCx = -1; self.BCy = -1; self.BCz = -1

        self.stretched = False
        self.beta = 0
        self.yp = None

        # Figure out how to do initialisation
        if len(arg) == 1:
            warn("You are using an old-style initialisation, the future is dynamic!", DeprecationWarning)
            self._init_fromjson(arg[0])
        else:
            self._init_new(*arg, **kwargs)

        # Finish off initialisation
        if self.Nx * self.Ny * self.Nz == 0:
            raise RuntimeError("Need to set Nx, Ny and Nz")
        elif self.Lx * self.Ly * self.Lz == 0:
            raise RuntimeError("Need to set Lx, Ly and Lz")
        elif (self.BCx == -1) or (self.BCy == -1) or (self.BCz == -1):
            raise RuntimeError("Need to set boundary conditions!")

        self.__init_decomposition()
        
    def _init_new(self, *args, **kwargs):

        for arg, val in kwargs.items():
            if arg == "n":
                self.Nx = val[0]
                self.Ny = val[1]
                self.Nz = val[2]
            elif arg == "l":
                self.Lx = val[0]
                self.Ly = val[1]
                self.Lz = val[2]
            elif arg == "bc":
                self.BCx = val[0]
                self.BCy = val[1]
                self.BCz = val[2]
            elif arg == "beta":
                self.beta = val
            elif arg == "stretched":
                self.stretched = val
            elif arg == "yp":
                self.yp = val
            else:
                warn("Unrecognised input to mesh: %s" % arg)

    def _init_fromjson(self, instance_dictionary):

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
        try:
            self.stretched = properties["stretched"]
            self.beta = properties["beta"]
        except:
            self.stretched = False
            self.beta = 0

        # Once we know the mesh layout we can set the derivative variables
        self.compute_derivvars()

        if self.stretched:
            try:
                self.yp = properties["yp"]
                with open(self.yp, "r") as ypfile:
                    j = 0
                    self.yp = np.zeros(self.Ny)
                    for row in ypfile:
                        self.yp[j] = float(row)
                        j += 1

                yp, yeta = self.calc_yp()
            except:
                self.yp, yeta = self.calc_yp()

            self.ppy = self.calc_ppy(yeta)
        else:
            self.yp = None

    def __init_decomposition(self):

        if py4incompact3d.HAVE_DECOMP2D:
            # Initialise 2decomp&fft
            prow = 0
            pcol = 0
            decomp2d.decomp4py.init_decomp4py(self.Nx, self.Ny, self.Nz, 0, 0)

            # Create local size arrays (one per pencil)
            self.NxLocal = [0, 0, 0]
            self.NyLocal = [0, 0, 0]
            self.NzLocal = [0, 0, 0]

            # Popoulate local size arrays
            for ax in ["x", "y", "z"]:
                self.___populate_local_sizes(ax)

            self.NxStart = [0, 0, 0]
            self.NyStart = [0, 0, 0]
            self.NzStart = [0, 0, 0]

            for ax in ["x", "y", "z"]:
                self.___populate_local_starts(ax)

        else:

            self.NxLocal = [self.Nx, self.Nx, self.Nx]
            self.NyLocal = [self.Ny, self.Ny, self.Ny]
            self.NzLocal = [self.Nz, self.Nz, self.Nz]

            self.NxStart = [0, 0, 0]
            self.NyStart = [0, 0, 0]
            self.NzStart = [0, 0, 0]
        
    def ___populate_local_sizes(self, ax):
        
        d = {"x":0, "y":1, "z":2}

        n = decomp2d.decomp4py.get_grid_size(ax)
        self.NxLocal[d[ax]] = n[0]
        self.NyLocal[d[ax]] = n[1]
        self.NzLocal[d[ax]] = n[2]
        
    def ___populate_local_starts(self, ax):
        
        d = {"x":0, "y":1, "z":2}

        n = decomp2d.decomp4py.get_grid_start(ax)
        self.NxStart[d[ax]] = n[0] - 1 # Use zero indexing 
        self.NyStart[d[ax]] = n[1] - 1 # Use zero indexing 
        self.NzStart[d[ax]] = n[2] - 1 # Use zero indexing 

    def global_index(self, pencil, i, j, k):
        """ Returns the global index of a point in a given pencil.

        Note - uses C/Python (i.e. row-major) ordering!
        """

        ig = i + self.NxStart[pencil]
        jg = j + self.NyStart[pencil]
        kg = k + self.NzStart[pencil]

        idxg = ig * self.Ny * self.Nz + jg * self.Nz + kg

        return idxg

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
            self.Nxm = self.Nx
        else:
            self.dx = self.Lx / float(self.Nx-1)
            self.Nxm = self.Nx - 1
            
        if (self.BCy==0):
            self.dy = self.Ly / float(self.Ny) # XXX This will not be correct for stretched grids
            self.Nym = self.Ny
        else:
            self.dy = self.Ly / float(self.Ny-1) # XXX This will not be correct for stretched grids
            self.Nym = self.Ny - 1
            
        if (self.BCz==0):
            self.dz = self.Lz / float(self.Nz)
            self.Nzm = self.Nz
        else:
            self.dz=self.Lz/float(self.Nz-1)
            self.Nzm = self.Nz - 1

        self.alpha = 1.0 / 3.0
        self.a = 14.0 / 9.0
        self.b = 1.0 / 9.0

    def calc_yp(self):

        self.compute_derivvars()
        
        yinf = -self.Ly / 2.0
        den = 2.0 * self.beta * yinf
        xnum = -(yinf + math.sqrt((math.pi * self.beta)**2 + yinf**2))
        alpha = abs(xnum / den)
        xcx = 1.0 / self.beta / alpha

        yp = np.zeros(self.Ny)
        yeta = np.zeros(self.Ny)
        
        if (alpha != 0.0):
            yp[0] = 0.0
            yeta[0] = -0.5
            for j in range(1, self.Ny):
                if (self.stretched == 1):
                    yeta[j] = (j - 1.0) / self.Nym
                elif (self.stretched == 2):
                    yeta[j] = (j - 1.0) / self.Nym - 0.5
                else:
                    yeta[j] = (j - 1.0) * 0.5 / self.Ny - 0.5

                den1 = math.sqrt(alpha * self.beta + 1)
                xnum = den1 / math.sqrt(alpha / math.pi) / math.sqrt(self.beta) \
                       / math.sqrt(math.pi)
                den = 2.0 * math.sqrt(alpha / math.pi) * math.sqrt(self.beta) \
                      * math.pi**1.5
                den3 = (math.sin(math.pi * yeta[j]))**2 / self.beta / math.pi \
                       + alpha / math.pi
                den4 = 2.0 * alpha * self.beta - math.cos(2.0 * math.pi * yeta[j]) + 1.0
                xnum1 = math.atan(xnum * math.tan(math.pi * yeta[j])) \
                         * den4 / den1 / den3 / den
                cst = math.sqrt(self.beta) * math.pi \
                      / (2.0 * math.sqrt(alpha) * math.sqrt(alpha * self.beta + 1.0))

                if (yeta[j] < 0.5):
                    if (self.stretched == 1):
                        yp[j] = xnum1 - cst - yinf
                    elif (self.stretched == 2):
                        yp[j] = xnum1 - cst + self.Ly
                    else:
                        yp[j] = 2 * (xnum1 - cst + self.Ly)
                elif (yeta[j] > 0.5):
                    if (self.stretched == 1):
                        yp[j] = xnum1 + cst - yinf
                    elif (self.stretched == 2):
                        yp[j] = xnum1 + cst + self.Ly
                    else:
                        yp[j] = 2 * (xnum1 - cst + self.Ly)
                else:
                    if (self.stretched == 1):
                        yp[j] = -yinf
                    elif (self.stretched == 2):
                        yp[j] = self.Ly
                    else:
                        yp[j] = 2 * self.Ly
        else:
            yp[0] = -1e10
            for j in range(1, self.Ny):
                yeta[j] = (j - 1.0) / float(self.Ny)
                yp[j] = -self.beta * math.cos(math.pi * yeta[j]) / math.sin(yeta[j] * math.pi)

        return yp, yeta

    def calc_ppy(self, yeta):

        ppy = np.zeros(self.Ny)
        alpha = self.calc_alpha()
        
        if (self.stretched == 3):
            yetai = self.calc_yetai(alpha)
            for j in range(self.Ny):
                ppy[j] = self.Ly * (alpha / math.pi \
                                    + (1.0 / math.pi / self.beta) * (math.sin(math.pi * yetai[j]))**2)
        else:
            for j in range(self.Ny):
                ppy[j] = self.Ly * (alpha / math.pi \
                                    + (1.0 / math.pi / self.beta) * (math.sin(math.pi * yeta[j]))**2)

        return ppy

    def calc_alpha(self):

        yinf = -self.Ly / 2.0
        den = 2.0 * self.beta * yinf
        xnum = -(yinf + math.sqrt(math.pi**2 * self.beta**2 + yinf**2))

        return abs(xnum / den)

    def calc_yetai(self, alpha):

        yetai = np.zeros(self.Ny)
        
        if (alpha != 0.0):
            for j in range(self.Ny):
                yetai[j] = (j - 0.5) * (0.5 / self.Nym) - 0.5

        else:
            for j in range(self.Ny):
                yetai[j] = (j - 1.0) / self.Ny

        return yetai

    def get_grid(self):
        """ Return the x,y,z arrays that describe the mesh. """

        x, y, z = np.zeros(self.Nx), np.zeros(self.Ny), np.zeros(self.Nz)

        for i in range(self.Nx):
            x[i] = i * self.dx

        use_yp = True

        if self.yp is None:
            use_yp = False
        elif isinstance(self.yp, np.ndarray):
            if not self.yp.any():
                use_yp = False
            
        if use_yp:
            for j in range(self.Ny):
                y[j] = self.yp[j]
        else:
            for j in range(self.Ny):
                y[j] = j * self.dy
                    
        for k in range(self.Nz):
            z[k] = k * self.dz

        return x, y, z
