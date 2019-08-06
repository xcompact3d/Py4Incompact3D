"""
 probe.py
"""

import numpy as np
from scipy.io import FortranFile

class Probe():

    def __init__(self, **kwargs):

        self.dtype = np.float64
        
        for arg, val in kwargs.items():
            if arg == "file":
                self.file = val
            elif arg == "variables":
                self.variables = val
            elif arg == "dtype":
                if val == "single":
                    self.dtype = np.float32
                else:
                    self.dtype = np.float64
            elif arg == "freq":
                self.freq = val

    def read(self, nx, ny, nz):

        with open(self.file, "rb") as bindat:
            fldat = np.fromfile(bindat, self.dtype)

        nfields = len(self.variables)
        if self.freq[0]:
            nx = nx // self.freq[0] + 1
        else:
            nx = 1
        if self.freq[1]:
            ny = ny // self.freq[1] + 1
        else:
            ny = 1
        if self.freq[2]:
            nz = nz // self.freq[2] + 1
        else:
            nz = 1

        ntime = len(fldat) // (nfields * nx * ny * nz + 1) # Note that Fortran adds a byte for new-line

        bindat = FortranFile(self.file, "r")
        fldat = []
        for t in range(ntime):
            fldat.append(bindat.read_reals(dtype=np.float64).reshape(nfields, nx, ny, nz))
        bindat.close()

        fldict = {}
        for field in range(nfields):
            fldict[self.variables[field]] = []

        for t in range(ntime):
            for f in range(nfields):
                fldict[self.variables[f]].append(fldat[t][f])

        return fldict
        
