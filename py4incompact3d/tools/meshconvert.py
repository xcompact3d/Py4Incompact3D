"""
       FILE: meshconvert.py
     AUTHOR: Paul Bartholomew <ptb08@ic.ac.uk>
DESCRIPTION: Converts a dataset from one mesh to another.
"""

import numpy as np
import numpy.mgrid
from scipy.interpolate import griddata

def meshconvert(meshin, meshout, infile, outfile, interpolation="cubic", dtype=np.float64):
    """"""

    numin = np.prod(np.array(meshin))
    datain = np.fromfile(infile, dtype=dtype)

    nfields = int(len(datain) / numin)
    assert((nfields * numin) == len(datain))

    numout = np.prod(np.array(meshout))
    dataout = np.zeros(nfields * numout)

    for f in range(nfields):
        dataout[f*numout:(f+1)*numout] = to_fortran(convertone(np.reshape(datain[f*numin:(f+1)*numin],
                                                                          (meshin[0], meshin[1], meshin[2]),
                                                                          "F"),
                                                               meshin,
                                                               meshout,
                                                               interpolation))

    np.tofile(outfile,  dataout)

def to_fortran(arr):
    return np.swapaxes(arr, 0, 2)

def convertone(var, meshin, meshout, interpolation):
    """"""

    return griddata(np.mgrid[0:meshin[0], 0:meshin[1], 0:meshin[2]],
                    var,
                    np.mgrid[0:meshout[0], 0:meshout[1], 0:meshout[2]],
                    method = interpolation)
    
    
