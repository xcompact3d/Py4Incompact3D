"""
       FILE: meshconvert.py
     AUTHOR: Paul Bartholomew <ptb08@ic.ac.uk>
DESCRIPTION: Converts a dataset from one mesh to another.
"""

import numpy.mgrid
from scipy.interpolate import griddata

def meshconvert(varlist, meshin, meshout, interpolation="cubic"):
    """"""
    
    newvarlist = []
    for var in varlist:
        newvarlist.append(convertone(var, meshin, meshout, interpolation))

    return newvarlist

def convertone(var, meshin, meshout, interpolation):
    """"""

    return griddata(np.mgrid[0:meshin[0], 0:meshin[1], 0:meshin[2]], var,
                    np.mgrid[0:meshout[0], 0:meshout[1], 0:meshout[2]],
                    method = interpolation)
    
    
