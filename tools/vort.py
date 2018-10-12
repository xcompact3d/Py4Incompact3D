"""
.. module:: vort
    :synopsis: Provides function for computing vorticity.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

import numpy as np

from Py4Incompact3D.deriv.deriv import deriv
from Py4Incompact3D.postprocess.fields import Field

def calc_vort(postprocess, time=-1):
    """ Computes the vorticity of the velocity field, assumes ux uy and uz have all been loaded. 

    :param postprocess: The postprocessing object. 
    :param time: The time to compute vorticity at, -1 means all times.

    :type postprocess: Py4Incompact3D.postprocess.postprocess.Postprocess
    :type time: int or list of int
    """

    if time == -1:
        time = postprocess.fields["ux"].data.keys()
    elif isinstance(time, int):
        time = [time]
    elif isinstance(time, list):
        pass
    else:
        raise RuntimeError

    for t in time:
        # Compute strain-rate tensor
        S = [[0, 0, 0],
             [0, 0, 0],
             [0, 0, 0]]
        for vel in ["ux", "uy", "uz"]:
            i = postprocess.fields[vel].direction[0]
            for j in range(3):
                S[i][j] = 0.5 * deriv(postprocess, vel, j, t)
        
        # Compute vorticity tensor
        vortxx = (S[0][0] - S[0][0])
        vortxy = (S[0][1] - S[1][0])
        vortxz = (S[0][2] - S[2][0])
        vortyx = -vortxy
        vortyy = (S[1][1] - S[1][1])
        vortyz = (S[1][2] - S[2][1])
        vortzx = -vortxz
        vortzy = -vortyz
        vortzz = (S[2][2] - S[2][2])

        prop_dict = {"name":"vortxx",
                     "description":"xx-component of vorticity",
                     "properties":{"filename":"vortxx",
                                   "direction":[0, 0],
                                   "precision":postprocess.fields["ux"].dtype,
                                   "fromfile":False}}
        postprocess.fields["vortxx"] = Field(prop_dict)
        postprocess.fields["vortxx"].data[t] = vortxx
        prop_dict = {"name":"vortxy",
                     "description":"xy-component of vorticity",
                     "properties":{"filename":"vortxy",
                                   "direction":[0, 1],
                                   "precision":postprocess.fields["ux"].dtype,
                                   "fromfile":False}}
        postprocess.fields["vortxy"] = Field(prop_dict)
        postprocess.fields["vortxy"].data[t] = vortxy
        prop_dict = {"name":"vortxz",
                     "description":"xz-component of vorticity",
                     "properties":{"filename":"vortxz",
                                   "direction":[0, 2],
                                   "precision":postprocess.fields["ux"].dtype,
                                   "fromfile":False}}
        postprocess.fields["vortxz"] = Field(prop_dict)
        postprocess.fields["vortxz"].data[t] = vortxz

        prop_dict = {"name":"vortyx",
                     "description":"xx-component of vorticity",
                     "properties":{"filename":"vortyx",
                                   "direction":[1, 0],
                                   "precision":postprocess.fields["ux"].dtype,
                                   "fromfile":False}}
        postprocess.fields["vortyx"] = Field(prop_dict)
        postprocess.fields["vortyx"].data[t] = vortyx
        prop_dict = {"name":"vortyy",
                     "description":"xy-component of vorticity",
                     "properties":{"filename":"vortyy",
                                   "direction":[1, 1],
                                   "precision":postprocess.fields["ux"].dtype,
                                   "fromfile":False}}
        postprocess.fields["vortyy"] = Field(prop_dict)
        postprocess.fields["vortyy"].data[t] = vortyy
        prop_dict = {"name":"vortyz",
                     "description":"xz-component of vorticity",
                     "properties":{"filename":"vortyz",
                                   "direction":[1, 2],
                                   "precision":postprocess.fields["ux"].dtype,
                                   "fromfile":False}}
        postprocess.fields["vortyz"] = Field(prop_dict)
        postprocess.fields["vortyz"].data[t] = vortyz

        prop_dict = {"name":"vortzx",
                     "description":"xx-component of vorticity",
                     "properties":{"filename":"vortzx",
                                   "direction":[2, 0],
                                   "precision":postprocess.fields["ux"].dtype,
                                   "fromfile":False}}
        postprocess.fields["vortzx"] = Field(prop_dict)
        postprocess.fields["vortzx"].data[t] = vortzx
        prop_dict = {"name":"vortzy",
                     "description":"xy-component of vorticity",
                     "properties":{"filename":"vortzy",
                                   "direction":[2, 1],
                                   "precision":postprocess.fields["ux"].dtype,
                                   "fromfile":False}}
        postprocess.fields["vortzy"] = Field(prop_dict)
        postprocess.fields["vortzy"].data[t] = vortzy
        prop_dict = {"name":"vortzz",
                     "description":"xz-component of vorticity",
                     "properties":{"filename":"vortzz",
                                   "direction":[2, 2],
                                   "precision":postprocess.fields["ux"].dtype,
                                   "fromfile":False}}
        postprocess.fields["vortzz"] = Field(prop_dict)
        postprocess.fields["vortzz"].data[t] = vortzz
