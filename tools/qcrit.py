"""
.. module:: qcrit
    :synopsis: Computes the q-criterion.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

import numpy as np

from Py4Incompact3D.tools.gradu import calc_gradu, get_gradu_tensor
from Py4Incompact3D.tools.vort import calc_vort
from Py4Incompact3D.deriv.deriv import deriv
from Py4Incompact3D.postprocess.fields import Field

def calc_qcrit(postprocess, time=-1):
    """ Computes the q-criterion of the velocit field, assumes ux uy uz vortx vorty vortz have all
    been loaded/computed.
    
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
        # Get gradu tensor
        if not "duxdx" in postprocess.fields.keys():
            calc_gradu(postprocess, t)
        gradu = get_gradu_tensor(postprocess, t)[t]

        # Get vorticity tensor
        if not "vortxx" in postprocess.fields.keys():
            calc_vort(postprocess, t)

        vort = [[0, 0, 0],
                [0, 0, 0],
                [0, 0, 0]]
        directions = ["x", "y", "z"]
        for i in range(3):
            for j in range(3):
                name = "vort" + directions[i] + directions[j]
                vort[i][j] = postprocess.fields[name].data[t]

        # Construct strain-rate tensor
        S = [[0, 0, 0],
             [0, 0, 0],
             [0, 0, 0]]
        for i in range(3):
            for j in range(3):
                S[i][j] = 0.5 * (gradu[i][j] + gradu[j][i])

        # Compute Q-criterion
        nx = postprocess.fields["ux"].data[t].shape[0]
        ny = postprocess.fields["ux"].data[t].shape[1]
        nz = postprocess.fields["ux"].data[t].shape[2]
        q = np.zeros([nx, ny, nz],
                     dtype=postprocess.fields["ux"].dtype)
        for i in range(3):
            for j in range(3):
                q += 0.5 * (vort[i][j] * vort[j][i] - S[i][j] * S[j][i])

        # Create Q field
        prop_dict = {"name":"Q",
                     "description":"Q-criterion",
                     "properties":{"filename":"Q",
                                   "direction":[-1],
                                   "precision":postprocess.fields["ux"].dtype,
                                   "fromfile":False}}
        postprocess.fields["Q"] = Field(prop_dict)
        postprocess.fields["Q"].data[t] = q

