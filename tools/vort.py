"""
.. module:: vort
    :synopsis: Provides function for computing vorticity.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

from Py4Incompact3D.deriv.deriv import deriv

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
                print(i, j)
                S[i][j] = deriv(postprocess, vel, j, t)
        
        # Compute vorticity field

