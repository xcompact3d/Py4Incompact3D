"""
.. module:: gradu
    :synopsis: Computes gradient of velocity field.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

from Py4Incompact3D.deriv.deriv import deriv
from Py4Incompact3D.postprocess.fields import Field

def calc_gradu(postprocess, time=-1):
    """ Computes the gradient of the velocity field, assumes ux uy uz have all been loaded.

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

    vel_list = ["ux", "uy", "uz"]
    grad_list = ["x", "y", "z"]
    for t in time:
        # Compute velocity-gradient tensor
        for vel in vel_list:
            i = postprocess.fields[vel].direction[0]
            for j in range(3):
                field_name = "d" + vel + d + grad_list[j]
                desc = "Gradient of " + vel + " wrt " + grad_list[j]
                prop_dict = {"name":field_name,
                             "description":desc,
                             "properties":{"filename":field_name,
                                           "direction":[i, j],
                                           "precision":postprocess.fields["ux"].dtype,
                                           "fromfile":False}}
                postprocess.fields[field_name] = Field(prop_dict)
                postprocess.fields[field_name].data[t] = deriv(postprocess, vel, j, t)

