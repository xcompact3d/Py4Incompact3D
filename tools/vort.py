"""
.. module:: vort
    :synopsis: Provides function for computing vorticity.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

from Py4Incompact3D.tools.gradu import calc_gradu, get_gradu_tensor
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

    grad_list = ["x", "y", "z"]
    for t in time:
        # Get gradu tensor
        if not "duxdx" in postprocess.fields.keys():
            calc_gradu(postprocess, t)
        gradu = get_gradu_tensor(postprocess, t)[t]
        
        # Compute vorticity tensor
        for i in range(3):
            for j in range(3):
                field_name = "vort" + grad_list[i] + grad_list[j]
                desc = grad_list[i] + grad_list[j] + "-component of vorticity"
                direction = [i, j]
                prop_dict = {"name":field_name,
                             "description":desc,
                             "properties":{"filename":field_name,
                                           "direction":[0, 0],
                                           "precision":postprocess.fields["ux"].dtype,
                                           "fromfile":False}}
                postprocess.fields[field_name] = Field(prop_dict)
                postprocess.fields[field_name].data[t] = 0.5 * (gradu[i][j] - gradu[j][i])

