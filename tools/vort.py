"""
.. module:: vort
    :synopsis: Provides function for computing vorticity.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

from Py4Incompact3D.tools.gradu import calc_gradu, get_gradu_tensor
from Py4Incompact3D.deriv.deriv import deriv
from Py4Incompact3D.postprocess.fields import Field

def get_vort_name(i, j):
    """ Get the name for the specified component of the vorticity tensor.

    :param i: The first component.
    :param j: The second component.

    :type i: int
    :type j: int

    :returns: The name of the specified component of the vorticity tensor.
    :rtype: str
    """

    dirlist = ["x", "y", "z"]
    return "vort" + dirlist[i] + dirlist[j]

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
                field_name = get_vort_name(i, j)
                desc = field_name.replace("vort", "") + "-component of vorticity"
                direction = [i, j]
                prop_dict = {"name":field_name,
                             "description":desc,
                             "properties":{"filename":field_name,
                                           "direction":[0, 0],
                                           "precision":postprocess.fields["ux"].dtype,
                                           "fromfile":False}}
                postprocess.fields[field_name] = Field(prop_dict)
                postprocess.fields[field_name].data[t] = 0.5 * (gradu[i][j] - gradu[j][i])

def get_vort_tensor(postprocess, time=-1):
    r""" Construct the vorticity tensor from the individual components.

    Returns the vorticity tensor in the form

    .. math::
        \frac{1}{2} \left( \frac{\partial u^i}{\partial x^j} - \frac{\partial u^j}{\partial x^i}
        \right)

    where :math:`i` is the first index and :math:`j` the second index, i.e.

    .. math::
        \Omega\left[1\right]left[2] = \frac{1}{2} \left( \frac{\partial v}{\partial z} -
        \frac{\partial w}{\partial z} \right)

    :param postprocess: The post processing object.
    :param time: The time(s) to get the vorticity tensor for.

    :type postprocess: Py4Incompact3D.postprocess.postprocess.Postprocess
    :type time: int or list of int

    :returns: :math`\boldsymbol{\Omega}` a time-keyed dictionary of the vorticity tensor.
    :rtype: dict
    """

    if time == -1:
        load_time = postprocess.fields["ux"].data.keys()
    elif isinstance(time, int):
        load_time = [time]
    elif isinstance(time, list):
        load_time = time
    else:
        raise RuntimeError

    vort = {}
    if not "duxdx" in postprocess.fields.keys():
        calc_gradu(postprocess, load_time)
    gradu = get_gradu_tensor(postprocess, load_time)
    for t in load_time:
        vort[t] = [[0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0]]
        for i in range(3):
            for j in range(3):
                vort[t][i][j] = 0.5 * (gradu[t][i][j] - gradu[t][j][i])

    return vort

