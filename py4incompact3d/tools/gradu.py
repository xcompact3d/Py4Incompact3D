"""
.. module:: gradu
    :synopsis: Computes gradient of velocity field.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

from py4incompact3d.deriv.deriv import deriv
from py4incompact3d.postprocess.fields import Field

def get_gradu_name(i, j):
    """ Determine the name for a component of the velocity gradient tensor.

    :param i: The velocity component.
    :param j: The gradient component.

    :type i: int
    :type j: int

    :returns: The name of the specified component of the velocity gradient tensor.
    :rtype: str
    """

    vel_list = ["ux", "uy", "uz"]
    grad_list = ["x", "y", "z"]

    return "d" + vel_list[i] + "d" + grad_list[j]

def calc_gradu(postprocess, time=-1):
    """ Computes the gradient of the velocity field, assumes ux uy uz have all been loaded.

    :param postprocess: The postprocessing object.
    :param time: The time to compute vorticity at, -1 means all times.

    :type postprocess: py4incompact3d.postprocess.postprocess.Postprocess
    :type time: int or list of int
    """

    if time == -1:
        load_time = postprocess.fields["ux"].data.keys()
    elif isinstance(time, int):
        load_time = [time]
    elif isinstance(time, list):
        load_time = time
    else:
        raise RuntimeError

    for t in load_time:
        # Compute velocity-gradient tensor
        for i in range(3):
            for j in range(3):
                field_name = get_gradu_name(i, j)
                desc = field_name
                prop_dict = {"name":field_name,
                             "description":desc,
                             "properties":{"filename":field_name,
                                           "direction":[i, j],
                                           "precision":postprocess.fields["ux"].dtype,
                                           "fromfile":False}}
                postprocess.fields[field_name] = Field(prop_dict)
                postprocess.fields[field_name].data[t] = deriv(postprocess, ["ux", "uy", "uz"][i], j, t)

def get_gradu_tensor(postprocess, time=-1):
    r""" Construct the gradient tensor from the individual components.

    Returns the gradient tensor in the form

    .. math::
        \frac{\partial u^i}{\partial x^j}

    where :math:`i` is the first index and :math:`j` the second index, i.e.

    .. math::
        grad\left( \boldsymbol{u} \right)\left[ 1 \right]\left[ 2 \right] = \frac{\partial
        v}{\partial z}

    :param postprocess: The post processing object
    :param time: The time(s) to get the gradient tensor for.

    :type postprocess: py4incompact3d.postprocess.postprocess.Postprocess
    :type time: int or list of int

    :returns: :math`\boldsymbol{\nabla}\boldsymbol{u}` a time-keyed dictionary of the gradient tensor.
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

    gradu = {}
    for t in load_time:
        gradu[t] = [[0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0]]
        for i in range(3):
            for j in range(3):
                gradu[t][i][j] = postprocess.fields[get_gradu_name(i, j)].data[t]

    return gradu

