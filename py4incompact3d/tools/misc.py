"""
.. module:: misc
    :synopsis: For tools which don't fit anywhere else

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

import numpy as np

def moving_avg(x, axis, nsample=16):
    """ Compute the moving average over nsamples.

    :param x: The input array to compute moving average of.
    :param axis: The axis to average over.
    :param nsample: How many samples to compute average over.

    :type x: numpy.ndarray
    :type axis: int
    :type nsample: int

    :returns: xavg - the moving average of x.
    :rtype: numpy.ndarray
    """

    x = np.swapaxes(x, 0, axis)
    xavg = np.zeros(x.shape)

    n = x.shape[0]
    for i in range(n):
        j0 = max(i - nsample // 2, 0)
        j1 = min(i + nsample // 2, n)
        for j in range(j0, j1):
            xavg[i] += x[j]
        xavg[i] /= (j1 - j0)

    x = np.swapaxes(x, 0, axis)
    return np.swapaxes(xavg, 0, axis)

def avg_over_axis(mesh, x, axis):
    """ Averages a field over an axis direction.

    :param mesh: The mesh the field is defined on.
    :param x: The input array to compute average of.
    :param axis: The axis to average over.

    :type mesh: py4incompact3d.postprocess.mesh.Mesh
    :type x: numpy.ndarray
    :type axis: int

    :returns: xavg - the rank-1 field (the averaged axis is collapsed).
    :rtype: numpy.ndarray
    """
    if 0 == axis:
        lx = mesh.Lx
    elif 1 == axis:
        lx = mesh.Ly
    else:
        lx = mesh.Lz

    return int_over_axis(mesh, x, axis) / lx

def int_over_axis(mesh, x, axis):
    """ Integrates a field over an axis direction.

    :param mesh: The mesh the field is defined on.
    :param x: The input array to integrate.
    :param axis: The axis to integrate over.

    :type mesh: py4incompact3d.postprocess.mesh.Mesh
    :type x: numpy.ndarray
    :type axis: int

    :returns: xavg - the rank-1 field (the integrated axis is collapsed).
    :rtype: numpy.ndarray
    """
    if 0 == axis:
        n = mesh.Nx
        dx = mesh.dx * np.ones(n)
        if 0 != mesh.BCx:
            dx[0] *= 0.5
            dx[-1] *= 0.5
    elif 1 == axis:
        n = mesh.Ny
        dx = mesh.dy * np.ones(n)
        if 0 != mesh.BCy:
            dx[0] *= 0.5
            dx[-1] *= 0.5
    else:
        n = mesh.Nz
        dx = mesh.dz * np.ones(n)
        if 0 != mesh.BCz:
            dx[0] *= 0.5
            dx[-1] *= 0.5

    lastaxis = len(x.shape) - 1
    xint = dx * np.swapaxes(x, axis, lastaxis)
    return np.swapaxes(xint, axis, lastaxis).sum(axis=axis)

