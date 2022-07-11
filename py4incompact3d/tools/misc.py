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
    import py4incompact3d
    from py4incompact3d import parallel

    # First transpose data
    if axis > 0:
        # Transpose to Y
        phiX = py4incompact3d.postprocess.Field()
        phiX.new(mesh, pencil=0)
        phiX.data[0][:,:,:] = x[:,:,:]

        phiY = py4incompact3d.postprocess.Field()
        phiY.new(mesh, pencil=1)
        phiY.data[0] = parallel.transpose(phiX.data[0], "xy", phiY.data[0])

        if axis > 1:
            # Transpose to Z
            phiZ = py4incompact3d.postprocess.Field()
            phiZ.new(mesh, pencil=2)
            phiZ.data[0] = parallel.transpose(phiY.data[0], "yz", phiZ.data[0])

            phi = phiZ.data[0]

        else:
            phi = phiY.data[0]
        
    else:
        phi = x

    n = phi.shape[axis]
    if 0 == axis:
        dx = mesh.dx
        bc = mesh.BCx
    elif 1 == axis:
        dx = mesh.dy
        bc = mesh.BCy
    else:
        dx = mesh.dz
        bc = mesh.BCz

    dx = dx * np.ones(n)
    if bc != 0:
        dx[0] *= 0.5
        dx[-1] *= 0.5

    lastaxis = len(phi.shape) - 1
    phi_int = np.swapaxes(dx * np.swapaxes(phi, axis, lastaxis),
                       axis,
                       lastaxis)
    s = phi_int.sum(axis=axis)

    if axis > 0:
        if axis > 1:
            # Transpose from Z
            for i in range(s.shape[0]):
                for j in range(s.shape[1]):
                        phiZ.data[0][i,j,:] = s[i,j]

            phiY.data[0] = parallel.transpose(phiZ.data[0], "zy", phiY.data[0])

        else:
            for i in range(s.shape[0]):
                for k in range(s.shape[1]):
                    phiY.data[0][i,:,k] = s[i,k]

        # Transpose from Y
        phiX.data[0] = parallel.transpose(phiY.data[0], "yx", phiX.data[0])
        phi = phiX.data[0]

    else:
        for j in range(s.shape[0]):
            for k in range(s.shape[1]):
                phi[:,j,k] = s[j,k]

    return phi

