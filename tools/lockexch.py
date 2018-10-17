"""
.. module:: lockexch
    :synopsis: Tools for post-processing lock-exchange problems

.. modulauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

import numpy as np
from scipy.signal import argrelextrema

def calc_h(postprocess, field="rho", gamma=0.998, time=-1):
    r""" Calculates the "height" of the gravity-current, assumes name field (default :math:`\rho`)
    is available.

    This is based on the technique proposed in Birman2005 where the height of the gravity current is
    defined as:

    .. math::
        h \left( x \right) = \frac{1}{L_y} \left( \frac{1}{1 - \gamma} \int^{L_y}_0 \overline{\rho}
        \left( x, y \right) dy - \frac{\gamma}{1 - \gamma} \right)

    where :math:`\overline{\rho}` is :math:`\rho` averaged over the z axis.

    :param postprocess: The postprocessing object.
    :param field: The name of the field to calculate the height by
    :param gamma: The density ratio, defined as :math:`\gamma = \frac{\rho_1}{\rho_2},\ 0 \leq
                  \gamma < 1`
    :param time: The time(s) to compute h for, -1 means all times.

    :type postprocess: Py4Incompact3D.postprocess.postprocess.Postprocess
    :type field: str
    :type gamma: float
    :type time: int or list of int

    :returns: h -- a time-keyed dictionary of :math:`h \left( x \right)`
    :rtype: dict

    .. note::
        In the Boussinesq limit, the appropriate field is a concentration field :math:`0 \leq c \leq
        1` for which, set :math:`\gamma = 0`.
    """

    assert(gamma < 1)
    assert(gamma >= 0)

    if time == -1:
        # Load all times
        load_times = postprocess.fields["rho"].data.keys()
    elif isinstance(time, int):
        load_times = [time]
    elif isinstance(time, list):
        load_times = time
    else:
        raise ValueError

    nx = postprocess.mesh.Nx
    ny = postprocess.mesh.Ny
    nz = postprocess.mesh.Nz

    # Setup the integration area arrays
    dy = postprocess.mesh.dy * np.ones(ny)
    if (1 == postprocess.mesh.BCy) or (2 == postprocess.mesh.BCy):
        dy[0] *= 0.5
        dy[-1] *= 0.5
    dz = postprocess.mesh.dz * np.ones(nz)
    if (1 == postprocess.mesh.BCz) or (2 == postprocess.mesh.BCz):
        dz[0] *= 0.5
        dz[-1] *= 0.5

    # Integrate at each time t
    h = {}
    for t in load_times:
        h[t] = [0] * nx
        rho = postprocess.fields["rho"].data[t]

        # Average in z
        h[t] = (rho * dz).sum(axis=2) / postprocess.mesh.Lz

        # Integrate in y
        h[t] *= dy
        h[t] = h[t].sum(axis=1)

        # Account for density ratio
        h[t] -= gamma * postprocess.mesh.Ly
        h[t] /= (1.0 - gamma)

        # Normalise by domain height
        h[t] /= postprocess.mesh.Ly

        # Clip to 0 <= h <= 1
        np.clip(h[t], 0, 1, h[t])

    return h

def get_frontloc_birman(h):
    """ Determines the front locations according to Birman2005. 

    :param h: Time-keyed dictionary of gravity-current height.
    :type h: dict

    :returns: idxr, idxw, idxf: time-keyed dictionaries containing the indices of the front
              locations.
    :rtype: dict, dict, dict

    .. note::
        In the case the front cannot be found, the index will have value None.
    """

    idxr = {}
    idxw = {}
    idxf = {}
    for t in h.keys():
        idxr[t] = None
        idxw[t] = None
        idxf[t] = None

        ## Get indices of local maxima/minima
        maxima = argrelextrema(h[t], np.greater)
        minima = argrelextrema(h[t], np.less)

        idxr[t] = minima[0]
        hw = h[t][0]
        for idx in minima:
            if h[t][idx] < hw:
                idxw[t] = idx
                hw = h[t][idx]
        hf = 0
        for idx in maxima:
            if (idx > idxw) and (h[t][idx] > hf):
                idxf[t] = idx
                hf = h[t][idx]

    return idxr, idxw, idxf
