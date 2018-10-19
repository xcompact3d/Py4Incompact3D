#!/usr/bin/python3
"""
.. module:: plot_h
    :synopsis: Calculates the height of the gravity current and plots as a function of time.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

####################################################################################################
#
## Imports:
#

import numpy as np
import matplotlib.pyplot as plt

from Py4Incompact3D.postprocess.postprocess import Postprocess
from Py4Incompact3D.tools.lockexch import calc_h, get_frontidx_birman
from Py4Incompact3D.tools.misc import moving_avg

####################################################################################################
#
## Globals:
#

# Load timesteps 0-T (add one as python 0-indexed)
T=20
T+=1

# Moving average size
NSMPL=16

# Figure dimensions.
# 5.0 x 3.5 is generally a nice ratio
# Y = (26/164) * X matches aspect ratio of Birman2005's plots
PLTX=5.0
PLTY=(26. / 164.) * PLTX

# Plot customisation
FIGTYPE=".eps"
XLBL=r"$x$"
YLBL=r"$h$"
CLR="black"
NAMERT="height"

####################################################################################################
#
## Script:
#
#  What follows should not need to be changed.
#

def main():

    ## Initialise the postprocessing object
    postprocess = Postprocess("input.json")

    ## Setup the x-array
    #
    #  np.arange() does not include the stop value, therefore we have to append it.
    #
    x = np.append(np.arange(start=0,
                            stop=postprocess.mesh.Lx,
                            step=postprocess.mesh.dx),
                  postprocess.mesh.Lx)

    ## Loop over time, compute h and plot
    h = {}
    for t in range(T):

        # Get h
        postprocess.load(time=t)
        h[t] = calc_h(postprocess)[t]
        postprocess.clear_data()

        # Smooth h using moving average
        h[t] = np.array(moving_avg(h[t], 0, NSMPL))

        # Plot
        plt.figure(figsize=(PLTX, PLTY))
        plt.plot(x, h[t], color=CLR)
        plt.xlabel(XLBL)
        plt.ylabel(YLBL)
        filename = NAMERT + "-t" + str(t) + FIGTYPE
        plt.savefig(filename, bbox_inches="tight")
        plt.close()

    ## Get fronts and compute velocity
    idxr, idxw, idxf = get_frontidx_birman(h)
    ur, uf = [], []
    n = len(h.keys())
    keys_ordered = sorted(h,keys())
    for t in range(1, n):
        idxr0 = idxr[keys_ordered[t - 1]]
        if idxr0:
            xr0 = idxr0 * postprocess.mesh.dx
        else:
            xr0 = 0

        idxrn = idxr[keys_ordered[t]]
        if idxrn:
            xrn = idxrn * postprocess.mesh.dx
        else:
            xrn = 0

        ur.append((xrn - xr0) / float(keys_ordered[t] - keys_ordered[t - 1]))

        idxf0 = idxf[keys_ordered[t - 1]]
        if idxf0:
            xf0 = idxf0 * postprocess.mesh.dx
        else:
            xf0 = 0

        idxfn = idxf[keys_ordered[t]]
        if idxfn:
            xfn = idxfn * postprocess.mesh.dx
        else:
            xfn = 0

        uf.append((xfn - xf0) / float(keys_ordered[t] - keys_ordered[t - 1]))
    plt.plot(ur)
    plt.plot(uf)
    plt.show()

####################################################################################################
#
# To run from command line: python3 plot_h.py
#                           or set executable and ./plot_h.py
#

if __name__ == "__main__":
    main()
