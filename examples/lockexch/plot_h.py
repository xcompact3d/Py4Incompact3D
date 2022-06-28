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

from py4incompact3d.postprocess.postprocess import Postprocess
from py4incompact3d.tools.lockexch import calc_h, get_frontidx_birman
from py4incompact3d.tools.misc import moving_avg

####################################################################################################
#
## Globals:
#

# What density ratio?
GAMMA=0.2

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
PLOTH=True # Do we even want to plot h?
FIGTYPE=".eps"
XLBL=r"$x$"
YLBL=r"$h$"
CLR="black"
NAMERT="height"
MARKFRNT=True # Should the fronts be indicated on the plots?

####################################################################################################
#
## Script:
#
#  What follows should not need to be changed.
#

def plot_h(postprocess, plot=True, mark=True):

    ## Print header
    print("+" + "=" * 70 + "+")
    header = "| Computing "
    if plot:
        header += "(and plotting) "
    header += "h"
    header += " " * (71 - len(header)) + "|"
    print(header)
    print("+" + "-" * 70 + "+")

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
    if mark:
        idxr, idxw, idxf = {}, {}, {}
    for t in range(T):

        msg = "| Processing time: " + str(t)
        msg += " " * (71 - len(msg)) + "|"
        print(msg)

        # Get h
        postprocess.load(time=t)
        h[t] = calc_h(postprocess, "rho", GAMMA)[t]
        postprocess.clear_data()

        # Smooth h using moving average
        h[t] = np.array(moving_avg(h[t], 0, NSMPL))

        ## Get fronts
        if mark:
            idxr, idxw, idxf = get_frontidx_birman(h)

        if plot:
            # Plot
            plt.figure(figsize=(PLTX, PLTY))
            plt.plot(x, h[t], color=CLR)
            if mark:
                if (idxr[t]):
                    plt.plot(x[idxr[t]], h[t][idxr[t]], label=r"$x_r$",
                             color=CLR,
                             marker="o",
                             ls="")
                if (idxw[t]):
                    plt.plot(x[idxw[t]], h[t][idxw[t]], label=r"$x_w$",
                             color=CLR,
                             marker="d",
                             ls="")
                if (idxf[t]):
                    plt.plot(x[idxf[t]], h[t][idxf[t]], label=r"$x_f$",
                             color=CLR,
                             marker="*",
                             ls="")
            plt.xlabel(XLBL)
            plt.ylabel(YLBL)
            filename = NAMERT + "-t" + str(t) + FIGTYPE
            plt.savefig(filename, bbox_inches="tight")
            plt.legend()
            plt.close()

    print("+" + "-" * 70 + "+")
    return h

def plot_vel(postprocess, h):

    ## Print header
    print("+" + "=" * 70 + "+")
    header = "| Computing velocities "
    header += " " * (71 - len(header)) + "|"
    print(header)
    print("+" + "-" * 70 + "+")

    ## Get fronts
    idxr, idxw, idxf = get_frontidx_birman(h)

    ## Compute velocity
    ur, uf = [], []
    n = len(h.keys())
    keys_ordered = sorted(h.keys())
    for t in range(1, n):

        msg = "| Processing time: " + str(t)
        msg += " " * (71 - len(msg)) + "|"
        print(msg)

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

    print("+" + "-" * 70 + "+")

    ## Plot
    plt.plot(np.absolute(ur), label="ur")
    plt.plot(np.absolute(uf), label="uf")
    plt.legend()
    # plt.xlim(xmin=0.5)
    # plt.ylim((0, 1))
    plt.show()

def main():

    ## Initialise the postprocessing object
    postprocess = Postprocess("input.json")

    ## Compute and plot h
    h = plot_h(postprocess, PLOTH, MARKFRNT)

    ## Compute and plot front velocities
    plot_vel(postprocess, h)

####################################################################################################
#
# To run from command line: python3 plot_h.py
#                           or set executable and ./plot_h.py
#

if __name__ == "__main__":
    main()
