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
from Py4Incompact3D.tools.lockexch import calc_h

####################################################################################################
#
## Globals:
#

# Load timesteps 0-T (add one as python 0-indexed)
T=20
T+=1

# Figure dimensions.
# 5.0 x 3.5 is generally a nice ratio
# Y = (26/164) * X matches Birman2005's plots
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
    for t in range(T):

        # Get h
        postprocess.load(time=t)
        h = calc_h(postprocess)
        postprocess.clear_data()

        # Smooth h

        # Plot
        plt.figure(figsize=(PLTX, PLTY))
        plt.plot(x, h[t], color=CLR)
        plt.xlabel(XLBL)
        plt.ylabel(YLBL)
        filename = NAMERT + "-t" + str(t) + FIGTYPE
        plt.savefig(filename, bbox_inches="tight")
        plt.close()

####################################################################################################
#
# To run from command line: python3 plot_h.py
#                           or set executable and ./plot_h.py
#

if __name__ == "__main__":
    main()
