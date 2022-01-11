""" mixing-layer.py

Reference data available from:

https://zenodo.org/record/2577851/files/mixing-layer.tar.xz?download=1
"""

import numpy as np

import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.rc("font", size=11)
plt.rc("axes", axisbelow=True)

import matplotlib.pyplot as plt

from py4incompact3d.postprocess.postprocess import Postprocess
from py4incompact3d.tools.misc import avg_over_axis

T=[182]
RHO1=0.5
RHO2=1.0
NCONTOURS=10
CONTOURS=[RHO1 + i * (RHO2 - RHO1) / (NCONTOURS - 1) for i in range(NCONTOURS)]
CONTOURS[0] *= 1.01
CONTOURS[-1] *= 0.99

def main():

    postprocess = Postprocess('mixing-layer.json')
    
    for t in T:
        msg = "-- timestep" + str(t)
        print(msg)

        postprocess.load(time=t)

        rho = postprocess.fields["rho"].data[t]
        rho = avg_over_axis(postprocess.mesh, rho, 2)

        plt.figure(figsize=(5, 3.5))
        plt.contour(rho.transpose(), CONTOURS, colors="black")
        plt.contour(rhof.transpose(), CONTOURS, colors="red", linestyles="dashed")
        plt.show()
        plt.close()

        postprocess.clear_data()

if __name__ == "__main__":
    main()
