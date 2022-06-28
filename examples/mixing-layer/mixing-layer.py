""" mixing-layer.py

Reference data available from:

https://zenodo.org/record/2577851/files/mixing-layer.tar.xz?download=1
"""

from mpi4py import MPI

import numpy as np

import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
#plt.rc("text", usetex=True)
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

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    postprocess = Postprocess('mixing-layer.json')
    
    for t in T:
        msg = "-- timestep" + str(t)
        print(msg)

        postprocess.load(time=t)

        rho = postprocess.fields["rho"].data[t]
        rho = avg_over_axis(postprocess.mesh, rho, 2)

        # Extract planar data
        rho_plane = np.zeros((postprocess.mesh.Nx, postprocess.mesh.Ny))
        if (postprocess.mesh.NzStart[0] == 0):
            for i in range(postprocess.mesh.NxLocal[0]):
                for j in range(postprocess.mesh.NyLocal[0]):

                    ii = i + postprocess.mesh.NxStart[0]
                    jj = j + postprocess.mesh.NyStart[0]

                    rho_plane[ii][jj] = rho[i][j][0]

        rho_plane_out = np.zeros(rho_plane.shape)
        comm.Reduce(rho_plane, rho_plane_out, op=MPI.SUM)

        if (rank == 0):
            plt.figure(figsize=(5, 3.5))
            plt.contour(rho_plane_out.transpose(), CONTOURS, colors="black")
            #plt.contour(rhof.transpose(), CONTOURS, colors="red", linestyles="dashed")
            plt.savefig("foo.png")
            plt.close()

        postprocess.clear_data()

if __name__ == "__main__":
    main()
