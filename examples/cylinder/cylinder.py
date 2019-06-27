"""
       FILE: cylinder.py
     AUTHOR: Paul Bartholomew
DESCRIPTION: Post processes the cylinder case with reference data.
"""

import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.rc("font", size=11)

from Py4Incompact3D.postprocess.postprocess import Postprocess
from Py4Incompact3D.tools.misc import avg_over_axis

REFPATH="/home/paul/DATA/benchmarking/cylinder/"
REFUMEAN="umean.csv"
REFUUMEAN="uumean.csv"

FIGX=5.0
FIGY=3.5

INPUT="input.json"

NTIME=120000

# Cylinder centre
X0=5.0
Y0=6.0

def main():

    postprocess = Postprocess(INPUT)
    mesh = postprocess.mesh
    t = "0160000"
    postprocess.load(time=[t])

    umean1 = postprocess.fields["umean1"].data[t] / float(NTIME)
    uumean1 = postprocess.fields["uumean1"].data[t] / float(NTIME)
    umean2 = postprocess.fields["umean2"].data[t] / float(NTIME)
    uumean2 = postprocess.fields["uumean2"].data[t] / float(NTIME)

    umean1 = avg_over_axis(mesh, umean1, 2)
    uumean1 = avg_over_axis(mesh, uumean1, 2)
    umean2 = avg_over_axis(mesh, umean2, 2)
    uumean2 = avg_over_axis(mesh, uumean2, 2)

    uprime1 = uumean1 - umean1**2
    uprime2 = uumean2 - umean2**2

    # Plot vs. reference data
    plt.figure(figsize=(FIGX, FIGY))
    
    umean_ref = (read_refdat(REFPATH+REFUMEAN))
    for key in umean_ref:
        y = get_nearest_profile(umean1, mesh, get_axialdist(key) + X0)
        y = get_subrange(y, mesh.dy, Y0, -2, +2)
        y = zero_offset(y)
        y = offset_by_axialdist(y, key)
        x = [j * mesh.dy - 0.5 * mesh.Ly for j in range(mesh.Ny)]
        x = get_subrange(x, mesh.dy, Y0, -2, +2)
        plt.plot(x, y, color="black",
                 ls="-")
        
        y = get_nearest_profile(umean2, mesh, get_axialdist(key) + X0)
        y = get_subrange(y, mesh.dy, Y0, -2, +2)
        y = zero_offset(y)
        y = offset_by_axialdist(y, key)
        x = [j * mesh.dy - 0.5 * mesh.Ly for j in range(mesh.Ny)]
        x = get_subrange(x, mesh.dy, Y0, -2, +2)
        plt.plot(x, y, color="red",
                 ls="--")
        
        x, y = get_xydata(umean_ref, key)
        y = zero_offset(y)
        y = offset_by_axialdist(y, key)
        plt.plot(x, y,
                 ls="", color="black",
                 marker="d")

    plt.xlim((-2, 2))
    plt.ylim((-3.5, -1))
    plt.xlabel(r"$y$")
    plt.ylabel(r"$\langle u \rangle$")
    plt.savefig("umean.eps", bbox_inches="tight")

    plt.figure(figsize=(FIGX, FIGY))

    uumean_ref = (read_refdat(REFPATH+REFUUMEAN))
    for key in uumean_ref:
        y = get_nearest_profile(uprime1, mesh, get_axialdist(key) + X0)
        y = get_subrange(y, mesh.dy, Y0, -2, +2)
        y = zero_offset(y)
        y = offset_by_axialdist(y, key)
        x = [j * mesh.dy - 0.5 * mesh.Ly for j in range(mesh.Ny)]
        x = get_subrange(x, mesh.dy, Y0, -2, +2)
        plt.plot(x, y, color="black",
                 ls="-")
        
        y = get_nearest_profile(uprime2, mesh, get_axialdist(key) + X0)
        y = get_subrange(y, mesh.dy, Y0, -2, +2)
        y = zero_offset(y)
        y = offset_by_axialdist(y, key)
        x = [j * mesh.dy - 0.5 * mesh.Ly for j in range(mesh.Ny)]
        x = get_subrange(x, mesh.dy, Y0, -2, +2)
        plt.plot(x, y, color="red",
                 ls="--")
        
        x, y = get_xydata(uumean_ref, key)
        y = zero_offset(y)
        y = offset_by_axialdist(y, key)
        plt.plot(x, y,
                 ls="", color="black",
                 marker="d")

    plt.xlim((-2, 2))
    # plt.ylim((-3.5, -0.5))
    plt.xlabel(r"$y$")
    plt.ylabel(r"$\langle u' u' \rangle$")
    plt.savefig("uumean.eps", bbox_inches="tight")

def zero_offset(y):

    offset = 0.5 * (y[0] + y[-1])
    
    zeroed = []
    for val in y:
        zeroed.append(val - offset)
        
    return zeroed

def offset_by_axialdist(y, key):

    offset = get_axialdist(key)

    offset_dat = []
    for val in y:
        offset_dat.append(val - offset)

    return offset_dat

def get_axialdist(key):

    return float(key.replace("x", ""))

def get_xydata(datadict, key):

    x = []
    y = []
    n = len(datadict[key][0])
    for i in range(n):
        x.append(datadict[key][0][i])
        y.append(datadict[key][1][i])

    return x, y

def get_subrange(data, dx, x0, xmin, xmax):

    subrange = []

    for i in range(len(data)):
        x = i * dx - x0
        if (x >= xmin) and (x <= xmax):
            subrange.append(data[i])

    return subrange

def read_refdat(datfile):

    refdat = {}
    with open(datfile, "r") as data:

        header1 = data.readline().split(",")
        keylist = []
        for key in header1:
            if key and not (key=="\n"):
                refdat[key] = [[], []]
                keylist.append(key)

        next(data) # Next line is redundant

        for row in data:
            vals = row.split(",")

            ctr = 0
            for key in keylist:
                refdat[key][0].append(float(vals[2 * ctr + 0]))
                refdat[key][1].append(float(vals[2 * ctr + 1]))
                ctr += 1

    return refdat

def get_nearest_profile(data, mesh, axialdist):

    i = axialdist / mesh.dx
    if i != int(i):
        if i > (int(i) + 0.5):
            i = int(i) + 1
        else:
            i = int(i)
    else:
        i = int(i)

    profile = []
    for j in range(mesh.Ny):
        profile.append(data[i][j])

    return profile

if __name__ == "__main__":
    main()
