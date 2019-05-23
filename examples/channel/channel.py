#!/usr/bin/python3
"""
       FILE: channel.py
     AUTHOR: Paul Bartholomew
DESCRIPTION: 
"""
import csv
import math

import numpy as np
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.rc("font", size=11)

from Py4Incompact3D.postprocess.postprocess import Postprocess
from Py4Incompact3D.deriv.deriv import deriv
from Py4Incompact3D.tools.misc import avg_over_axis

INPUT_FILE="input.json"
RE=4200.0 # Bulk Reynolds number
NTIME=60000
HDR="=" * 72
LINE="-" * 72

# Location of reference data
REFLEE="/home/paul/DATA/benchmarking/channel-flow/data_lee_retau180.txt"
REFVREM="/home/paul/DATA/benchmarking/channel-flow/data_vreman_retau180.txt"

def main ():

    print(HDR)
    print("Post-processing channel flow.")
    print(LINE)
    
    # Load data
    postprocess = Postprocess(INPUT_FILE)
    mesh = postprocess.mesh
    yp = mesh.get_grid()[1]

    t = "0100000" # This is the timestamp of the latest statistic output
    postprocess.load(time=[t])

    # Convert to time mean
    umean = postprocess.fields["umean"].data[t] / float(NTIME)
    vmean = postprocess.fields["vmean"].data[t] / float(NTIME)
    wmean = postprocess.fields["wmean"].data[t] / float(NTIME)

    uumean = postprocess.fields["uumean"].data[t] / float(NTIME)
    vvmean = postprocess.fields["vvmean"].data[t] / float(NTIME)
    wwmean = postprocess.fields["wwmean"].data[t] / float(NTIME)

    # Get dUdy
    dUdy = deriv(postprocess, "umean", 1, t)

    # Average in x and z
    umean = avg_over_axis(mesh, avg_over_axis(mesh, umean, 2), 0)
    vmean = avg_over_axis(mesh, avg_over_axis(mesh, vmean, 2), 0)
    wmean = avg_over_axis(mesh, avg_over_axis(mesh, wmean, 2), 0)

    uumean = avg_over_axis(mesh, avg_over_axis(mesh, uumean, 2), 0)
    vvmean = avg_over_axis(mesh, avg_over_axis(mesh, vvmean, 2), 0)
    wwmean = avg_over_axis(mesh, avg_over_axis(mesh, wwmean, 2), 0)

    dUdy = avg_over_axis(mesh, avg_over_axis(mesh, dUdy, 2), 0)

    # Compute friction velocity
    dUdy = (abs(dUdy[0]) + abs(dUdy[-1])) / 2
    dUdy = ((umean[1] - umean[0]) / (mesh.yp[1] - mesh.yp[0]) \
            + abs((umean[-2] - umean[-1]) / (mesh.yp[-2] - mesh.yp[-1]))) / 2
    tauw = dUdy / RE
    utau = math.sqrt(tauw)
    Retau = RE * utau

    msg = " Achieved u_tau = " + str(utau) + "; Re_tau = " + str(Retau)
    print(msg)

    # Re-scale to wall units
    yp *= utau * RE
    
    umean /= utau
    vmean /= utau
    wmean /= utau

    uumean /= utau**2
    vvmean /= utau**2
    wwmean /= utau**2

    # Compute u'
    uprime = uumean - umean**2
    vprime = vvmean - vmean**2
    wprime = wwmean - wmean**2

    # # Solutions are symmetric
    # umean = apply_symmetry(umean)
    # uprime = apply_symmetry(uprime)
    # vprime = apply_symmetry(vprime)
    # wprime = apply_symmetry(wprime)

    # yp = yp[:len(umean)]

    # Plot
    print("Plotting...")
    
    plt.plot(yp, umean, label="X3D",
             color="red")

    yplee, ulee = read_lee(REFLEE)
    plt.plot(yplee, ulee, label="LEE",
             color="blue",
             ls="", marker="*")

    ypvrem, uvrem = read_vreman(REFVREM)
    plt.plot(ypvrem, uvrem, label="VREMAN",
             color="black",
             ls="", marker="+")
    
    plt.xlabel(r"$y_+$")
    plt.ylabel(r"$\left< u_+ \right>$")
    plt.xscale("log")
    plt.xlim((yp[1], 200))
    plt.ylim(ymin=0)
    plt.legend(loc="upper left",
               numpoints=1)
    plt.savefig("umean.eps", bbox_inches="tight")

    # Save to file
    outfile = "channel.csv"
    msg = "Writing to file: " + outfile
    print(msg)

    with open(outfile, "w") as csvfile:
        writer = csv.writer(csvfile, delimiter=" ")
        writer.writerow(["yp", "u+", "v+", "w+", "u'", "v'", "w'"])
        for j in range(len(yp)):
            writer.writerow([yp[j],
                             umean[j], vmean[j], wmean[j],
                             uprime[j], vprime[j], wprime[j]])

    print(LINE)

def read_lee(filepath):

    u = []
    yp = []
    with open(filepath, "r") as datfile:
        for row in datfile:
            if "%" not in row:
                words = row.split()
                if len(words):
                    u.append(float(words[2]))
                    yp.append(float(words[1]))

    return yp, u

def read_vreman(filepath):

    u = []
    yp = []
    with open(filepath, "r") as datfile:
        for row in datfile:
            if "%" not in row:
                words = row.split()
                if len(words):
                    u.append(float(words[1]))
                    yp.append(float(words[0]))

    return yp, u
    
if __name__ == "__main__":
    main()
