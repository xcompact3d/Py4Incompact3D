#!/usr/bin/python3
"""
       FILE: channel.py
     AUTHOR: Paul Bartholomew
DESCRIPTION: 
"""
import csv
import math

import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.rc("font", size=11)

from Py4Incompact3D.postprocess.postprocess import Postprocess
from Py4Incompact3D.deriv.deriv import deriv
from Py4Incompact3D.tools.misc import avg_over_axis

#INPUT_FILE="input.json"
INPUT_FILE="input-2ndorder.json"
RE=4200.0 # Bulk Reynolds number
NTIME=150000 #400000
HDR="=" * 72
LINE="-" * 72

# Location of reference data
REFLEE="/home/paul/DATA/benchmarking/channel-flow/data_lee_retau180.txt"
REFLEEPRIME="/home/paul/DATA/benchmarking/channel-flow/data_lee_fluct_retau180.txt"
REFVREMU="/home/paul/DATA/benchmarking/channel-flow/data_vreman_u_retau180.txt"
REFVREMV="/home/paul/DATA/benchmarking/channel-flow/data_vreman_v_retau180.txt"
REFVREMW="/home/paul/DATA/benchmarking/channel-flow/data_vreman_w_retau180.txt"

# Start and end of arrays of interest
FIRST=0       # Where does the IBM start (N.B. counts from zero)
LAST=FIRST+129 # Where does the channel end?

def main ():

    print(HDR)
    print("Post-processing channel flow.")
    print(LINE)
    
    # Load data
    postprocess = Postprocess(INPUT_FILE)
    mesh = postprocess.mesh
    yp = mesh.get_grid()[1]

    # t = "0170000" # This is the timestamp of the latest statistic output
    t = "0600000" # This is the timestamp of the latest statistic output
    # t = "_full" # This is the "timestamp" of the full statistics
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

    # dUdy = avg_over_axis(mesh, avg_over_axis(mesh, dUdy, 2), 0)

    # Extract non-IBM subets
    umean = umean[FIRST:LAST]
    vmean = vmean[FIRST:LAST]
    wmean = wmean[FIRST:LAST]

    uumean = uumean[FIRST:LAST]
    vvmean = vvmean[FIRST:LAST]
    wwmean = wwmean[FIRST:LAST]

    # dUdy = dUdy[FIRST:LAST]
    yp = mesh.yp[FIRST:LAST]
    yp -= yp[0]

    # Solutions are symmetric
    umean = apply_symmetry(umean)
    vmean = apply_symmetry(vmean)
    wmean = apply_symmetry(wmean)

    uumean = apply_symmetry(uumean)
    vvmean = apply_symmetry(vvmean)
    wwmean = apply_symmetry(wwmean)

    print(umean[-1])

    yp = yp[:len(umean)]

    # Compute friction velocity
    # dUdy = (abs(dUdy[0]) + abs(dUdy[-1])) / 2
    dUdy = (umean[1] - umean[0]) / (yp[1] - yp[0])
    tauw = dUdy / RE
    utau = math.sqrt(tauw)
    yp1 = utau * yp[1] / RE
    Retau = RE * utau

    msg = " Achieved u_tau = " + str(utau) + "; Re_tau = " + str(Retau)
    print(msg)

    # Re-scale to wall units
    yp *= utau * RE
    msg = "y+ = " + str(yp[1])
    print(msg)
    
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
    for i in range(len(uprime)):
        uprime[i] = math.sqrt(uprime[i])
        vprime[i] = math.sqrt(vprime[i])
        wprime[i] = math.sqrt(wprime[i])

    # # Limit to yp <= Re_tau
    # umean = limit_to_retau(umean, yp, Retau)
    # uprime = limit_to_retau(uprime, yp, Retau)
    # vprime = limit_to_retau(vprime, yp, Retau)
    # wprime = limit_to_retau(wprime, yp, Retau)
    # yp = limit_to_retau(yp, yp, Retau)

    # Plot
    print("Plotting...")

    plt.figure(figsize=(5.0, 3.5))
    plt.plot(yp, umean, label="X3D",
             color="red")

    yplee, ulee = read_lee(REFLEE)
    plt.plot(yplee, ulee, label="LEE",
             color="black",
             ls="", marker="x",
             markevery=2)

    ypvrem, uvrem = read_vreman(REFVREMU)
    plt.plot(ypvrem, uvrem, label="VREMAN",
             color="black",
             ls="", marker="+",
             markevery=3)
    
    plt.xlabel(r"$y_+$")
    plt.ylabel(r"$U_+$")
    plt.xscale("log")
    # plt.xlim((yp[1], 200))
    plt.xlim((1, 200))
    plt.ylim(ymin=0)
    plt.legend(loc="upper left",
               numpoints=1)
    plt.savefig("umean" + t + ".eps", bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(5.0, 3.5))
    plt.plot(yp, uprime, label="u'",
             color="black")
    plt.plot(yp, vprime, label="v'",
             color="blue")
    plt.plot(yp, wprime, label="w'",
             color="red")

    yplee, uprime_lee, vprime_lee, wprime_lee = read_lee_prime(REFLEEPRIME)
    plt.plot(yplee, uprime_lee,
             color="black",
             ls="", marker="x",
             markevery=2)
    plt.plot(yplee, vprime_lee,
             color="blue",
             ls="", marker="x",
             markevery=2)
    plt.plot(yplee, wprime_lee,
             color="red",
             ls="", marker="x",
             markevery=2)
    ypvrem, uprime_vrem, vprime_vrem, wprime_vrem = read_vreman_prime(REFVREMU, REFVREMV, REFVREMW)
    plt.plot(ypvrem, uprime_vrem,
             color="black",
             ls="", marker="+",
             markevery=3)
    plt.plot(ypvrem, vprime_vrem,
             color="blue",
             ls="", marker="+",
             markevery=3)
    plt.plot(ypvrem, wprime_vrem,
             color="red",
             ls="", marker="+",
             markevery=3)

    plt.xlabel(r"$y_+$")
    plt.ylabel(r"$\langle u'_+ \rangle$")
    plt.xlim((yp[1], 180))
    plt.legend()
    plt.savefig("velprime" + t + ".eps", bbox_inches="tight")
    plt.close()

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

def read_lee_prime(filepath):

    u = []
    v = []
    w = []
    yp = []
    with open(filepath, "r") as datfile:
        for row in datfile:
            if "%" not in row:
                words = row.split()
                if len(words):
                    u.append(math.sqrt(float(words[2])))
                    v.append(math.sqrt(float(words[3])))
                    w.append(math.sqrt(float(words[4])))
                    yp.append(float(words[1]))

    return yp, u, v, w

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

def read_vreman_prime(filepathu, filepathv, filepathw):

    u = []
    v = []
    w = []
    yp = []
    with open(filepathu, "r") as datfile:
        for row in datfile:
            if "%" not in row:
                words = row.split()
                if len(words):
                    u.append(float(words[2]))
                    yp.append(float(words[0]))
    with open(filepathv, "r") as datfile:
        for row in datfile:
            if "%" not in row:
                words = row.split()
                if len(words):
                    v.append(float(words[2]))
    with open(filepathw, "r") as datfile:
        for row in datfile:
            if "%" not in row:
                words = row.split()
                if len(words):
                    w.append(float(words[2]))

    return yp, u, v, w

def limit_to_retau(phi, yp, retau):

    phi_limited = []
    for i in range(len(yp)):
        if yp[i] <= retau:
            phi_limited.append(phi[i])
    return phi_limited

def apply_symmetry(prof):

    N = len(prof)
    n = int(N / 2)

    if (N / 2) > n:
        symprof = np.zeros(n + 1)
    else:
        symprof = np.zeros(n)

    for i in range(n):
        symprof[i] = 0.5 * (prof[i] + prof[-(i + 1)])

    if (N / 2) > n:
        symprof[n] = prof[n]

    return symprof
    
if __name__ == "__main__":
    main()
