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

from py4incompact3d.postprocess.postprocess import Postprocess
from py4incompact3d.deriv.deriv import deriv
from py4incompact3d.tools.misc import avg_over_axis

INPUT_FILES={"CS6":"input.json",
             # "FD2":"input-2ndorder.json",
             # "FD2-2x":"input-2ndorder-2x.json"
}
RE=4200.0 # Bulk Reynolds number
NTIME={"CS6":100000}
T={"CS6":"0100000"}
HDR="=" * 72
LINE="-" * 72
LINES={"CS6":"-"}

# Location of reference data
REFLEE="./data_lee_retau180.txt"
REFLEEPRIME="./data_lee_fluct_retau180.txt"
REFVREMU="./data_vreman_u_retau180.txt"
REFVREMV="./data_vreman_v_retau180.txt"
REFVREMW="./data_vreman_w_retau180.txt"

# Start and end of arrays of interest
FIRST=0       # Where does the IBM start (N.B. counts from zero)
LAST=FIRST+129 # Where does the channel end?

def load_dataset(input_file, t):

    postprocess = Postprocess(input_file)
    mesh = postprocess.mesh
    yp = mesh.get_grid()[1]

    postprocess.load(time=[t])

    return postprocess.fields["umean"].data[t], \
        postprocess.fields["vmean"].data[t], \
        postprocess.fields["wmean"].data[t], \
        postprocess.fields["uumean"].data[t], \
        postprocess.fields["vvmean"].data[t], \
        postprocess.fields["wwmean"].data[t], \
        mesh

def calc_mean(phi, ntime):
    return phi #/ float(ntime)

def calc_uprime(umean, uumean):
    return (uumean - umean**2)**0.5

def main ():

    print(HDR)
    print("Post-processing channel flow.")
    print(LINE)

    # Load data
    umean = {}; vmean = {}; wmean = {}
    uprime = {}; vprime = {}; wprime = {}
    yp={}
    for input_file in INPUT_FILES.keys():
        t = T[input_file]
        u, v, w, uu, vv, ww, mesh = load_dataset(INPUT_FILES[input_file], t)
        
        # Convert to time mean[input_file]
        u = calc_mean(u, NTIME[input_file])
        v = calc_mean(v, NTIME[input_file])
        w = calc_mean(w, NTIME[input_file])
        uu = calc_mean(uu, NTIME[input_file])
        vv = calc_mean(vv, NTIME[input_file])
        ww = calc_mean(ww, NTIME[input_file])

        # Average in x and z
        u = avg_over_axis(mesh, avg_over_axis(mesh, u, 2), 0)
        v = avg_over_axis(mesh, avg_over_axis(mesh, v, 2), 0)
        w = avg_over_axis(mesh, avg_over_axis(mesh, w, 2), 0)

        uu = avg_over_axis(mesh, avg_over_axis(mesh, uu, 2), 0)
        vv = avg_over_axis(mesh, avg_over_axis(mesh, vv, 2), 0)
        ww = avg_over_axis(mesh, avg_over_axis(mesh, ww, 2), 0)

        # Extract non-IBM subets
        u = u[FIRST:LAST]
        v = v[FIRST:LAST]
        w = w[FIRST:LAST]

        uu = uu[FIRST:LAST]
        vv = vv[FIRST:LAST]
        ww = ww[FIRST:LAST]

        yp[input_file] = mesh.yp[FIRST:LAST]
        yp[input_file] -= yp[input_file][0]

        # Solutions are symmetric
        u = apply_symmetry(u)
        v = apply_symmetry(v)
        w = apply_symmetry(w)

        uu = apply_symmetry(uu)
        vv = apply_symmetry(vv)
        ww = apply_symmetry(ww)

        yp[input_file] = yp[input_file][:len(u)]

        # Compute friction velocity
        # dUdy = (abs(dUdy[0]) + abs(dUdy[-1])) / 2
        dUdy = (u[1] - u[0]) / (yp[input_file][1] - yp[input_file][0])
        tauw = dUdy / RE
        print(tauw)
        utau = math.sqrt(tauw)
        Retau = RE * utau

        msg = " Achieved u_tau = " + str(utau) + "; Re_tau = " + str(Retau)
        print(msg)

        # Re-scale to wall units
        yp[input_file] *= utau * RE
        msg = "y+ = " + str(yp[input_file][1])
        print(msg)
    
        u /= utau
        v /= utau
        w /= utau

        uu /= utau**2
        vv /= utau**2
        ww /= utau**2

        umean[input_file] = u
        vmean[input_file] = v
        wmean[input_file] = w

        # Compute u'
        print(u[0], u[-1])
        uprime[input_file] = calc_uprime(u, uu)
        vprime[input_file] = calc_uprime(v, vv)
        wprime[input_file] = calc_uprime(w, ww)

        # # Limit to yp <= Re_tau
        # umean[input_file] = limit_to_retau(umean[input_file], yp, Retau)
        # uprime = limit_to_retau(uprime, yp, Retau)
        # vprime = limit_to_retau(vprime, yp, Retau)
        # wprime = limit_to_retau(wprime, yp, Retau)
        # yp = limit_to_retau(yp, yp, Retau)

    ## Plot
    print("Plotting...")

    # Umean
    plt.figure(figsize=(5.0, 3.5))
    for input_file in INPUT_FILES.keys():
        plt.plot(yp[input_file], umean[input_file], label=input_file,
                 color="red",
                 ls=LINES[input_file])

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
    plt.savefig("umean" + t + ".png", bbox_inches="tight")
    plt.close()

    # Uprime
    plt.figure(figsize=(5.0, 3.5))
    plt.plot(yp[input_file], uprime[input_file], label="u'",
             color="black",
             ls=LINES[input_file])
    plt.plot(yp[input_file], vprime[input_file], label="v'",
             color="blue",
             ls=LINES[input_file])
    plt.plot(yp[input_file], wprime[input_file], label="w'",
             color="red",
             ls=LINES[input_file])

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
    plt.xlim((yp[input_file][1], 180))
    plt.legend()
    plt.savefig("velprime" + t + ".png", bbox_inches="tight")
    plt.close()

    # # Save to file
    # outfile = "channel.csv"
    # msg = "Writing to file: " + outfile
    # print(msg)

    # with open(outfile, "w") as csvfile:
    #     writer = csv.writer(csvfile, delimiter=" ")
    #     writer.writerow(["yp", "u+", "v+", "w+", "u'", "v'", "w'"])
    #     for j in range(len(yp)):
    #         writer.writerow([yp[j],
    #                          umean[j], vmean[j], wmean[j],
    #                          uprime[j], vprime[j], wprime[j]])

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
        symprof[i] = 0.5 * (prof[0][i][0] + prof[0][-(i + 1)][0])

    if (N / 2) > n:
        symprof[n] = prof[0][n][0]

    return symprof
    
if __name__ == "__main__":
    main()
