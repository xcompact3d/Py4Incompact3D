#!/usr/bin/python3
"""
       FILE: channel.py
     AUTHOR: Paul Bartholomew
DESCRIPTION: 
"""
import csv
import math

import matplotlib.pyplot as plt

from Py4Incompact3D.postprocess.postprocess import Postprocess
from Py4Incompact3D.deriv.deriv import deriv

INPUT_FILE="input.json"
RE=4200.0 # Bulk Reynolds number
NTIME=60000
HDR="=" * 72
LINE="-" * 72

def main ():

    print(HDR)
    print("Post-processing channel flow.")
    print(LINE)
    
    # Load data
    postprocess = Postprocess(INPUT_FILE)
    mesh = postprocess.mesh
    yp = mesh.get_grid()[1]

    t = 0
    postprocess.load(time=t)

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
    dUdy = (dUdy[0] + dUdy[-1]) / 2
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

    # Plot
    print("Plotting...")
    
    plt.plot(umean, yp)
    plt.xlabel(r"$\left< u_+ \right>$")
    plt.ylabel(r"$y_+$")
    plt.show()
    
    plt.plot(uprime, yp)
    plt.xlabel(r"$\left< u'_+$ \right>")
    plt.ylabel(r"$y_+$")
    plt.show()
    
    plt.plot(vprime, yp)
    plt.xlabel(r"$\left< v'_+$ \right>")
    plt.ylabel(r"$y_+$")
    plt.show()
    
    plt.plot(wprime, yp)
    plt.xlabel(r"$\left< w'_+$ \right>")
    plt.ylabel(r"$y_+$")
    plt.show()

    # Save to file
    outfile = "channel.csv"
    msg = "Writing to file: " + outfile
    print(msg)

    with open(outfile, "w") as csvfile:
        writer = csv.writer(csvfile, delimeter=" ")
        writer.writerow(["yp", "u+", "v+", "w+", "u'", "v'", "w'"])
        for j in range(len(yp)):
            writer.writerow(yp[j],
                            umean[j], vmean[j], wmean[j],
                            uprime[j], vprime[j], wprime[j])

    print(LINE)
    
if __name__ == "__main__":
    main()
