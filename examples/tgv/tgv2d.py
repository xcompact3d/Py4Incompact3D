"""
       FILE: tgv2d.py
     AUTHOR: Paul Bartholomew
DESCRIPTION:
"""

import math

import numpy as np
import matplotlib.pyplot as plt

from Py4Incompact3D.postprocess.postprocess import Postprocess
from Py4Incompact3D.tools.misc import avg_over_axis

CASES=["CS6-33",
       "CS6-65",
       "CS6-129",
       "CS6-257",
       "FD2-33",
       "FD2-65",
       "FD2-129",
       "FD2-257"
]
INPUT={"CS6-33":"input-cs633.json",
       "CS6-65":"input-cs665.json",
       "CS6-129":"input-cs6129.json",
       "CS6-257":"input-cs6257.json",
       "FD2-33":"input-fd233.json",
       "FD2-65":"input-fd265.json",
       "FD2-129":"input-fd2129.json",
       "FD2-257":"input-fd2257.json"
}
RE=1600.0
HDR="=" * 72
LINE="-" * 72

NSTEP=10 + 1

def main():

    print(HDR)
    print("Post-processing TGV2D.")
    print(LINE)

    err = {}
    for case in CASES:
        msg = "- Case: " + case
        print(msg)
        
        case_description = case.split("-")
        scheme_key = case_description[0]
        resolution = int(case_description[1])

        if not (scheme_key in err.keys()):
            err[scheme_key] = {}

        postprocess = Postprocess(INPUT[case])

        time = range(NSTEP)
        u = {}
        v = {}
        p = {}
        for t in time:
            msg = "-- timestep " + str(t)
            print(msg)
            
            postprocess.load(time=t)

            u[t] = postprocess.fields["ux"].data[t]
            v[t] = postprocess.fields["uy"].data[t]
            p[t] = postprocess.fields["pp"].data[t]

            postprocess.clear_data()

        err[scheme_key][resolution] = calc_err(u, v, p, time, postprocess.mesh)

    print(LINE)

    print("- Plotting...")

    plt.figure(figsize=(5.0, 3.5))
    for key in err.keys():
        e = []
        n = []
        for r in err[key].keys():
            e.append(err[key][r])
            n.append(r)

        plt.plot(n, e, label=key,
                 ls="",
                 marker="o")

    plt.plot([32, 256], [err["CS6"][33], err["CS6"][33] / (8**6)],
             color="black")
    plt.plot([32, 256], [err["FD2"][33], err["FD2"][33] / (8**2)],
             color="black")
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("N")
    plt.ylabel(r"$\varepsilon_{RMS}$")
    plt.show()

def calc_err(u, v, p, time, mesh):

    err = 0.0

    N = mesh.Nx * mesh.Ny
    X, Y, Z = mesh.get_grid()

    for t in time:
        print(t)
        
        u[t] = avg_over_axis(mesh, u[t], 2)
        v[t] = avg_over_axis(mesh, v[t], 2)
        p[t] = avg_over_axis(mesh, p[t], 2)

        usol, vsol, psol = calc_tgv2d(mesh, t)


        erru = math.sqrt(((u[t] - usol)**2).sum() / N)
        errv = math.sqrt(((v[t] - vsol)**2).sum() / N)
        errp = math.sqrt(((p[t] - psol)**2).sum() / N)
        
        # err += (erru + errv + errp) /  3
        err += (erru + errv) /  2

    err /= len(time)

    return err

def calc_tgv2d(mesh, t):

    u = np.zeros((mesh.Nx, mesh.Ny))
    v = np.zeros((mesh.Nx, mesh.Ny))
    p = np.zeros((mesh.Nx, mesh.Ny))

    F = math.exp(-2 * t / RE)

    for i in range(mesh.Nx):
        x = i * mesh.dx
        
        for j in range(mesh.Ny):
            y = j * mesh.dy

            u[i][j] = math.sin(x) * math.cos(y) * F
            v[i][j] = -math.cos(x) * math.sin(y) * F
            p[i][j] = -(math.cos(2 * x) + math.cos(2 * y)) * F**2 / 4

    return u, v, p

if __name__ == "__main__":
    main()
