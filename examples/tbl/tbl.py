#!/usr/bin/python3
"""
       FILE: tbl.py
     AUTHOR: Paul Bartholomew
DESCRIPTION: 
"""

import numpy as np

def calc_uprime(u, uu):
    up = uu - u**2
    up[up < 0] = 0 # Remove -ve values
    return up**0.5

def friction_vel(u, Re, dy, ppy):

    u_tau = -(11.0 / 5.0) * u[:, 0] + 3 * u[:, 1] - 1.5 * u[:, 2] + (1.0 / 3.0) * u[:, 3]
    u_tau /= dy
    u_tau *= ppy[0]

    return (u_tau / Re)**0.5

def friction_factor(u_tau):
    return 2 * u_tau**2

def friction_reynolds(Re, u_tau, delta_99):
    return Re * u_tau * delta_99

def momentum_reynolds(Re, theta):
    return Re * theta

def shape_factor(delta_s, theta):
    return delta_s / theta

def wake_factor(u_tau, delta_99, Re):
    return 0.5 * 0.41 * (0.99 / u_tau - 5.2 - np.log(delta_99 * u_tau * Re) / 0.41)

def displacement_thickness(u, us):
    return 1 - u / us

def momentum_thickness(u, us):
    return (1 - u / us) * (u / us)

def integrate(u, y, j0, je):

    u0 = u[:, j0:je-1]
    u1 = u[:, j0+1:je]
    y0 = y[j0:je-1]
    y1 = y[j0+1:je]
    
    return (0.5 * (u1 + u0) * (y1 - y0)).sum(axis=1)

def boundary99_thickness(u, y, us, j0, je):

    u0 = u[:, j0:je-1]
    u1 = u[:, j0+1:je]
    y0 = y[j0:je-1]
    y1 = y[j0+1:je]

    den = u1 - u0
    den[den == 0] = 1
    delta_99 = ((u1 - 0.99 * us) * y0 + (0.99 * us - u0) * y1) / den
    delta_99[den == 0] = 0

    # Apply masking operation
    mask = u1 / us
    mask[mask < 0.99] = 0
    mask[mask != 0] = 1
    delta_99 = mask * delta_99
    
    # Zero out values below threshold
    threshold = 0.0001
    delta_99[delta_99 < threshold] = 0 

    # Take first non-zero value at every i
    non_zero_indices = (delta_99 != 0).argmax(axis=1)
    nrows = delta_99.shape[0]
    return delta_99[[i for i in range(nrows)], non_zero_indices]
    
# Read reference data
def refdat(reffile):

    # Read y+ (col 1), U+ (col 2), u'RMS (col 3)
    # Reference data file has been edited to mark header/comment lines beggining with '#'
    yp = []; u = []; uprime = []
    with open(reffile, "r") as data:
        for row in data:
            if row[0] != "#":
                words = row.split()
                yp.append(float(words[1]))
                u.append(float(words[2]))
                uprime.append(float(words[3]))

    return yp, u, uprime

def main():

    from py4incompact3d.postprocess.postprocess import Postprocess
    from py4incompact3d.tools.misc import avg_over_axis

    Re = 2000.0
    
    # Load data
    t = "600000"
    postprocess = Postprocess("input.json")
    postprocess.load(time=[t])

    u = postprocess.fields["umean"].data[t] 
    uu = postprocess.fields["uumean"].data[t]

    uprime = calc_uprime(u, uu)

    mesh = postprocess.mesh
    
    # Average in Z
    u = avg_over_axis(mesh, u, 2)[:,:,0]
    uprime = avg_over_axis(mesh, uprime, 2)[:,:,0]

    # Get freestream
    us = u[:, mesh.Ny - 5]

    # Compute flow quantities
    dya = mesh.Ly / (mesh.Ny - 1.0)
    utau = friction_vel(u, Re, dya, mesh.ppy)
    cf = friction_factor(utau)

    j0 = 0
    je = mesh.Ny - 10
    delta_s = integrate(displacement_thickness(u, us[:, None]), mesh.yp, j0, je)
    theta = integrate(momentum_thickness(u, us[:, None]), mesh.yp, j0, je)
    delta_99 = boundary99_thickness(u, mesh.yp, us[:, None], j0, je)

    Re_tau = friction_reynolds(Re, utau, delta_99)
    Re_theta = momentum_reynolds(Re, theta)
    H_12 = shape_factor(delta_s, theta)
    Wk = wake_factor(utau, delta_99, Re)

    n = Re_theta.shape[0]
    i_670 = 0; i_1000 = 0

    def print_stats(Re_tgt, Re_theta, i, H12, cf, utau, Re_tau, Re, nx, nz):
        print(f"Selected theta {Re_tgt} {Re_theta} {i} {H12} {cf} {utau} {Re_tau} {(300 / (nx - 1)) * utau * Re} {(5 / nz) * utau * Re}")
        
    for i in range(n):

        if (Re_theta[i] > 668) and (Re_theta[i] < 672):
            i_670 = i
            print_stats(670, Re_theta[i], i, H_12[i], cf[i], utau[i], Re_tau[i], Re, mesh.Nx, mesh.Nz)
        if (Re_theta[i] > 998) and (Re_theta[i] < 1002):
            i_1000 = i
            print_stats(1000, Re_theta[i], i, H_12[i], cf[i], utau[i], Re_tau[i], Re, mesh.Nx, mesh.Nz)

    u_670 = u[i_670, :]
    up_670 = uprime[i_670, :]
    u_1000 = u[i_1000, :]
    up_1000 = uprime[i_1000, :]

    # First find location of this file, the reference data should be stored alongside.
    import os
    mypath = os.path.dirname(os.path.abspath(__file__))
    reffile = os.path.join(mypath, "vel_1000.prof.txt")

    # Plotting
    def plot(u, uprime, utau, Re_theta, reffile=None):
        import matplotlib.pyplot as plt

        def plotone(u, xlabel, varname, yp_ref=None, u_ref=None):
            plt.plot(mesh.yp * utau * Re, u / utau,
                     ls="",
                     marker="+",
                     label="X3D")
            if ((yp_ref != None) and (u_ref != None)):
                plt.plot(yp_ref, u_ref,
                         label="Ref")
            plt.xscale("log")
            plt.xlabel(xlabel)
            plt.ylabel("Y+")
            plt.legend()
            figname = varname + "_" + str(Re_theta) + ".pdf"
            plt.savefig(figname)
            plt.close()

        if (reffile != None):
            yp_ref, u_ref, uprime_ref = refdat(reffile)
            plotone(u, "U+", "u", yp_ref, u_ref)
            plotone(uprime, "U'", "uprime", yp_ref, uprime_ref)
        else:
            plotone(u, "U+", "u")
            plotone(uprime, "U'", "uprime")

    plot(u_670, up_670, utau[i_670], 670)
    plot(u_1000, up_1000, utau[i_1000], 1000, reffile)

    def data_to_file(filename, datasets):

        n = len(datasets[0])
        with open(filename, "w") as output:
            for i in range(n):
                line = ""
                for d in datasets:
                    line += str(d[i]) + " "
                line += "\n"
                output.write(line)

    data_to_file("u_prof670.dat",
                 [mesh.yp * utau[i_670] * Re,
                  u_670 / utau[i_670],
                  up_670 / utau[i_670]])
    data_to_file("u_prof1000.dat",
                 [mesh.yp * utau[i_1000] * Re,
                  u_1000 / utau[i_1000],
                  up_1000 / utau[i_1000]])
    data_to_file("evol.dat",
                 [[i for i in range(mesh.Nx)],
                  cf,
                  delta_99,
                  delta_s,
                  theta,
                  Re_theta,
                  H_12,
                  Re_tau])

if __name__ == "__main__":
    main()
