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
    for i in range(n):

        if (Re_theta[i] > 668) and (Re_theta[i] < 672):
            i_670 = i
            print(f"Selected theta 670 {Re_theta[i]} {i_670} {H_12[i]} {cf[i]} {utau[i]} {Re_tau[i]} {(300 / (mesh.Nx-1)) * utau[i] * Re} {(5 / mesh.Nz) * utau[i] * Re}")
        if (Re_theta[i] > 998) and (Re_theta[i] < 1002):
            i_1000 = i
            print(f"Selected theta 1000 {Re_theta[i]} {i_1000} {H_12[i]} {cf[i]} {utau[i]} {Re_tau[i]} {(300 / (mesh.Nx-1)) * utau[i] * Re} {(5 / mesh.Nz) * utau[i] * Re}")

    u_670 = u[i_670, j0:je]
    up_670 = uprime[i_670, j0:je]
    u_1000 = u[i_670, j0:je]
    up_1000 = uprime[i_670, j0:je]
    
if __name__ == "__main__":
    main()
