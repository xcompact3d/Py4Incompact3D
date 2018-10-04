"""
.. module:: deriv
    :synopsis: Computes the derivatives of data fields using compact finite differences.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

import numpy as np

def tdma(a, b, c, rhs):
    """ The Tri-Diagonal Matrix Algorithm.

    Solves tri-diagonal matrices using TDMA where the matrices are of the form
    [b0 c0
     a1 b1 c1
        a2 b2   c2

           an-2 bn-2 cn-1
                an-1 bn-1]

    :param a: The `left' coefficients.
    :param b: The diagonal coefficients. (All ones?)
    :param c: The 'right' coefficients.
    :param rhs: The right-hand side vector.

    :type a: numpy.ndarray
    :type b: numpy.ndarray
    :type c: numpy.ndarray
    :type rhs: numpy.ndarray

    :returns: rhs -- the rhs vector overwritten with derivatives.
    :rtype: numpy.ndarray
    """
    
    for i in range(rhs.shape[0]):
        for j in range(rhs.shape[1]):
            # Forward elimination
            for k in range(1, rhs.shape[2]):
                m = a[k] / b[k - 1]
                b[k] -= m * c[k - 1]
                rhs[i][j][k] -= m * c[k - 1]

            # Backward substituion
            rhs[i][j][-1] /= b[-1]
            for k in range(rhs.shape[2] - 2, -1, -1):
                rhs[i][j][k] -= c[k] * rhs[k + 1]
                rhs[i][k][k] /= b[k]

    return rhs

def compute_rhs_0(mesh, field, axis):
    """ Compute the rhs for the derivative for periodic BCs.

    :param mesh: The mesh on which derivatives are taken.
    :param field: The field for the variable who's derivative we want.
    :param axis: A number indicating direction in which to take derivative: 0=x; 1=y; 2=z.

    :type mesh: Py4Incompact3D.postprocess.mesh.Mesh
    :type axis: int

    :returns: rhs -- the right-hand side vector.
    :rtype: numpy.ndarray
    """

    # Setup
    rhs = np.zeros([field.shape[0], field.shape[1], field.shape[2]])

    if axis == 0:
        invdx = 1.0 / mesh.dx
    elif axis == 1:
        invdx = 1.0 / mesh.dy
    else:
        invdx = 1.0 / mesh.dz

    a = mesh.a * invdx / 2.0
    b = mesh.b * invdx / 4.0

    # Compute RHS
    for i in range(field.shape[0]):
        for j in range(field.shape[1]):
            # XXX Due to python's negative indices, BC @ k = 0 automatically applied
            for k in range(field.shape[2] - 2):
                rhs[i][j][k] = a * (field[i][j][k + 1] - field[i][j][k - 1]) \
                               + b * (field[i][j][k + 2] - field[i][j][k - 2])
                
            # BCs @ k = n
            k = field.shape[2] - 2
            rhs[i][j][k] = a * (field[i][j][k + 1] - field[i][j][k - 1]) \
                           + b * (field[i][j][0] - field[i][j][k - 2])
            k = field.shape[2] - 1
            rhs[i][j][k] = a * (field[i][j][0] - field[i][j][k - 1]) \
                           + b * (field[i][j][1] - field[i][j][k - 2])

    return rhs

def compute_rhs_1(mesh, field, axis):
    """ Compute the rhs for the derivative for free slip BCs.

    :param mesh: The mesh on which derivatives are taken.
    :param field: The field for the variable who's derivative we want.
    :param axis: A number indicating direction in which to take derivative: 0=x; 1=y; 2=z.

    :type mesh: Py4Incompact3D.postprocess.mesh.Mesh
    :type axis: int

    :returns: rhs -- the right-hand side vector.
    :rtype: numpy.ndarray
    """

    # Setup
    rhs = np.zeros([field.shape[0], field.shape[1], field.shape[2]])

    if axis == 0:
        invdx = 1.0 / mesh.dx
    elif axis == 1:
        invdx = 1.0 / mesh.dy
    else:
        invdx = 1.0 / mesh.dz

    a = mesh.a * invdx / 2.0
    b = mesh.b * invdx / 4.0

    for i in range(field.shape[0]):
        for j in range(field.shape[1]):
            #BCs @ k = 0
            if axis in field.vidx:
                # npaire = 0
                k = 0
                rhs[i][j][k] = 0.0
                k = 1
                rhs[i][j][k] = a * (field[i][j][k + 1] - field[i][j][k - 1]) \
                               + b * (field[i][j][k + 2] - field[i][j][k])
            else:
                #npaire = 1
                k = 0
                rhs[i][j][k] = 2 * (a * field[i][j][k + 1] + b * field[i][j][k + 2])
                k = 1
                rhs[i][j][k] = a * (field[i][j][k + 1] - field[i][j][k - 1]) \
                               + b * (field[i][j][k + 2] + field[i][j][k])

            # Internal nodes
            for k in range(2, field.shape[2] - 2):
                rhs[i][j][k] = a * (field[i][j][k + 1] - field[i][j][k - 1]) \
                               + b * (field[i][j][k + 2] - field[i][j][k - 2])

            # BCs @ k = n
            k = field.shape[2] - 2
            rhs[i][j][k] = a * (field[i][j][k + 1] - field[i][j][k - 1]) \
                           + b * (field[i][j][k] - field[i][j][k - 2])
            k = field.shape[2] - 1
            if axis in field.vidx:
                # npaire = 0
                rhs[i][j][k] = 0.0
            else:
                # npaire = 1
                rhs[i][j][k] = 2 * (a * field[i][j][k - 1] + b * field[i][j][k - 2])

    return rhs

def compute_rhs_2(mesh, field, axis):
    """ Compute the rhs for the derivative for Dirichlet BCs.

    :param mesh: The mesh on which derivatives are taken.
    :param field: The field for the variable who's derivative we want.
    :param axis: A number indicating direction in which to take derivative: 0=x; 1=y; 2=z.

    :type mesh: Py4Incompact3D.postprocess.mesh.Mesh
    :type axis: int

    :returns: rhs -- the right-hand side vector.
    :rtype: numpy.ndarray
    """

    # Setup
    rhs = np.zeros([field.shape[0], field.shape[1], field.shape[2]])

    if axis == 0:
        invdx = 1.0 / mesh.dx
    elif axis == 1:
        invdx = 1.0 / mesh.dy
    else:
        invdx = 1.0 / mesh.dz

    a = mesh.a * invdx / 2.0
    b = mesh.b * invdx / 4.0

    for i in range(field.shape[0]):
        for j in range(field.shape[1]):
            # BCs @ k = 0
            k = 0
            rhs[i][j][k] = -(5.0 * field[i][j][k] - 4.0 * field[i][j][k + 1] - field[i][j][k + 2]) \
                           * (0.5 * invdx)
            k = 1
            rhs[i][j][k] = 1.5 * (field[i][j][k + 1] - field[i][j][k - 1]) * (0.5 * invdx)
            
            # Internal nodes
            for k in range(2, field.shape[2] - 2):
                rhs[i][j][k] = a * (field[i][j][k + 1] - field[i][j][k - 1]) \
                               + b * (field[i][j][k + 2] - field[i][j][k - 2])
                
            # BCs @ k = n
            k = field.shape[2] - 2
            rhs[i][j][k] = 1.5 * (field[i][j][k + 1] - field[i][j][k - 1]) * (0.5 * invdx)
            k = field.shape[2] - 1
            rhs[i][j][k] = (field[i][j][k] + 4.0 * field[i][j][k - 1] - 5.0 * field[i][j][k - 2]) \
                           * (0.5 * invdx)

    return rhs

def compute_rhs(mesh, field, axis):
    """ Compute the rhs for the derivative.

    :param mesh: The mesh on which derivatives are taken.
    :param field: The field for the variable who's derivative we want.
    :param axis: A number indicating direction in which to take derivative: 0=x; 1=y; 2=z.

    :type mesh: Py4Incompact3D.postprocess.mesh.Mesh
    :type axis: int

    :returns: rhs -- the right-hand side vector.
    :rtype: numpy.ndarray
    """

    if axis == 0:
        if mesh.BCx == 0:
            return compute_rhs_0(mesh, field, axis)
        elif mesh.BCx == 1:
            return compute_rhs_1(mesh, field, axis)
        else:
            return compute_rhs_2(mesh, field, axis)
    elif axis == 1:
        if mesh.BCy == 0:
            return compute_rhs_0(mesh, field, axis)
        elif mesh.BCy == 1:
            return compute_rhs_1(mesh, field, axis)
        else:
            return compute_rhs_2(mesh, field, axis)
    else:
        if mesh.BCz == 0:
            return compute_rhs_0(mesh, field, axis)
        elif mesh.BCz == 1:
            return compute_rhs_1(mesh, field, axis)
        else:
            return compute_rhs_2(mesh, field, axis)

def deriv(postproc, mesh, phi, axis):
    """ Take the derivative of field 'phi' along axis.

    :param postproc: The basic Postprocess object.
    :param mesh: The mesh on which derivatives are taken.
    :param phi: The name of the variable who's derivative we want.
    :param axis: A number indicating direction in which to take derivative: 0=x; 1=y; 2=z.

    :type postproc: Py4Incompact3D.postprocess.postprocess.Postprocess
    :type mesh: Py4Incompact3D.postprocess.mesh.Mesh
    :type phi: str
    :type axis: int

    :returns: dphidx -- the derivative
    :rtype: numpy.ndarray
    """

    mesh.compute_derivvars()

    # Transpose the data to make loops more efficient
    postproc.fields[phi] = np.swapaxes(postproc.fields[phi], axis, 2)

    rhs = compute_rhs(mesh, postproc.fields[phi])

    # Transpose back to normal orientation and return
    postproc.fields[phi] = np.swapaxes(postproc.fields[phi], 2, axis)
    rhs = np.swapaxes(rhs, 2, axis)
    return tdma(rhs)
