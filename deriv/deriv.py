"""
.. module:: deriv
    :synopsis: Computes the derivatives of data fields using compact finite differences.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

import numpy as np

def tdma(a, b, c, rhs):
    """ The Tri-Diagonal Matrix Algorithm.

    Solves tri-diagonal matrices using TDMA.

    :param a:
    :param b:
    :param c:
    :param rhs:

    :type a: numpy.ndarray

    :returns: dphidx -- the derivative
    :rtype: numpy.ndarray
    """
    pass

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

    rhs = np.zeros([field.shape[0], field.shape[1], field.shape[2]])

    if axis in field.vidx:
        # npaire = 0
        pass
    else:
        # npaire = 1
        pass
    
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

    rhs = compute_rhs(mesh, postproc.fields[phi])
    return tdma(rhs)
