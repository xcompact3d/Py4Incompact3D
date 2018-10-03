"""
.. module:: deriv
    :synopsis: Computes the derivatives of data fields using compact finite differences.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

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

def compute_rhs(mesh, field):
	  """ Compute the rhs for the derivative.
	  
	  :param mesh: The mesh on which derivatives are taken.
	  :param field: The field for the variable who's derivative we want.
    
	  :type mesh: Py4Incompact3D.postprocess.mesh.Mesh
    
	  :returns: rhs -- the right-hand side vector.
	  :rtype: numpy.ndarray
	  """
    
	  pass

def deriv(postproc, mesh, phi, axis):
	  """ Take the derivative of field 'phi' along axis. 
    
	  :param postproc: The basic Postprocess object.
	  :param mesh: The mesh on which derivatives are taken.
	  :param phi: The name of the variable who's derivative we want.
	  :param axis: A number indicating direction in which to take derivative: 1=x; 2=y; 3=z.
    
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
