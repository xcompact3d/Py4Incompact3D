"""
.. module:: deriv
    :synopsis: Computes the derivatives of data fields using compact finite differences.

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

def tdma(a, b, c, rhs):
	""" The Tri-Diagonal Matrix Algorithm. 
	
	Solves tri-diagonal matrices using TDMA.

	:param a:
	:type a: numpy.ndarray 
	:param b:
	:param c:
	:param rhs:
	:returns: dphidx -- the derivative
	:rtype: numpy.ndarray
	"""
	pass

def compute_rhs():
	return rhs

def compute_derivvars():
	return derivvars

def deriv(postproc, mesh, phi, axis):
	""" Take the derivative of field 'phi' along axis. 

	:param postproc:
	:param mesh:
	:param phi: The name of the variable who's derivative we want.
	:type phi: str
	:param axis:
	"""

	derivvars = compute_derivvars()

	# Transpose the data to make loops more efficient
	rhs = compute_rhs()
	return tdma(rhs)
