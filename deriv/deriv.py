"""
       FILE: deriv.py
     AUTHOR: Paul Bartholomew
DESCRIPTION: Compact derivatives - this is essentially a translation from incompact to python,
             allowing for post-processing data requiring derivatives.
"""

def derx_0(phi, nx, ny, nz, dx):
	""" X-derivative for periodic BCs. """
	
	tx = [[[0 for k in range(nz)] for j in range(ny)] for i in range(nx)]
	rx = [[[0 for k in range(nz)] for j in range(ny)] for i in range(nx)]

	# Constants
	afix = (7.0 / 9.0) / dx
	bfix = (1.0 / 36.0) / dx
	alfaix = 1.0 / 3.0

	# BC @ x = 0
	for j in range(ny):
		for k in range(nz):
			tx[0][j][k] = afix * (phi[1][j][k] - phi[-1][j][k]) + bfix * (phi[2][j][k] - phi[-2][j][k])
			rx[0][j][k] = -1.0
			tx[1][j][k] = afix * (phi[2][j][k] - phi[0][j][k]) + bfix * (phi[3][j][k] - phi[-1][j][k])
			rx[1][j][k] = 0.0

	# Loop over internal nodes
	for i in range(2, nx - 2):
		for j in range(ny):
			for k in range(nz):
				tx[i][j][k] = afix * (phi[i + 1][j][k] - phi[i - 1][j][k]) \
											+ bfix * (phi[i + 2][j][k] - phi[i - 2][j][k])
				rx[i][j][k] = 0.0

	# BC @ x = nx
	for j in range(ny):
		for k in range(nz):
			tx[-2][j][k] = afix * (phi[-1][j][k] - phi[-3][j][k]) + bfix * (phi[0][j][k] - phi[-4][j][k])
			rx[-2][j][k] = 0.0
			tx[-1][j][k] = afix * (phi[0][j][k] - phi[-2][j][k]) + bfix * (phi[1][j][k] - phi[-3][j][k])
			rx[-1][j][k] = alfaix

	# Tridiagonal solve
	for i in range(1, nx):
		for j in range(ny):
			for k in range(nz):
				tx[i][j][k] -= tx[i - 1][j][k] * fsx[i]
				rx[i][j][k] -= rx[i - 1][j][k] * fsx[i]
	for j in range(ny):
		for k in range(nz):
			tx[-1][j][k] *= fwx[-1]
			rx[-1][j][k] *= fwx[-1]
	for i in range(nx - 2, -1, -1):
		for j in range(ny):
			for k in range(nz):
				tx[i][j][k] -= tx[i + 1][j][k] * ffx[i]
				tx[i][j][k] *= fwx[i]
				rx[i][j][k] -= rx[i + 1][j][k] * ffx[i]
				rx[i][j][k] *= fwx[i]
	for j in range(ny):
		for k in range(nz):
			sx[j][k] = (tx[1][j][k] - alfaix * tx[-1][j][k]) / (1.0 + rx[1][j][k] - alfaix * rx[-1][j][k])
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				tx[i][j][k] -= sx[j][k] * rx[i][j][k]

	return tx

def derx_1(phi, nx, ny, nz, npaire):
	""" X-derivative for Neumann BCs. """

	tx = [[[0 for k in range(nz)] for j in range(ny)] for i in range(nx)]
	
	if npaire == 1:

		# BCs @ x = 0
		for j in range(ny):
			for k in range(nz):
				tx[0][j][k] = 0.0
				tx[1][j][k] = afix * (phi[2][j][k] - phi[0][j][k]) + bfix * (phi[3][j][k] - phi[1][j][k])

		# Internal nodes
		for i in range(2, nx - 2):
			for j in range(ny):
				for k in range(nz):
					tx[i][j][k] = afix * (phi[i + 1][j][k] - phi[i - 1][j][k]) \
												+ bfix * (phi[i + 2][j][k] - phi[i - 2][j][k])

		# BCs @ x = nx
		for j in range(ny):
			for k in range(nz):
				tx[-2][j][k] = afix * (phi[-1][j][k] - phi[-3][j][k]) \
											 + bfix * (phi[-2][j][k] - phi[-4][j][k])
				tx[-1][j][k] = 0.0

		# Tri-diagonal solve
		for i in range(1, nx):
			for j in range(ny):
				for k in range(nz):
					tx[i][j][k] -= fsx[i] * tx[i - 1][j][k]
		for j in range(ny):
			for k in range(nz):
				tx[-1][j][k] *= fwx[-1]
		for i in range(nx - 2, -1, -1):
			tx[i][j][k] -= ffx[i] * tx[i + 1][j][k]
			tx[i][j][k] *= fwx[i]
	else:

		# BCs @ x = 0
		for j in range(ny):
			for k in range(nz):
				tx[0][j][k] = 2 * (afix * phi[1][j][k] + bfix * phi[2][j][k])
				tx[1][j][k] = afix * (phi[2][j][k] - ux[0][j][k]) + bfix * (phi[3][j][k] + phi[1][j][k])

		# Internal nodes
		for i in range(2, nx - 2):
			for j in range(ny):
				for k in range(nz):
					tx[i][j][k] = afix * (phi[i + 1][j][k] - phi[i - 1][j][k]) \
												+ bfix * (phi[i + 2][j][k] - phi[i - 2][j][k])

		# BCs @ x = nx
		for j in range(ny):
			for k in range(nz):
				tx[-2][j][k] = afix * (phi[-1][j][k] - phi[-3][j][k]) \
											 - bfix * (phi[-2][j][k] + phi[-4][j][k])
				tx[-1][j][k] = -2 * (afix * phi[-2][j][k] + bfix * phi[-3][j][k])

		# Tri-diagonal solve
		for i in range(2, nx):
			for j in range(ny):
				for k in range(nz):
					tx[i][j][k] -= fsx[i] * tx[i - 1][j][k]
		for j in range(ny):
			for k in range(nz):
				tx[-1][j][k] *= fwx[-1]
		for i in range(nx - 2, -1, -1):
			for j in range(ny):
				for k in range(nz):
					tx[i][j][k] -= ffx[i] * tx[i + 1][j][k]
					tx[i][j][k] *= fwx[i]
					
	return tx
		
def derx_2():
	""" X-derivative for Dirichlet BCs. """

	tx = [[[0 for k in range(nz)] for j in range(ny)] for i in range(nx)]

	# BCs @ x = 0
	for j in range(ny):
		for k in range(nz):
			tx[0][j][k] = af1x * phi[0][j][k] + bf1x * phi[1][j][k] + cf1x * phi[2][j][k]
			tx[1][j][k] = af2x * (phi[2][j][k] - phi[0][j][k])

	# Internal nodes
	for i in range(2, nx - 2):
		for j in range(ny):
			for k in range(nz):
				tx[i][j][k] = afix * (phi[i + 1][j][k] - phi[i - 1][j][k]) \
											+ bfix * (phi[i + 2][j][k] - phi[i - 2][j][k])

	# BCs @ x = nx
	for j in range(ny):
		for k in range(nz):
			tx[-2][j][k] = afmx * (phi[-1][j][k] - phi[-3][j][k])
			tx[-1][j][k] = -(afnx * phi[-1][j][k] + bfnx* phi[-2][j][k] + cfnx * phi[-3][j][k])

	# Tri-diagonal solve
	for i in range(1, nx):
		for j in range(ny):
			for k in range(nz):
				tx[i][j][k] -= fsx[i] * tx[i - 1][j][k]
	for j in range(ny):
		for k in range(nz):
			tx[-1][j][k] *= fwx[-1]
	for i in range(nx - 2, -1, -1):
		tx[i][j][k] -= ffx[i] * tx[i + 1][j][k]
		tx[i][j][k] *= fwx[i]

	return tx
	
def derx(phi, nx, ny, nz, nclx, npaire=0):
	""" Take derivatives of phi in the x direction by calling out to derivative functions determined
	by BCs. """

	if nclx == 0:
		derx_0(phi, nx, ny, nz)
	elif nclx == 1:
		derx_1(phi, nx, ny, nz, npaire)
	else:
		derx_2(phi, nx, ny, nz)
