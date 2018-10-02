"""
       FILE: mesh.py
     AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
DESCRIPTION: Provides tools for interacting with xcompact3d rectilinear meshes.
"""

from xcompact3d_tools.tools.errchck import check_domain

def create_xyz(nx, ny, nz, lx, ly, lz, cx=0, cy=0, cz=0, beta=0):
	"""Given size of a (cuboid) computational domain, centred at offset (cx, cy, cz), create the x, y
	and z arrays."""

	check_domain(lx, ly, lz, cx, cy, cz)

	##-----------------------------------------------------------------------------------------------
	# Function start
	dx = lx / (nx - 2.)
	dy = ly / (ny - 2.)
	dz = lz / (nz - 2.)

	x = []
	y = []
	z = []
	for i in range(nx):
		x.append([])
		y.append([])
		z.append([])
		for j in range(ny):
			x[-1].append([])
			y[-1].append([])
			z[-1].append([])
			for k in range(nz):
				x[-1][-1].append((i + 1) * dx - cx)
				y[-1][-1].append((j + 1) * dy - cy)
				z[-1][-1].append((k + 1) * dz - cz)

	return x, y, z, dx, dy, dz
