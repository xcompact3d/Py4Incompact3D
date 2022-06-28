#!/usr/bin/python3
"""
         FILE: plot_qcrit.py
       AUTHOR: Paul Bartholomew
  DESCRIPTION: Compute the q-criterion using Py4Incompact3D and plot it at a plane.
"""

import sys
sys.path.insert(0,'../../..')
import matplotlib.pyplot as plt

from py4incompact3d.postprocess.postprocess import Postprocess
from py4incompact3d.tools.gradu import calc_gradu
from py4incompact3d.tools.qcrit import calc_qcrit

T=59 # End time
T+=1 # i.e. run up to t=End time
def main():

	# Load data
	input_file = "input.json"
	postprocess = Postprocess(input_file)
	# postprocess.mpi_init()

	for t in range(T):
		print("Postprocessing t = " + str(t))
		postprocess.load(time=t)
		
		# Compute Q-criterion
		print("Computing grad(u)")
		calc_gradu(postprocess)
		# print("Computing vorticity")
		# calc_vort(postprocess)
		print("Computing Q")
		calc_qcrit(postprocess)
		
		# Write to data file
		print("Writing data")
		vel_list = ["ux", "uy", "uz"]
		directions = ["x", "y", "z"]
		# for i in range(3):
		# 	for j in range(3):
		# 		# name = "d" + vel_list[i] + "d" + directions[j]
		# 		# postprocess.write(vars=[name])
		# 		name = "vort" + directions[i] + directions[j]
		# 		postprocess.write(vars=[name])
		postprocess.write(vars=["Q"])

		# Cleanup (save memory)
		postprocess.clear_data()
	
	# # Plot at plane
	# plane = postprocess.plane([postprocess.mesh.Lx / 2,
	# 													 postprocess.mesh.Ly / 2,
	# 													 postprocess.mesh.Lz / 2],
	# 													[1, 0, 0])
	# qcrit = postprocess.interpolate_to(plane, "qcrit")

	# postprocess.mpi_finalise()

if __name__ == "__main__":
	main()

