#!/usr/bin/python3
"""
         FILE: plot_qcrit.py
       AUTHOR: Paul Bartholomew
  DESCRIPTION: Compute the q-criterion using Py4Incompact3D and plot it at a plane.
"""

import sys
sys.path.insert(0,'../../..')
import matplotlib.pyplot as plt

from Py4Incompact3D.postprocess.postprocess import Postprocess
from Py4Incompact3D.tools.vort import calc_gradu
from Py4Incompact3D.tools.vort import calc_vort
from Py4Incompact3D.tools.qcrit import calc_qcrit


T=35 # End time
T+=1 # i.e. run up to t=End time
T=-1
def main():

    # Load data
    input_file = "input.json"
    postprocess = Postprocess(input_file)
    t=8
    # postprocess.mpi_init()

    print("Postprocessing t = " + str(t))
    postprocess.load(time=t)

    # Compute Q-criterion
    print("Computing grad(u)")
    calc_gradu(postprocess)
    print("Computing vorticity")
    calc_vort(postprocess)
    print("Computing Q")
    calc_qcrit(postprocess)

    # Write to data file
    print("Writing data")
    vel_list = ["ux", "uy", "uz"]
    directions = ["x", "y", "z"]
    postprocess.write(vars=["Q"])

    # Cleanup (save memory)
    postprocess.clear_data()

if __name__ == "__main__":
	main()

