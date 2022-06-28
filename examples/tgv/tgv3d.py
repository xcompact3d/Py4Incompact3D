"""
       FILE: tgv3d.py
     AUTHOR: Paul Bartholomew
DESCRIPTION:
"""

import numpy as np

from py4incompact3d.postprocess.postprocess import Postprocess
from py4incompact3d.tools.vort import calc_gradu
from py4incompact3d.tools.vort import calc_vort
from py4incompact3d.tools.qcrit import calc_qcrit

T=20+1

def main():

    postprocess = Postprocess("input3d.json")

    for t in range(1, T):

        print("Postprocessing t = " + str(t))
        
        postprocess.load(time=t)

        calc_gradu(postprocess)
        calc_vort(postprocess)
        calc_qcrit(postprocess)
        postprocess.write(vars=["Q"])
        
        postprocess.clear_data()

if __name__ == "__main__":
    main()
