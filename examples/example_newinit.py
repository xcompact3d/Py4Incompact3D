"""
This file is to test the new json free initialisation
"""

from py4incompact3d.postprocess.postprocess import Postprocess

NX=1441
NY=225
NZ=209
LX=18.0
LY=2.0
LZ=2.0
BCX=1
BCY=2
BCZ=1
DATAPATH="/home/paul/DATA/Incompact3d/lock-exchange/boussinesq"

postproc = Postprocess(n=[NX, NY, NZ],
                       l=[LX, LY, LZ],
                       bc=[BCX, BCY, BCZ])

fieldpath = DATAPATH + "/pp"
postproc.add_field("pp", fieldpath)

postproc.load(time=0)
