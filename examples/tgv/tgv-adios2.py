"""
       FILE: tgv-adios2.py
     AUTHOR: Paul Bartholomew
DESCRIPTION: Tests reading and processing TGV with ADIOS2.
"""

import math

from py4incompact3d.postprocess.postprocess import Postprocess

HDR="=" * 72
LBR="-" * 72

NX=65
NY=NX
NZ=NX

LX=math.pi
LY=LX
LZ=LX

BCX=1
BCY=BCX
BCZ=BCX

ADIOS2_FILE="/home/paul/src/Xcompact3d/Xcompact3d/examples/Taylor-Green-Vortex/data.bp4"

def ke(u, v, w):
    """Given a velocity field, computue the kinetic energy."""
    return 0.5 * (u**2 + v**2 + w**2)

def integrate_func(f, postprocess, time):
    """Given some function f, apply to a field and integrate."""
    pass    
    
def main():

    print(HDR)
    print("Post-processing TGV.")
    print(LBR)

    io_name = "solution-io"

    postprocess = Postprocess(n=[NX, NY, NZ],
                              l=[LX, LY, LZ],
                              bc=[BCX, BCY, BCZ])

    postprocess.init_io(io_name)
    
    postprocess.add_field("ux", ADIOS2_FILE, direction=0)
    postprocess.add_field("uy", ADIOS2_FILE, direction=1)
    postprocess.add_field("uz", ADIOS2_FILE, direction=2)

    postprocess.load(time=1)

    print(postprocess.fields["ux"].data[1])

if __name__ == "__main__":
    main()
