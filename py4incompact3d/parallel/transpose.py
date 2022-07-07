import py4incompact3d

import numpy as np

if py4incompact3d.HAVE_DECOMP2D:
    from decomp2d import decomp4py as d4py

def transpose(arr_in, op, arr_to):

    if py4incompact3d.HAVE_DECOMP2D:
        farr_to = np.asfortranarray(arr_to)
        d4py.transpose(np.asfortranarray(arr_in), op, farr_to)
        arr_to = np.ascontiguousarray(farr_to)

    else:
        arr_to = np.copy(arr_in)

    # For some reason need to explicitly return and set
    return arr_to
