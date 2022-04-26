import numpy as np

from decomp2d import decomp4py as d4py

def transpose(arr_in, op, arr_to):

    farr_to = np.asfortranarray(arr_to)
    d4py.transpose(np.asfortranarray(arr_in), op, farr_to)
    arr_to = np.ascontiguousarray(farr_to)

    # For some reason need to explicitly return and set
    return arr_to
