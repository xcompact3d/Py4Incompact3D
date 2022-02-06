"""
        FILE: transpose_global_arr.py
 DESCRIPTION: Tests/demonstrates the use of transposing data between parallel pencils.
              For an initial (x) orientation, a field is created so that its values are given by the
              global index at each point, i.e.
                \phi( x_i ) = i
              this can then be transposed into different pencils and values checked for correctness
              based on the global indexing in the current pencil.
"""

import py4incompact3d

N=64
L=1.0
BC=0

def init():

    nx = N
    ny = N
    nz = N
    n = [nx, ny, nz]

    l = [L, L, L]

    bc = [BC, BC, BC]
    
    mesh = py4incompact3d.postprocess.Mesh(
        n = n, l = l, bc = bc
    )

    return mesh

def main():

    mesh = init()

    # Pre-transpose
    pencil = 0
    for i in range(mesh.NxLocal[pencil]):
        for j in range(mesh.NyLocal[pencil]):
            for k in range(mesh.NzLocal[pencil]):

                idxg = mesh.global_index(pencil, i, j, k)
                # TODO: Set field value
                
    # Transpose X->Y
    pencil = 1
    for i in range(mesh.NxLocal[pencil]):
        for j in range(mesh.NyLocal[pencil]):
            for k in range(mesh.NzLocal[pencil]):

                idxg = mesh.global_index(pencil, i, j, k)
                # TODO: Check field value

    # Transpose Y->Z
    pencil = 2
    for i in range(mesh.NxLocal[pencil]):
        for j in range(mesh.NyLocal[pencil]):
            for k in range(mesh.NzLocal[pencil]):

                idxg = mesh.global_index(pencil, i, j, k)
                # TODO: Check field value
                
    # Transpose Z->Y
    pencil = 1
    for i in range(mesh.NxLocal[pencil]):
        for j in range(mesh.NyLocal[pencil]):
            for k in range(mesh.NzLocal[pencil]):

                idxg = mesh.global_index(pencil, i, j, k)
                # TODO: Check field value
                
    # Transpose Y->X
    pencil = 0
    for i in range(mesh.NxLocal[pencil]):
        for j in range(mesh.NyLocal[pencil]):
            for k in range(mesh.NzLocal[pencil]):

                idxg = mesh.global_index(pencil, i, j, k)
                # TODO: Check field value
                
if __name__ == "__main__":
    main()
