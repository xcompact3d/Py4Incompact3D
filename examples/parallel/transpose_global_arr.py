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
from py4incompact3d import parallel

N=2
L=1.0
BC=0

def init():

    nx = N
    ny = 2 * N
    nz = 4 * N
    n = [nx, ny, nz]

    l = [L, L, L]

    bc = [BC, BC, BC]
    
    mesh = py4incompact3d.postprocess.Mesh(
        n = n, l = l, bc = bc
    )

    phiX = py4incompact3d.postprocess.Field()
    phiY = py4incompact3d.postprocess.Field()
    phiZ = py4incompact3d.postprocess.Field()
    phiX.new(mesh, pencil=0)
    phiY.new(mesh, pencil=1)
    phiZ.new(mesh, pencil=2)
    
    return mesh, phiX, phiY, phiZ

def main():

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    parallel.print_par("="*72)
    parallel.print_par("Running transpose example")
    parallel.print_par("-"*72)
    
    mesh, phiX, phiY, phiZ = init()

    # Pre-transpose
    pencil = 0
    for i in range(mesh.NxLocal[pencil]):
        for j in range(mesh.NyLocal[pencil]):
            for k in range(mesh.NzLocal[pencil]):

                idxg = mesh.global_index(pencil, i, j, k)
                phiX.data[0][i][j][k] = idxg
                
    # Transpose X->Y
    parallel.print_par("Testing X-Y transpose")
    pencil = 1
    phiY.data[0] = parallel.transpose(phiX.data[0], "xy", phiY.data[0])
    for i in range(mesh.NxLocal[pencil]):
        for j in range(mesh.NyLocal[pencil]):
            for k in range(mesh.NzLocal[pencil]):

                idxg = mesh.global_index(pencil, i, j, k)
                assert(phiY.data[0][i][j][k] == idxg)
    comm.Barrier()
    parallel.print_par("+ PASSED")
    
    # Transpose Y->Z
    parallel.print_par("Testing Y-Z transpose")
    pencil = 2
    phiZ.data[0] = parallel.transpose(phiY.data[0], "yz", phiZ.data[0])
    for i in range(mesh.NxLocal[pencil]):
        for j in range(mesh.NyLocal[pencil]):
            for k in range(mesh.NzLocal[pencil]):

                idxg = mesh.global_index(pencil, i, j, k)
                assert(phiZ.data[0][i][j][k] == idxg)
    comm.Barrier()
    parallel.print_par("+ PASSED")
                
    # Transpose Z->Y
    parallel.print_par("Testing Z-Y transpose")
    pencil = 1
    phiY.data[0][:][:][:] = 0 # Explicitly zero to ensure data transposed back
    phiY.data[0] = parallel.transpose(phiZ.data[0], "zy", phiY.data[0])
    for i in range(mesh.NxLocal[pencil]):
        for j in range(mesh.NyLocal[pencil]):
            for k in range(mesh.NzLocal[pencil]):

                idxg = mesh.global_index(pencil, i, j, k)
                assert(phiY.data[0][i][j][k] == idxg)
    comm.Barrier()
    parallel.print_par("+ PASSED")
                
    # Transpose Y->X
    parallel.print_par("Testing Y-X transpose")
    pencil = 0
    phiX.data[0][:][:][:] = 0 # Explicitly zero to ensure data transposed back
    phiX.data[0] = parallel.transpose(phiY.data[0], "yx", phiX.data[0])
    for i in range(mesh.NxLocal[pencil]):
        for j in range(mesh.NyLocal[pencil]):
            for k in range(mesh.NzLocal[pencil]):

                idxg = mesh.global_index(pencil, i, j, k)
                assert(phiX.data[0][i][j][k] == idxg)
    comm.Barrier()
    parallel.print_par("+ PASSED")

if __name__ == "__main__":
    main()
