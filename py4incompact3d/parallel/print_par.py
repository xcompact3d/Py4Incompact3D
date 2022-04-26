from mpi4py import MPI

def print_par(msg, comm=MPI.COMM_WORLD, print_rank=0):

    rank = comm.Get_rank()
    if ((rank == print_rank) or (print_rank < 0)):
        print(msg)
