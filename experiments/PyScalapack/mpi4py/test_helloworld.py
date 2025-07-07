from mpi4py import MPI

import helloworld

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print(f"({rank}/{size})")
fcomm = MPI.COMM_WORLD.py2f()
helloworld.sayhello(fcomm)
