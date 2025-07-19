from jobscript import JobScript, cores_per_node
import math
import sys

if __name__ == "__main__":
    must_submit = '-s' in sys.argv or '--submit' in sys.argv
    print(f"{must_submit=}")

    job = JobScript()
    na = 4096
    job.set('na', na)
    nev = na
    job.set('nev', nev)
    nblk = 32
    job.set('nblk', nblk)
    for nprow in range(4,11):
        nranks = nprow*nprow
        nprow0 = int(math.sqrt(nranks))
        if not nprow==nprow0:
            raise ValueError
        job.set('nranks', nranks)
        job.set('nprows', nprow0)
        job.set('npcols', nprow0)

        nnodes = int(math.ceil(nranks/cores_per_node()))
        job.set('nnodes', nnodes)
        
        job.write(submit=must_submit)
