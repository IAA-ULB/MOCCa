from jobscript import JobScript, get_cores_per_node, available_backends
import math
import sys

sqrt2 = math.sqrt(2)

if __name__ == "__main__":
    must_submit = '-s' in sys.argv or '--submit' in sys.argv
    print(f"{must_submit=}")

    job = JobScript()
    na0 = 4096
    naprev = na0
    
    for i in range(4):
        if i%2 == 0: # odd times
            na = int(round(naprev*sqrt2))
        else: # even times
            na = naprev*2
            naprev = na

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

            nnodes = int(math.ceil(nranks/get_cores_per_node()))
            job.set('nnodes', nnodes)
            
            for backend in available_backends():
                job.set('backend', backend)
                job.write(submit=must_submit)
