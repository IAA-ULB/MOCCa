import subprocess
import math
import sys

if len(sys.argv) > 1:
    SUBMIT = bool(int(sys.argv[1]))
else:
    SUBMIT = True
print(f"{SUBMIT=}")

template = \
"""#!/bin/bash
#SBATCH --job-name {jobname}
#SBATCH --output=slurm-%x.%j.out
#SBATCH --time=10:00:00
#SBATCH --mem=0
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes={nnodes} --tasks-per-node=64 --cpus-per-task=1

. ../../../vaughan/ml.sh -p -v

mpirun -np {nranks} python pdsyev.py {n} {bf} 
"""

values = {'jobname': ""
    , 'nnodes': 0
    , 'nranks': 0
    , 'n': 0
    , 'bf': 8
    }

def set_jobname(val: str):
    global values
    values['jobname'] = val

def set_nnodes(val: int):
    global values
    values['nnodes'] = val

def set_nranks(val: int):
    global values
    values['nranks'] = val

def set_n(val: int):
    global values
    values['n'] = val

def set_bf(val: int):
    global values
    values['bf'] = val

def make_job_script(submit=True):
    # put the values in the script
    script = template.format(**values)
    # write the job script
    job_sh = f"{values['jobname']}.sh"
    with open(job_sh, 'w') as f:
        f.write(script)
    if submit:
        # submit it
        print(f"About to submit job script: {job_sh}")
        subprocess.run(["sbatch", job_sh])
    else:
        print(f"Wrote job script: {job_sh} (submit==False)")

if __name__ == "__main__":
    bf = 8
    for n1node in [2048,4096,8192]:
        for nnodes in [1,4,9,16]:
            # every node holds n1n0de x n1n0de matrix elements
            n = int(n1node * math.sqrt(nnodes))
            nranks = nnodes*64
            set_nranks(nranks)
            nprow = int(math.sqrt(nranks))
            if not nprow*nprow == nranks:
                raise ValueError
            jobname = f'job-{n1node}xnn-{nprow}x{nprow}'
            set_jobname(jobname)
            set_nnodes(nnodes)
            set_n(n)
            set_bf(bf)

            make_job_script(submit=SUBMIT)