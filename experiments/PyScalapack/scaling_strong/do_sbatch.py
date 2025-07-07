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
#SBATCH --time=4:00:00
#SBATCH --mem=0
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes={nnodes} --tasks-per-node=64 --cpus-per-task=1

. ../../../vaughan/ml.sh -p -v

mpirun -np {nranks} python scaling_small.py {n} {bf} 
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
    n = 8192
    bf = 8
    for nprow in [9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31]:
        nranks = nprow*nprow
        nprow0 = int(math.sqrt(nranks))
        if not nprow==nprow0:
            raise ValueError
        set_nranks(nranks)
        jobname = f'job-{nprow}x{nprow}-{n}'
        set_jobname(jobname)
        nnodes = int(math.ceil(nranks/64))
        if not nnodes*64 >= nranks:
            raise RuntimeError(f"check nnodes for {jobname}")
        set_nnodes(nnodes)
        set_n(n)
        set_bf(bf)

        make_job_script(submit=SUBMIT)