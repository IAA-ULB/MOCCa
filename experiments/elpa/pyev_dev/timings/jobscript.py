import subprocess
import os

walltime_halfhour = '00:30:00' 
walltime_4hours   =  '4:00:00' 
walltime_long     = '12:00:00'
walltime_day      = '24:00:00'

backend_pyelpa = '-e'
backend_pyscalapack = '-s'

def cores_per_node():
    cpn = {
        'breniac' : 28,
        'vaughan' : 64,
    }
    cluster = os.environ['VSC_INSTITUTE_CLUSTER']
    return cpn[cluster]


class JobScript:

    template = \
"""#!/bin/bash
#SBATCH --job-name {jobname}
#SBATCH --output=slurm-%x.%j.out
#SBATCH --time={walltime}
#SBATCH --mem=0
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes={nnodes} --tasks-per-node={cores_per_node} --cpus-per-task=1
#SBATCH --reservation=rocky9

. /data/antwerpen/201/vsc20170/tantalus_full/experiments/env/ml.sh -p -v

srun -n {nranks} python ev.py {na} {nev} {nblk} {backend} 
"""
    jobname_template = "job-(na={na},nev={nev},nblk={nblk})-(nnodes={nnodes},nranks={nranks}={nprows}x{npcols}).sh"
    def __init__(self):
        self.template_values = {
            'jobname' : '',
            'walltime': walltime_4hours,
            'nnodes'  :  0,
            'nranks'  :  0,
            'nprows'  :  0,
            'npcols'  :  0,
            'na'      :  0,
            'nev'     :  0,
            'nblk'    : 32,
            'backend' : backend_pyelpa,
            'cores_per_node' : cores_per_node(),
        }

    def set(self, key, val=None):
        if not key in self.template_values:
            raise KeyError(f"No {key=} in self.template_values.")
        if key == 'jobname':
            val = self.jobname_template.format(**self.template_values)
        if val is None:
            raise ValueError(f"Jobscript.set({key=},{val=}) val must not be `None`.")
        self.template_values[key] = val
 
    
    def jobname(self):
        if not self.template_values['jobname']:
            self.template_values['jobname'] = self.jobname_template.format(**self.template_values)
        return self.template_values['jobname']

    def write(self, submit=False):
        self.set('jobname')
        # put the values in the script
        script = self.template.format(**self.template_values)

        # write the job script
        with open(self.jobname(), 'w') as f:
            f.write(script)
            print(f"Wrote job script: {self.jobname()}")
        
        self.template_values['jobname'] = ''
        
        if submit:
            subprocess.run(["sbatch", self.jobname()])
