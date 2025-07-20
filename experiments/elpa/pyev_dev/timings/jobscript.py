import subprocess
import os

walltime_hh  = '00:30:00' 
walltime_4h  =  '4:00:00' 
walltime_12h = '12:00:00'
walltime_24h = '24:00:00'
walltime_72h = '72:00:00'

backend_pyelpa = 'e'
backend_pyelpa_gpu = 'g'
backend_pyscalapack = 's'
backends_cpu = [
    backend_pyelpa, 
    # backend_pyelpa_gpu, 
    backend_pyscalapack,
]
backends_gpu = [
    backend_pyelpa, 
    backend_pyelpa_gpu, 
    backend_pyscalapack,
]

def has_gpu():
    # needs work on lumi
    return False

def available_backends():
    return backends_gpu if has_gpu() else backends_cpu

def get_cluster():
    try:
        cluster = os.environ['VSC_INSTITUTE_CLUSTER']
    except:
        # need implementation for lumi.
        raise NotImplemented('get_cluster(): Undefined $VSC_INSTITUTE_CLUSTER')
    return cluster
    
def get_cores_per_node():
    cpn = {
        'breniac' : 28,
        'vaughan' : 64,
    }
    return cpn[get_cluster()]


_TEMPLATE_VALUES = {
    'jobname' : '',
    'walltime': walltime_hh,
    'nnodes'  :  0,
    'nranks'  :  0,
    'nprows'  :  0,
    'npcols'  :  0,
    'na'      :  0,
    'nev'     :  0,
    'nblk'    : 32,
    'backend' : backend_pyelpa,
    # template variables below depend on the cluster and must only be set on import.
    'cluster' : get_cluster(),
    'cores_per_node' : get_cores_per_node(),
}

_TEMPLATE_JOBSCRIPT = \
"""#!/bin/bash
#SBATCH --job-name {jobname}
#SBATCH --output=slurm-%x.%j.out
#SBATCH --time={walltime}
#SBATCH --mem=0
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes={nnodes} --tasks-per-node={cores_per_node} --cpus-per-task=1
#SBATCH --reservation=rocky9

. /data/antwerpen/201/vsc20170/tantalus_full/experiments/env/ml.sh -p -v

srun -n {nranks} python ev.py {na} {nev} {nblk} -{backend} 
"""

_TEMPLATE_JOBNAME = "job-(na={na},nev={nev},nblk={nblk})-(nnodes={nnodes},nranks={nranks}={nprows}x{npcols})-(cluster={cluster},backend={backend}).sh"

class JobScript:

    def __init__(self):
        self.template_values = _TEMPLATE_VALUES

    def set(self, key, val=None):
        if not key in self.template_values:
            raise KeyError(f"No {key=} in self.template_values.")
        if key == 'jobname':
            val = _TEMPLATE_JOBNAME.format(**self.template_values)
        if val is None:
            raise ValueError(f"Jobscript.set({key=},{val=}) val must not be `None`.")
        self.template_values[key] = val
 
    def jobname(self):
        if not self.template_values['jobname']:
            self.template_values['jobname'] = _TEMPLATE_JOBNAME.format(**self.template_values)
        return self.template_values['jobname']

    def write(self, submit=False):
        self.set('jobname')
        # put the values in the script
        script = _TEMPLATE_JOBSCRIPT.format(**self.template_values)

        # write the job script
        with open(self.jobname(), 'w') as f:
            f.write(script)
            print(f"Wrote job script: {self.jobname()}")
        
        self.template_values['jobname'] = ''
        
        if submit:
            subprocess.run(["sbatch", self.jobname()])
