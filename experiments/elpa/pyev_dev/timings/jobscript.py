import envtools
import re
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


_TEMPLATE_JOBNAME = "job-(na={na},nev={nev},nblk={nblk})-(nnodes={nnodes},nranks={nranks}={nprows}x{npcols})-(cluster={cluster},backend={backend}).sh"

class JobScript:

    def __init__(self):
        self.script_template = [
            '#!/bin/bash',
            '#SBATCH --job-name {jobname}',
            '#SBATCH --output=slurm-%x.%j.out',
            '#SBATCH --time={walltime}',
            '#SBATCH --mem=0',
            '#SBATCH --account=ap_calcua_epicure',
            '#SBATCH --nodes={nnodes} --tasks-per-node={cores_per_node} --cpus-per-task=1',
            '',
            '. {workspace}/tantalus_full/experiments/env/ml.sh -p -v',
            '',
            'srun -n {nranks} python ev.py {na} {nev} {nblk} -{backend}',
        ]
        # find the index of the first line after the #SBATCH directives.
        self._l = 1
        for line in self.script_template[1:]:
            if line.startswith('#SBATCH'):
                self._l += 1
            else:
                break

        self.jobname_template = "job-(na={na},nev={nev},nblk={nblk})-(nnodes={nnodes},nranks={nranks}={nprows}x{npcols})-(cluster={cluster},backend={backend}).sh"
        
        self.template_parameters = {
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
        }
        # template variables below depend on the cluster and must only be set on import.
        cluster = envtools.get_cluster()
        self.template_parameters['cluster'] = cluster
        self.template_parameters['cores_per_node'] = envtools.get_cpus_per_compute_node(cluster)

    def add_SBATCH_line(self, option, values=[]):
        """
        Add a SBATCH line f'#SBATCH {option}', where the str option may contain 
        template_parameter keys. The values argument may provide a value for each key.
        """
        pattern = re.compile(r"\{\w+\}")
        s = option
        parameter_counter = 0
        m = True
        while m:
            m = re.search(pattern, s)
            if m:
                span = m.span()
                key = s[span[1:-1]]
                if values:
                    self.template_parameters[key] = values[parameter_counter]
                else:
                    self.template_parameters[key] = ''
                s = s[span[1]:]
                parameter_counter += 1
            else:
                break
        
        self.script_template.insert(self._l, f"#SBATCH {option}")

 
    def get_jobname(self):
        self.template_parameters['jobname'] = self.jobname_template.format(**self.template_parameters)
        return self.template_parameters['jobname']


    def get_script(self):
        return "\n".join(self.script_template).format(**self.template_parameters) 


    def write(self, submit=False):
        jobname = self.get_jobname() 
        
        script  = self.get_script()

        with open(jobname, 'w') as f:
            f.write(script)
            print(f"Wrote job script: {jobname}")
        
        if submit:
            subprocess.run(["sbatch", jobname])