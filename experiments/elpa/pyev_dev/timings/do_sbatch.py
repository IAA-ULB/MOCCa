import envtools
import jobscript
import math
import sys

sqrt2 = math.sqrt(2)

if __name__ == "__main__":
    must_submit = '-s' in sys.argv or '--submit' in sys.argv
    print(f"{must_submit=}")

    job = jobscript.JobScript()

    if job.template_parameters['cluster'] == 'lumi':
        # check https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/partitions/
        # ['standard', 'standard-g']
        job.add_SBATCH_line('--account=project_465000095')
        job.add_SBATCH_line('--partition=standard')
        job.template_parameters['workspace'] = '/pfs/lustrep4/projappl/project_465000095/entijske'
    
    elif job.template_parameters['cluster'] in ['vaughan', 'breniac']:
        job.add_SBATCH_line('--account=ap_calcua_epicure')
        job.values['workspace'] = '/data/antwerpen/201/vsc20170'
        if job.template_parameters['cluster'] == 'vaughan':
            # check https://docs.vscentrum.be/antwerp/tier2_hardware/vaughan_hardware.html
            # ['zen2','zen3', 'zen3_512', 'ampere_gpu', 'arcturus_gpu']
            job.template_parameters['partition'] = 'zen3' 
    
    job.template_parameters['walltime'] = jobscript.walltime_12h

    na0 = 4096
    naprev = na0
    
    for i in range(1):
        # increase na by sqrt2 every time
        if i%2 == 0: # odd times
            na = int(round(naprev*sqrt2))
        else: # even times
            na = naprev*2
            naprev = na

        job.template_parameters['na'] = na
        nev = na
        job.template_parameters['nev'] = nev
        nblk = 32
        job.template_parameters['nblk'] = nblk
        for nprow in range(4,11):
            nranks = nprow*nprow
            nprow0 = int(math.sqrt(nranks))
            if not nprow==nprow0:
                raise ValueError
            job.template_parameters['nranks'] = nranks
            job.template_parameters['nprows'] = nprow0
            job.template_parameters['npcols'] = nprow0

            nnodes = int(math.ceil(nranks/envtools.get_cpus_per_compute_node(job.template_parameters['cluster'])))
            job.template_parameters['nnodes'] = nnodes
            
            for backend in jobscript.available_backends():
                if backend is jobscript.backend_pyelpa_gpu:
                    # the job script must be adapted.
                    raise NotImplementedError()
                job.template_parameters['backend'] = backend
                job.write(submit=must_submit)
