import envtools
import jobscript
import math
import sys

sqrt2 = math.sqrt(2)

def get_max_nprows(nnodes:int, cpus_per_node:int|None = None, verbose:bool = False) -> int:
    """Compute the largest int that if squared is <= nnodes*ncores"""
    cpn = cpus_per_node if cpus_per_node else envtools.get_cpus_per_compute_node()
    ncores = nnodes*cpn
    nprows = int(math.floor(math.sqrt(ncores)))
    npr = nprows
    while True:
        npr += 1
        if npr*npr > ncores:
            break
        nprows = npr
    if verbose:
        print(f"{nprows*nprows} <= {ncores} < {npr*npr}")
    return nprows

# for nnodes in [2,8,18,32,50,72, 2*n^2] all cores are in use on lumi (as it has 128 cores per node )
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
    
    nnodes = int(sys.argv[1])

    na = 4096 if nnodes <= 8 else 4096*4

    if nnodes==1:
        try:
            npr = int(sys.argv[2])
        except:
            npr = 4
            
    naprev = na
    for i in range(16):
        
        job.template_parameters['na'] = na
        nev = na
        job.template_parameters['nev'] = nev
        nblk = 32
        job.template_parameters['nblk'] = nblk

        if nnodes == 1:
            # On a single node we go from 4x4, 5x5, ...
            npr_max = get_max_nprows(nnodes=1)
            for nprow in range(npr,npr_max+1):
                nranks = nprow*nprow
                nprow0 = int(math.sqrt(nranks))
                if not nprow==nprow0:
                    raise ValueError
                job.template_parameters['nranks'] = nranks
                job.template_parameters['nprows'] = nprow0
                job.template_parameters['npcols'] = nprow0
                job.template_parameters['nnodes'] = nnodes
                
                for backend in jobscript.available_backends():
                    if backend is jobscript.backend_pyelpa_gpu:
                        # the job script must be adapted.
                        raise NotImplementedError()
                    job.template_parameters['backend'] = backend
                    job.write(submit=must_submit)
        else:
            nprows = get_max_nprows(nnodes, cpus_per_node=envtools.get_cpus_per_compute_node(), verbose=True)
            nranks = nprows*nprows
            job.template_parameters['nranks'] = nranks
            job.template_parameters['nprows'] = nprows
            job.template_parameters['npcols'] = nprows
            job.template_parameters['nnodes'] = nnodes
            for backend in jobscript.available_backends():
                if backend is jobscript.backend_pyelpa_gpu:
                    # the job script must be adapted.
                    raise NotImplementedError()
                job.template_parameters['backend'] = backend
                job.write(submit=must_submit)
            
        # increase na by sqrt2 every time
        if i%2 == 0: # odd times
            na = int(round(naprev*sqrt2))
        else: # even times
            na = naprev*2
            naprev = na
        if na > 50000*nnodes:
            break
