import os,sys


class _USE_PYELPA: # for internal use onlyx
    
    def __init__(self, value=None):
        """
        class _USE_PYELPA encapsulates choosing a backend for pyev. 
        
        Args:
            value: If None, select pyscalapack as a backend, unless the environment variable USE_PYELPA is '1'. If True|False selects pyelpa|pyscalapack regardlessof the environment variable USE_PYELPA.

        You can influence the choice programmatically as:

        >...from helpers import use_pyelpa
        >...use_pyelpa.value = False|True # pyscalapack|pyelpa
        """
        if not value is None:
            self.value = bool(value)
            # we neglect sys.argv if given
        else:
            try:
                # print(f"{os.environ['USE_PYELPA']}")
                self.value = os.environ['USE_PYELPA'] == '1'
            except KeyError:
                self.value = False

    def from_argv(self):
        """
        Initialize self.value from the command line arguments. Valid options are:

            --use-pyelpa=0|1 -> False|True (pyev/pyscalapack|pyelpa)
            --use-pyelpa 0|1 -> False|True (pyev/pyscalapack|pyelpa)
            -s               -> False      (pyev/pyscalapack       ) 
            -e               -> True       (                 pyelpa) 
        """
        if len(sys.argv) == 1:
            raise Warning("_USE_PYELPA.from_argv got empty commandline argument list.")
        for i in range(1,len(sys.argv)):
            if sys.argv[i].startswith('--use-pyelpa'):
                if sys.argv[i].startswith('--use-pyelpa='):
                    # '--use-pyelpa=0|1'
                    v = sys.argv[i][13:]
                else:
                    # '--use-pyelpa 0|1'
                    v = sys.argv[i+1]
                if v in ['0','1']:
                    self.value = (v == '1')
                else:
                    raise ValueError(f"--use-pyelpa expecting 0|1, got '{v}'")
                break
            elif sys.argv[i] == '-s': # s for scalapack
                self.value = False
                break
            elif sys.argv[i] == '-e': # s for elpa
                self.value = True
                break
        else:
            print(f"Warning: sys.argv did not set use_pyelpa. Current backend is {self}.")

    def __bool__(self):
        return self.value

    def __str__(self):
        """Return the chosen backend as a string."""
        return 'pyelpa' if self else 'pyev/pyscalapack'
            
    
# create a default _USE_PYELPA object which reads the USE_PYELPA environment variable
# If environment variable USE_PYELPA is not set to 1 pyscalapack is the default backend.
use_pyelpa = _USE_PYELPA()
# If uou want to use pyelpa as the default backend replace the line above with 
#   use_pyelpa = _USE_PYELPA(value=True) 
# of for pyscalapack as the default backend 
#   use_pyelpa = _USE_PYELPA(value=False) 
# Alternatively, you can modify the calling script as
#   from helpers import use_pyelpa
#   use_pyelpa.value = True|False # pyelpa|pyscalapack

def comm2str(comm):
    """
    Returns "rank[{rank}/{size}]" for the mpi4py communicator `comm`.
    """
    return f"rank[{comm.Get_rank()}/{comm.Get_size()}]"

# Do not expose class _USE_PYELPA
__all__ = [
    'use_pyelpa',
    'comm2str',
]