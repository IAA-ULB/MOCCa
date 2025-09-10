#!/usr/bin/env python3

# move job script with corresponding .out file to folder done

from pathlib import Path
import re
import shutil
import subprocess
import sys


def get_jobid(p:Path) -> str:
    """If the file pointed to by p is a .out file, returns the jobid from the filename,
    otherwise returns None.
    """
    if p.suffix == '.out':
        return str(Path(p.stem).suffix)[1:]

def is_still_running(jobid:str) -> bool:
    result = subprocess.run(['squeue', '-u,' 'entijske'], capture_output=True)
    stdout = result.stdout.decode()
    return jobid in stdout


if __name__ == "__main__":
    """move job script and corresponding .out file to folder 'done'.
    - if '-d' in sys.argv and the corresponding .out file contains the string 'error' (case insensitive)
      the '.out' file is deleted
    - if '-s' in sys.argv and the corresponding .out file does not exist the job script is re-submitted
    - if the job is not yet finished, the output file is left as is.
    """

    error_pattern                       = re.compile(r'error', re.IGNORECASE)
    cancelled_due_to_time_limit_pattern = re.compile(r'DUE TO TIME LIMIT')
    out_of_memory_pattern               = re.compile(r'Out Of Memory')
    path_done = Path('done')
    for path_sh in Path(".").glob("*.sh"):
        path_out = list(Path('.').glob('slurm-' + str(path_sh) + '.*.out'))
        moved = False

        # print(path_sh)
        # print(path_out)
        for p in path_out:
            jobid = get_jobid(p)
            print(f'{jobid=}')
            if is_still_running(jobid):
                print(f"Job {jobid} still running")
                continue
            else:
                with p.open() as f:
                    contents = f.read()
                    if re.search(error_pattern, contents):
                        if re.search(cancelled_due_to_time_limit_pattern, contents):
                            print(f"Job cancelled due to time limit -> deleting '{p}'")
                            p.unlink()
                            print(f"Job cancelled due to time limit -> deleting       '{path_sh}'")
                            path_sh.unlink()
                        if re.search(out_of_memory_pattern, contents):
                            print(f"Job cancelled due to out of memory -> deleting '{p}'")
                            p.unlink()
                            print(f"Job cancelled due to out of memory -> deleting       '{path_sh}'")
                            path_sh.unlink()
                        if '-d' in sys.argv: # delete the file as it contains errors
                            print(f"deleting '{p}' (error)")
                            p.unlink()
                    else:
                        print(f"moving '{p}'")
                        shutil.move(p, path_done)
                        moved = True
                    if moved:
                        print(f"moving '{path_sh}'")
                        shutil.move(path_sh, path_done)
        else:
            if '-s' in sys.argv:
                print(f"resubmitting {path_sh}")
                subprocess.run(["sbatch", str(path_sh)])

                