import numpy as np
import sys
from pathlib import Path
import json
import re

if __name__ == '__main__':
    if len(sys.argv) == 1:
        slurmdotout = Path('slurm.out')
    else:
        slurmdotout = Path(sys.argv[1])

    results = []
    pathlist = slurmdotout.glob('*.out')
    # pattern = r"slurm-job-(\d+)x(\d+)-(\d+)(-TEST)?\.\d\d\d\d\d\d\d\.out"
    for file in pathlist:
        print(file)
        # m = re.match(pattern,str(file.name))
        # if m:
        #     nprow = int(m[1])
        #     n = int(m[3])
        #     nranks = nprow*nprow
        # else:
        #     raise ValueError
        with open(file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                # print(f"{file=}:>{line}")
                if line.startswith(">>"):
                    results.append(line[2:-1])

    print(results)

    with open(slurmdotout/'results.json', 'w') as f:
        json.dump(results, f, indent=4)