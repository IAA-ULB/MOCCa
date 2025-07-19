# post-process

import sys
from pathlib import Path
import re

def collect_timings():
    print("collecting:")
    p = Path()
    timings = {}
    pathlist = p.glob('**/*.out')
    for file in pathlist:
        with open(file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                # print(f"{file=}:>{line}")
                if line.startswith("cputime="):
                    cputime = float(line[8:-2])
                    break
        print(f"{file} {cputime=}s")
        timings[str(file)] = cputime

    return timings

def select(timings, kwargs):
    criteria = []
    for key,value in kwargs.items():
        if key=='x':
            s = f"{value[0]}x{value[1]}"
        else:
            s = f"{key}={value}"
        criteria.append(s)

    print(f"Selecting: {criteria}")
    selection = {}
    for filename,cputime in timings.items():
        for c in criteria:
            if not c in filename:
                break
        else:
            selection[filename] = cputime
            print(filename,cputime)

    return selection


if __name__ == '__main__':
    if '--collect' in sys.argv:
        timings = collect_timings()
    selection = select(timings,{'nnodes':1, 'x':(4,4)})
        
        
