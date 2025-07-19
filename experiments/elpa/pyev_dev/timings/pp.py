# post-process

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
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

    return selection, criteria

def plot(selection, criteria=None, x=None, title=None):
    if x is None:
        raise ValueError("ERROR: you must provide x.")
    if not criteria and not title:
        raise ValueError("ERROR: criteria and title must not both be None.")
    
    _title = f"{criteria}" if not title else title
    
    pattern = re.compile(f"{x}=(\d+)")
    xd = []
    yd = []
    for filename,cputime in selection.items():
        m = pattern.findall(filename)
        xd.append(int(m[0]))
        yd.append(cputime)
    xd = np.array(xd)
    yd = np.array(yd)
    idx = np.argsort(xd)
    xd = xd[idx]
    yd = yd[idx]

    plt.figure()
    plt.title(title)
    plt.xlabel(x)
    plt.ylabel(f"walltime [s]")
    plt.plot(xd,yd,'o--')
    plt.savefig(f"{_title}.png")
    plt.close()
    
if __name__ == '__main__':
    if '--collect' in sys.argv:
        timings = collect_timings()
    selection, criteria = select(timings,{'nnodes':1, 'x':(4,4)})
    # plot(timings, x='nnodes')        
    plot(timings, x='nranks', title='walltime=f(nranks)')    

