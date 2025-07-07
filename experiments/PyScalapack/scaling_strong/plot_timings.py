import json

import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path

if __name__ == '__main__':
    if len(sys.argv) == 1:
        slurmdotout = Path('slurm.out')
    else:
        slurmdotout = Path(sys.argv[1])

    with open(slurmdotout/'results.json') as f:
        results = json.load(f)

    timings = results['timings']
    # n_doubles = results['n_doubles']

    timings_n = {}
    for k,v in timings.items():
        subkeys = [int(e) for e in k.split(' ')]
        nr= subkeys[0]
        n = subkeys[1]
        if n not in timings_n:
            timings_n[n] = {}
        timings_n[n][nr] = v

    for n, d in timings_n.items():
        l = len(d)
        nr = []
        dt = []
        avg = 1e16
        for k,v in d.items():
            if v < 10*avg:
                nr.append(k)
                dt.append(v)
                avg = sum(dt)/len(dt)
            else:
                print(f"outlier excluded ({k},{v})")
        nr = np.array(nr)
        dt = np.array(dt)
        idx = np.argsort(nr)
        nr = nr[idx]
        dt = dt[idx]
        timings_n[n] = (nr,dt)

    for n,d in timings_n.items():
        plt.figure()
        title = f"{n=}"
        plt.title(title)
        plt.xlabel(f"#processes")
        plt.ylabel(f"walltime (s)")
        plt.plot(d[0],d[1],'o--')
        plt.savefig(slurmdotout/f"{title}.png")
        plt.close()