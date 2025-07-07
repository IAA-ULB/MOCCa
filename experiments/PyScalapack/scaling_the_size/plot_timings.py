import json

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import sys
from pathlib import Path


def model2(x, a, b):
    return a * x * x + b * x


def model3(x, a, b, c):
    return a * x * x * x + b * x * x + c * x


if __name__ == '__main__':
    if len(sys.argv) == 1:
        slurmdotout = Path('slurm.out')
    else:
        slurmdotout = Path(sys.argv[1])

    with open(slurmdotout / 'results.json') as f:
        results = json.load(f)

    timings = results['timings']
    # n_doubles = results['n_doubles']

    timings_nranks = {}
    for k, v in timings.items():
        subkeys = [int(e) for e in k.split(' ')]
        n = subkeys[0]
        nranks = subkeys[1]
        nnodes = nranks / 64
        if nranks not in timings_nranks:
            timings_nranks[nranks] = {}
        timings_nranks[nranks][n] = v

    for n1node, d in timings_nranks.items():
        l = len(d)
        nr = []
        dt = []
        for k, v in d.items():
            nr.append(k)
            dt.append(v)
        nr = np.array(nr)
        dt = np.array(dt)
        idx = np.argsort(nr)
        nr = nr[idx]
        dt = dt[idx]
        timings_nranks[n1node] = (nr, dt)

    for nranks, d in timings_nranks.items():
        nprow = int(np.sqrt(nranks))
        title = f"{nprow}x{nprow} processes"
        try:
            x = d[0]
            y = d[1]

            # (a2, b2), pcov = curve_fit(model2, x, y, p0=[0.05, 0.05])
            (a3, b3, c3), pcov = curve_fit(model3, x, y, p0=[0.05, 0.05, 0.05])

            m = 10000*(1+int(np.max(x)/10000))
            xm = np.linspace(np.min(x), m)
            # y2 = model2(xm, a2, b2)
            y3 = model3(xm, a3, b3, c3)

            fig, ax1 = plt.subplots()

            plt.title(title)
            plt.xlabel(f"#processes")
            plt.ylabel(f"walltime (s)")
            plt.plot(x , y, 'bo')
            # plt.plot(xm, y2, 'r--')
            plt.plot(xm, y3, 'g--')
            plt.savefig(slurmdotout / f"{title}.png")
            plt.close()
        except:
            print(f"Failed to create {(slurmdotout / f"{title}.png")}")
            try:
                # Remove the figure as it is not corresponding to the current data
                (slurmdotout / f"{title}.png").unlink()
            except:
                continue
            continue