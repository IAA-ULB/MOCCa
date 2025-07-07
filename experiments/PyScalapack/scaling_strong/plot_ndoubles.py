import json

import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path
from scipy.optimize import curve_fit, root_scalar


if __name__ == '__main__':
    if len(sys.argv) == 1:
        slurmdotout = Path('slurm.out')
    else:
        slurmdotout = Path(sys.argv[1])

    with open(slurmdotout/'results.json') as f:
        results = json.load(f)

    n_doubles = results['n_doubles']

    n = []
    nd= []
    for k,v in n_doubles.items():
        n.append(float(k))
        nd.append(v)
    n  = np.array(n)
    nd = np.array(nd)
    idx = np.argsort(n)
    x  = n [idx]
    y = nd[idx]


    def lq(x, a, b):
        return a * x * x + b * x

    def lqz(x, a, b):
        return a * x * x + b * x - 536870912


    (a, b), pcov = curve_fit(lq, x, y, p0=[0.05, 0.05])

    xmodel = 1000 * np.arange(30)
    ymodel = lq(xmodel, a, b)
    xlimit = [0, 30000]
    ylimit = [536870912, 536870912]

    fig, ax1 = plt.subplots()
    title = "memory_use"
    plt.title(title)
    plt.xlabel(f"system size nxn")
    plt.ylabel(f"memory use (#doubles)")
    ax1.plot(x, y, "bo-", label="Experiment")
    ax1.plot(xmodel, ymodel, "g--", label="Model")
    ax1.plot(xlimit, ylimit, "r-", label="Limits")
    rr = root_scalar(lqz, args=(a, b), x0=22000, x1=23000)
    z = rr.root
    ax1.plot([z, z], [0, 9e8], 'r-')
    plt.savefig(slurmdotout/f"{title}.png")
    plt.close()