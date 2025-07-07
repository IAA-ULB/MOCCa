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

def get_rank(entries):
    return int(entries[0])

def get_nranks(entries):
    return int(entries[1])

def get_nprow(entries):
    return int(entries[2])

def get_npcol(entries):
    return int(entries[3])

def get_n(entries):
    return int(entries[4])

def get_bf(entries):
    return int(entries[5])

def get_wt(entries):
    return float(entries[6])

def get_bDP(entries):
    return int(entries[7])

def line2dict(line):
    entries = [e.strip() for e in line.split(',')]

    # print(entries)
    return {
        'rank'  :   int(entries[0]),
        'nranks':   int(entries[1]),
        'nprow' :   int(entries[2]),
        'npcol' :   int(entries[3]),
        'n'     :   int(entries[4]),
        'bf'    :   int(entries[5]),
        'wt'    : float(entries[6]),
        'nDP'   : int(float(entries[7])),
    }

def tuples2array(tuples):
    x = []
    y = []
    for t in tuples:
        x.append(t[0])
        y.append(t[1])
    x = np.array(x)
    y = np.array(y)
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]
    return (x,y)

def select(results, conditions, kx, ky):
    """Select all data from results that satisfy conditions and return two arrays with kx and ky entry
    """
    tuples = []
    for line in results:
        result = line2dict(line)
        keep = True
        for k,v in conditions.items():
            if not result[k] == v:
                keep = False
                break
        if keep:
            tuples.append((result[kx],result[ky]))

    return tuples2array(tuples)



if __name__ == '__main__':
    if len(sys.argv) == 1:
        slurmdotout = Path('slurm.out')
    else:
        slurmdotout = Path(sys.argv[1])

    with open(slurmdotout / 'results.json') as f:
        results = json.load(f)
    
    for nprow in [4,8,16,24,32,40]:
        nranks = nprow*nprow
        print(f"{nranks=}")
        
        kx = 'n'
        ky = 'wt'
        x,y = select(results, {'nranks':nranks}, kx=kx, ky=ky)
        (a3, b3, c3), pcov = curve_fit(model3, x, y, p0=[0.05, 0.05, 0.05])

        xm = np.linspace(np.min(x), np.max(x))
        # y2 = model2(xm, a2, b2)
        y3 = model3(xm, a3, b3, c3)

        fig, ax1 = plt.subplots()
        title = f"{nprow}x{nprow} = {nranks} processes"
        plt.title(title)
        plt.xlabel(f"{kx}")
        plt.ylabel(f"walltime (s)")
        plt.plot(x , y, 'bo')
        # plt.plot(xm, y2, 'r--')
        plt.plot(xm, y3, 'g--')
        savefile = slurmdotout / f"{title}.png"
        print(f"saving figure {savefile}")
        plt.savefig(savefile)
        plt.close()

    fig, ax1 = plt.subplots()
    for nprow in [4,8,16,24,32,40]:
        nranks = nprow*nprow
        print(f"{nranks=}")
        
        kx = 'n'
        ky = 'wt'
        x,y = select(results, {'nranks':nranks}, kx=kx, ky=ky)
        (a3, b3, c3), pcov = curve_fit(model3, x, y, p0=[0.05, 0.05, 0.05])

        xm = np.linspace(np.min(x), np.max(x))
        # y2 = model2(xm, a2, b2)
        y3 = model3(xm, a3, b3, c3)

        plt.title(title)
        plt.xlabel(f"{kx}")
        plt.ylabel(f"walltime (s)")
        plt.plot(x , y, 'bo')
        # plt.plot(xm, y2, 'r--')
        plt.plot(xm, y3, 'g--')

    savefile = slurmdotout / f"all.png"
    print(f"saving figure {savefile}")
    plt.savefig(savefile)
    plt.close()
