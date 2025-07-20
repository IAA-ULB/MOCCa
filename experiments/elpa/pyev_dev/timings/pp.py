# post-process

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import re
import random
import json

folder = Path('slurm.out')

def read_timings():
    timings_json = folder/'timings.json'
    print(f"Reading timings from {timings_json}.")
    with open(timings_json,mode='r') as f:
        timings = json.load(f)
    print(f" -> {len(timings)} entries read.")
    return timings

def collect_timings():
    print("collecting:")
    timings = {}
    pathlist = folder.glob('**/*.out')
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
    
    with open(folder/'timings.json',mode='w') as f:
        json.dump(timings, f, indent=4)

    return timings

def select(timings, criteria:dict, verbose=True):
    """
    Args:
        timings: dict {dotoutfilename:cputime}, data from which to select.
        criteria: {varname:value} selects all timings with f"{varname}={value}" in the dotoutfilename.
    """
    assert criteria, "Error: select(): criteria must not be empty"
    criteria_ = []
    for key,value in criteria.items():
        if key=='x':
            s = f"{value[0]}x{value[1]}"
        else:
            s = f"{key}={value}"
        criteria_.append(s)

    if verbose and criteria_:
        print("Selecting entries satisfying:")
        for c in criteria_:
            print(f"    {c}")

    selection = {}
    for filename,cputime in timings.items():
        for c in criteria_:
            if not c in filename:
                break
        else:
            selection[filename] = cputime
            if verbose:
                print(f" -> {filename} : {cputime}")

    return selection, criteria_

def plot(selections, x=None, title='', legend=None):
    """
    Args:
        selections: list of selection tuples `(selection,criteria)` to plot in a single figure
        x: name of the feature that must appear in the x-axis 
        title: title of the figure
        legend: list of labels for the selections. Well be composed from selection criteria if None.
    """
    if x is None:
        raise ValueError("ERROR: you must provide x.")
    
    if isinstance(selections, tuple):
        selections = [selections]
    
    if len(selections) > 1:
        legend_ = legend
        if legend_ is None:
            # no legend provided, create from criteria
            legend_ = [" & ".join(criteria) for selection, criteria in selections]
                
    fig, ax = plt.subplots()
    plt.title(title)
    for i,tpl in enumerate(selections):
        selection = tpl[0]
        pattern = re.compile(f"{x}=(\\d+)")
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

        plt.xlabel(x)
        plt.ylabel(f"walltime [s]")
        plt.plot(xd,yd,'o--',label=legend_[i])

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    figure_name = folder/f"{title}.png" if title else f"{random.randrange(0,1_000_000):06}.png"
    plt.savefig(figure_name)
    plt.close()
    
if __name__ == '__main__':
    timings = collect_timings() if '--collect' in sys.argv else read_timings()
    # selection, criteria = select(timings,{'nnodes':1, 'x':(4,4)})
    # plot(timings, x='nnodes')        
    backend_s = select(timings,{'backend':'s'})
    backend_e = select(timings,{'backend':'e'})
    plot([backend_e,backend_s], x='nranks', title='walltime=f(nranks)')    

