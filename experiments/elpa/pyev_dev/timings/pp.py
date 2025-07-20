# post-process

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import re
import random
import json
from jobscript import assert_exist
from multiindex import MultiIndex

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

def get_values(timings, varname):
    """
    Returns:
        a set of all values found in timings for variable varname.
    """
    assert_exist(varname)

    values = set()
    if varname == 'x':
        pattern = re.compile("(\\d+x\\d+)")
    elif varname == 'cluster':
        pattern = re.compile("cluster=(\\w+),")
    elif varname == 'backend':
        pattern = re.compile("backend=(\\w)\\)")
    else:
        pattern = re.compile(f"{varname}=(\\d+)")
    for filename in timings:
        m = pattern.findall(filename)
        values.add(m[0])
    return values
    
def expand(timngs,criteria):

    def next(mi, criteria_list_of_tuples):
        mi.increment()
        crit = {}
        for d in range(mi.ndims):
            crit[criteria_list_of_tuples[d][0]] = criteria_list_of_tuples[d][1][mi.indx[d]]
        return crit 

    # step 1: replace the '*' with a set of possible values and replace concrete values with {concrete_value}
    for k,v in criteria.items():       
        if v == '*':
            criteria[k] = list(get_values(timings, k))
        else:
            criteria[k] = [v]
    # expand
    # now our criteria have all a list of values as value
    # how do we iterate over something of which we don not now the dimensions? there may be any number of criteria -> MultiIndex.py

    list_tuples = list(criteria.items())
    dims = []
    for d in range(len(list_tuples)):
        dims.append(len(list_tuples[d][1]))
    mi = MultiIndex(dims)
    expanded_criteria = []
    for i in range(len(mi)):
        expanded_criteria.append(next(mi, list_tuples))
    # print(f"{expanded_criteria=}")
    return expanded_criteria
    
def select(timings, criteria:dict, verbose=True, _already_expanded=False):
    """
    Args:
        timings: dict {dotoutfilename:cputime}, data from which to select.
        criteria: {varname:value} selects all timings with f"{varname}={value}" in the dotoutfilename. 
            If value=='*' a list of selections is returned, with one selection for each possible value
    """
    assert criteria, "Error: select(): criteria must not be empty"

    if not _already_expanded:
        # check if we must expand:
        expanded_criteria = []
        for key,value in criteria.items():
            if value == "*":
                expanded_criteria = expand(timings, criteria)
                break
        if expanded_criteria:
            selections = []
            for c in expanded_criteria:
                tpl = select(timings, c, verbose=verbose, _already_expanded=True)
                selections.append(tpl)
            return selections
        
    # already expanded
    criteria_ = []
    for key,value in criteria.items():
        if key=='x':
            s = f"{value[0]}x{value[1]}"
        else:
            assert_exist(key)
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

    nranks_backends = select(timings,{'backend':'*','nranks':'*'})
    plot(nranks_backends, x='na')
    na_backends = select(timings,{'backend':'*','na':'*'})
    plot(na_backends, x='nranks')
    
    # for i in range(2):
    #     for e in backends[i]:
    #         print(e)
    # print(backends)
    # backend_e = select(timings,{'backend':'e'})
    # plot([backend_e,backend_s], x='nranks', title='walltime=f(nranks)')    
    # backends = get_values(timings,'backend')

    # print(f"{get_values(timings,'x')=}")
    # print(f"{get_values(timings,'cluster')=}")
    # print(f"{get_values(timings,'backend')=}")
    # print(f"{get_values(timings,'nranks')=}")

    # criteria = expand(timings, {'nranks':'*', 'backend':'*'})
    # for c in criteria:
    #     print(c)

    

