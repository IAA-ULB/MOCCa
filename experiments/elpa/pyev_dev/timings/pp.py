# post-process

import json
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
import random
import re
import sys

from do_sbatch import get_max_nprows
from multiindex import MultiIndex

# where the output files are:
folder = Path('done')

next_char = {'nblk':')'
            ,'nranks':'='
            ,'x':')'
            ,'backend':')'
            }
def get_next_char(varname):
    """We must include the next character when searching for files that match a criterion. E.g.
    a file with 'nnodes=18' would otherwise match criterion 'nnodes=1'.
    """
    return next_char.get(varname, ',')

def get_walltime(files: Path|list[Path]) -> list[float]:
    """Read the walltime directly from the file."""
    if isinstance(files, Path):
        files = [files]
    result = []
    for file in files:
        with open(file) as f:
            lines = f.readlines()
            for line in reversed(lines):
                if line.startswith("cputime="):
                    result.append (float(line[8:-2]))
                    break
            else:
                raise ValueError(f"Could not find cputime in {file}")
    return result


pattern_nprows  = re.compile("(\\d+x\\d+)")         # obtain nprows from substring '={nprows}x{nprows}'

def get_possible_values(varname: str, files: list[Path]):
    """
    Parse all filenames and retrieve the possible values for variable varname, i.e. the values that
    occur in the filenames.

    Args:
        files: list of filepaths
        varname: name of a variable occurring in the filename. varnames appear
            as 'varname=value' in the filename, with the exception of 'x' which indicates 
            nprows and appears as '{nprows}x{nprows}'

    Returns:
        a list of all values found in timings for variable varname.
    """
    values = set()
    pattern = pattern_nprows if varname == 'x' \
        else re.compile(f"{varname}=(\\w+)") # obtain varname from substring 'varname={value}'

    for file in files:
        m = pattern.findall(str(file.name))
        values.add(m[0])

    # convert the set into a list
    return list(values)


def expand(criteria: dict, files: list[Path], ) -> list[dict]:
    """
    Check if the criteria must be expanded, i.e. if one or more criteria contain
      - a list of length greater than 1
      - an asterisk, which is shorthand for all currently occuring values in the filenames.
    Expand if needed.

    Args:
        criteria: { varname1: value                  # no expansion required
                  , varname2: '*'                    # expansion required
                  , varname3: [value1, value2, ...]  # expansion required
                  , ...
                  }]
        files: list of filepaths, optional, only needed if one of the criteria uses '*', which stands
            for all currently occuring values in the filenames.
    """
    expand = False
    for varname,value in criteria.items():
        if value == '*':
            criteria[varname] = get_possible_values(varname, files)
            expand = True
        elif isinstance(value, list):
            if len(value) > 1:
                expand = True
        else:
            # put the concrete values also in a list with a single item
            criteria[varname] = [value]

    if expand:
        # At least one of the criteria needs expansion. criteria is a dict where where the values are lists with one
        # or more items.
        # We must turn that into a list of dicts where each value is a single value yielding all possible combinations
        varnames = list(criteria.keys())
        values   = list(criteria.values())
        mi = MultiIndex(values)
        list_of_criteria = []
        for i in range(len(mi)):
            c = dict(zip(varnames, mi.value()))
            list_of_criteria.append(c)
            mi.increment()
    else:
        list_of_criteria = [criteria]

    return list_of_criteria

def select(criteria: dict[str,str|int], files: list[Path], verbose=True) -> list[Path]:
    """
    Return the subset of files satisfying criteria.

    Args:
        criteria: {varname:value} selects all timings with f"{varname}={value}" in the dotoutfilename.
            If value=='*' a list of selections is returned, with one selection for each possible value.
            The process of generating that list of selections is called expanding.
        files: list of files from which to select.

    Returns: a subset of files satisfying criteria as a list. As a side effect the values in the criteria
        are transformed into the corresponding str that must occur in the file name, in order to be selected.
    """
    assert criteria, "Error: select(): criteria must not be empty"

    if not isinstance(files, list):
        files = list(files)

    list_of_criteria = expand(criteria, files)

    list_of_selections = []
    for criteria in list_of_criteria:
        # Transform (key,value)-pair criteria into the corresponding str that must occur in the filename
        # We do this by overwriting the value in the dict
        for key,value in criteria.items():
            criteria[key] = f"{value}" if key=='x' else f"{key}={value}"

        if verbose:
            print("Selecting entries satisfying:")
            for c in criteria.values():
                print(f"    {c}")

        # perform the selection
        selection = []
        for file in files:
            str_filename = str(file.name)
            for varname,c in criteria.items():
                c_ = c + get_next_char(varname)
                if not c_ in str_filename:
                    break
            else:
                selection.append(file)
                if verbose:
                    print(f" -> {file} : {get_walltime(file)}")
        # It may happen that the selection is empty (because sometimes not all combinations are there)
        if selection:
            list_of_selections.append((selection, criteria))

    return list_of_selections

def plot(list_of_selection_tuples: list[tuple[list[Path],dict]], xAxis, title='', legend=None, plotfun=plt.plot):
    """
    Produce a plot

    Args:
        list_of_selection_tuples: list of (selection, criteria) tuples, Every tuple yields a graph in the figurs
        xAxis: varname of the feature that must appear in the x-axis
        title: title of the figure
        legend: list of labels for the selections. Well be composed from selection criteria if None.
    """

    # if isinstance(selection[0], list):
    #     assert isinstance(criteria, list) and len(selection) == len(criteria)
    # else:
    #     selection = [selection]
    #     criteria  = [criteria]

    if len(selection) > 1:
        legend_ = legend
        if legend_ is None:
            # no legend provided, create from criteria
            legend_ = [" & ".join(tpl[1].values()) for tpl in list_of_selection_tuples]
                
    fig, ax = plt.subplots()
    plt.title(title)
    for i,tpl in enumerate(selection):
        pattern = re.compile(f"{xAxis}=(\\d+)")
        xd = []
        for filepath in tpl[0]:
            m = pattern.findall(str(filepath.name))
            xd.append(int(m[0]))
        xd = np.array(xd)
        yd = np.array(get_walltime(tpl[0]))
        idx = np.argsort(xd)
        xd = xd[idx]
        yd = yd[idx]

        plt.xlabel(xAxis)
        plt.ylabel(f"walltime [s]")

        if len(selection) > 1:
            plotfun(xd,yd,'o--', label=legend_[i])
        else:
            plotfun(xd,yd,'o--')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    figure_name = folder/f"{title}.png" if title else f"{random.randrange(0,1_000_000):06}.png"
    print(f"Writing figure {figure_name}.")
    plt.savefig(figure_name)
    plt.close()


def add_napn():
    """
    Check all .out files for the appearance of 'napn=', if not add it (nanp = na/nnodes).
    (We add a computed variable napn to the file names to be able to select on the computed variable)
    """
    files = list(folder.glob("**/*.out"))
    napn_pattern = re.compile("napn=(\\d+)")
    na_pattern = re.compile("na=(\\d+)")
    nnodes_pattern = re.compile("nnodes=(\\d+)")
    for file in files:
        str_filename = str(file.name)
        m = napn_pattern.findall(str_filename)
        if not m:
            m = na_pattern.findall(str_filename)
            na = int(m[0])
            m = nnodes_pattern.findall(str_filename)
            nnodes = int(m[0])
            napn = na // nnodes
            str_filename_new = str_filename[0:11]+f"{napn=}," + str_filename[11:]
            print(str_filename, '->')
            print(str_filename_new)
            os.rename(file, folder/str_filename_new)
        else:
            print(str_filename)
        print()

def sort_selection(selection: list[tuple[list[Path],dict]],varname:str):
    """"""
    if varname == 'x':
        def sort_func(tpl):
            value = int(tpl[1][varname].split('x')[0])
            return value
        sorted_selection = sorted(selection, key=sort_func)
    elif varname == 'nnodes':
        def sort_func(tpl):
            value = int(tpl[1][varname][7:])
            return value
        sorted_selection = sorted(selection, key=sort_func)
    elif varname == 'napn':
        def sort_func(tpl):
            value = int(tpl[1][varname][5:])
            return value
        sorted_selection = sorted(selection, key=sort_func)
    elif varname == 'backend':
        def sort_func(tpl):
            value = tpl[1][varname]
            return value
        sorted_selection = sorted(selection, key=sort_func, reverse=True)
    else:
        raise ValueError(f"Don't know how to sort on {varname=}.")

    return sorted_selection


if __name__ == '__main__':
    add_napn()

    files = list(folder.glob("**/*.out")) # convert to list so it can be used more than once.

    do_all = False
    if do_all:
        title = 'single node scaling (elpa)'
        criteria = {'nnodes':1, 'x':'*', 'backend':'e'}
        selection = select(criteria, files)
        selection = sort_selection(selection, 'x')
        legend = ['4x4 cores'
                 ,'5x5 cores'
                 ,'6x6 cores'
                 ,'7x7 cores'
                 ,'8x8 cores'
                 ,'9x9 cores'
                 ,'10x10 cores'
                 ,'11x11 cores'
                 ]
        plot(selection, xAxis='na', title=title, legend=legend)

        title = 'single node scaling (scalapack)'
        criteria = {'nnodes':1, 'x':'*', 'backend':'s'}
        selection = select(criteria, files)
        selection = sort_selection(selection, 'x')
        plot(selection, xAxis='na', title=title, legend=legend)

        title = 'single node scaling (scalapack vs elpa)'
        criteria = {'nnodes': 1, 'x': ['4x4','8x8','11x11'], 'backend': '*'}
        selection = select(criteria, files)
        selection = sort_selection(selection, 'x')
        selection = sort_selection(selection, 'backend')
        legend = ['scalapack, 4x4 cores', 'scalapack, 8x8 cores', 'scalapack, 11x11 cores'
                 ,'elpa, 4x4 cores'     , 'elpa, 8x8 cores',      'elpa, 11x11 cores'
                 ]
        plot(selection, xAxis='na', title=title, legend=legend)

        title ='scaling elpa'
        criteria = {'nnodes': [2,3,4,8,18,32], 'backend': 'e'}
        selection = select(criteria, files)
        selection = sort_selection(selection, 'nnodes')
        legend = [ get_max_nprows(int(s[1]['nnodes'].split('=')[1]), cpus_per_node=128) for s in selection ]
        legend = [ f'{l}x{l} cores' for l in legend]
        plot(selection, xAxis='na', title=title, plotfun=plt.semilogy, legend=legend)

        title ='scaling scalapack'
        criteria = {'nnodes': [2,3,4,8,18,32], 'backend': 's'}
        selection = select(criteria, files)
        selection = sort_selection(selection, 'nnodes')
        plot(selection, xAxis='na', title=title, plotfun=plt.semilogy, legend=legend)

        title ='scaling elpa vs scalapack'
        criteria = {'nnodes': [2,3,4,8,18,32], 'backend': '*'}
        selection = select(criteria, files)
        selection = sort_selection(selection, 'nnodes')
        selection = sort_selection(selection, 'backend')
        legend = []
        for s in selection:
            backend = 'elpa' if s[1]['backend']=='e' else 'scalapack'
            nprows = get_max_nprows(int(s[1]['nnodes'].split('=')[1]), cpus_per_node=128)
            legend.append(f'{backend}, {nprows}x{nprows} cores')
        # legend = [ , cpus_per_node=128) for s in selection ]
        # legend = [ f'{l}x{l} cores' for l in legend]
        plot(selection, xAxis='na', title=title, plotfun=plt.semilogy, legend=legend)

    title = 'strong scaling elpa all'
    criteria = {'napn': '*', 'backend': 'e'}
    selection = select(criteria, files)
    selection = sort_selection(selection, 'napn')
    plot(selection, xAxis='nranks', title=title, plotfun=plt.semilogy)

    title = 'strong scaling elpa'
    criteria = {'napn': ['4096','8192','16384','23170','32768','65536'], 'backend': 'e'}
    selection = select(criteria, files)
    selection = sort_selection(selection, 'napn')
    plot(selection, xAxis='nranks', title=title, plotfun=plt.semilogy)

