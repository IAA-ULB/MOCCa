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




class SimpleSelection:
    def __init__(self
            , criteria: dict[str, str | int]
            , files: list[Path]
            , order:str|list[str]|None = None
            ):
        """
        Args:
            criteria: {varname:value} selects all files satisfying f"{varname}={value}" from the `files`
                argument, for all varname/value pairs in the dict
                The process of generating that list of selections is called expanding.
            files: list of files from which to select.

        Returns: a subset of files satisfying criteria as a list. As a side effect the values in the criteria
            are transformed into the corresponding str that must occur in the file name, in order to be selected.
        """
        self.criteria = criteria
        self.expanded = False
        self.from_files = files if isinstance(files, list) else list(files)

        if not order:
            self.order = []
        elif isinstance(order, str):
            self.order = [order]
        else:
            self.order = order

    def is_simple(self) -> bool:
        """We allow for more general criteria where the value corresponding to the varname is actually a list
        or '*' which acts as a wildcard for all occurring values in the filenames.

        Returns:
            False if the criteria are in the above form, True otherwise.
        """
        for varname,value in self.criteria.items():
            if value == '*' or isinstance(value,list):
                return False
        return True

    def select(self, verbose:bool = False):
        assert self.is_simple()
        if verbose:
            print("Selecting entries satisfying:")
            for c in self.criteria.values():
                print(f"    {c}")
        # perform the selection
        self.selection = []
        for file in self.from_files:
            str_filename = str(file.name)
            for key, value in self.criteria.items():
                c_  = (f"{value}" if key == 'x' else f"{key}={value}") + get_next_char(key)
                if not c_ in str_filename:
                    break
            else:
                self.selection.append(file)
                if verbose:
                    print(f" -> {file} : {get_walltime(file)}")

        return self # allow chaining


    def expand(self, verbose=False) -> list:
        if self.is_simple():
            return [self]
        else:
            self._expand_criteria()
            self._sort_criteria()
            expanded = []
            for criteria in self.criteria:
                s = SimpleSelection(criteria, self.from_files).select(verbose=verbose)
                if s:
                    expanded.append(s)
            return expanded

    def __bool__(self):
        try:
            if self.selection:
                return True
            else:
                return False
        except AttributeError:
            return False

    def _expand_criteria(self):
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
                for all currently occurring values in the filenames.
        """
        expand = False
        for varname,value in self.criteria.items():
            if value == '*':
                self.criteria[varname] = get_possible_values(varname, self.from_files)
                expand = True
            elif isinstance(value, list):
                if len(value) > 1:
                    expand = True
            else:
                # put the concrete values also in a list with a single item
                self.criteria[varname] = [value]

        if expand:
            # At least one of the criteria needs expansion.
            # Criteria is a dict where where the values are lists with one
            # or more items.
            # We must turn that into a list of dicts where each value is a single value yielding all possible combinations
            varnames = list(self.criteria.keys())
            values   = list(self.criteria.values())
            mi = MultiIndex(values)
            list_of_criteria = []
            for i in range(len(mi)):
                c = dict(zip(varnames, mi.value()))
                list_of_criteria.append(c)
                mi.increment()

        else:
            list_of_criteria = [self.criteria]


        self.criteria = list_of_criteria

    def is_expanded(self) -> bool:
        return isinstance(self.criteria, dict)

    def _sort_criteria(self):
        """"""
        for varname in self.order:

            if varname == 'x':
                def sort_func(criteria):
                    value = int(criteria[varname].split('x')[0])
                    return value

                self.criteria = sorted(self.criteria, key=sort_func)

            elif varname == 'nnodes':
                def sort_func(criteria):
                    value = int(criteria[varname])
                    return value

                self.criteria = sorted(self.criteria, key=sort_func)

            elif varname == 'napn':
                def sort_func(criteria):
                    value = int(criteria[varname])
                    return value

                self.criteria = sorted(self.criteria, key=sort_func)

            elif varname == 'backend':
                def sort_func(criteria):
                    value = criteria[varname]
                    return value

                self.criteria = sorted(self.criteria, key=sort_func, reverse=True)

            else:
                raise NotImplementedError(f"Don't know how to sort on {varname=}.")


    def __repr__(self):
        try:
            files = f", {self.selection}"
        except AttributeError:
            files = ""
        return f"{self.__class__.__name__}({self.criteria}{files})"

    def get_data_arrays(self,xAxis:str):
        assert self.is_simple()
        if not (hasattr(self, 'xd') and hasattr(self,'yd')):
            pattern = re.compile(f"{xAxis}=(\\d+)")
            self.xd = np.array([ int(pattern.findall(str(filepath.name))[0]) for filepath in self.selection ])
            self.yd = np.array(get_walltime(self.selection))
            # Sort according to x
            idx = np.argsort(self.xd)
            self.xd = self.xd[idx]
            self.yd = self.yd[idx]


identity = lambda x: x

class Expression:
    """
    Expr apply a mathematical function to one or more variables of type list[SimpleSelection].

    """
    def __init__(self, variables:list[SimpleSelection], fun=identity):
        """
        Args:
            variables (list[SimpleSelection)
            expr : a lambda function
        """
        for var in variables:
            assert len(var) == len(variables[0])
        self.variables = variables
        self.fun = fun

    def __call__(self, i:int, xAxis:str):
        pattern = re.compile(f"{xAxis}=(\\d+)")
        xd = np.array([ pattern.findall(str(filepath.name))[0] for filepath in self.variables[0][i][0] ])
        yd = [ np.array(get_walltime(s[i][0])) for s in self.variables ]
        l = [ len(y) for y in yd ]
        lm = min(len(xd), *l)
        lM = max(len(xd), *l)
        if not lm == lM:
            xd = xd[:lm]
            yd = [ y[:lm] for y in yd ]
        yd = self.expr(*yd)
        return xd, yd

    def nGraphs(self):
        """
        Returns:
            the number graphs .
        """
        return len(self.variables[0])

    def nVariables(self):
        """
        Returns:
            the number of variables.
        """
        return len(self.variables)


def plot( arg: Expression|list[SimpleSelection]|SimpleSelection
        , xAxis:str
        , title:str = ''
        , legend:list[str]|None = None
        , plotfun=plt.plot
        ):
    """
    Produce a plot

    Args:
        arg: Expr or list[SimpleSelection]
        xAxis: varname of the feature that must appear in the x-axis
        title: title of the figure
        legend: list of labels for the selections. Well be composed from selection criteria if None.
    """
    if isinstance(arg, SimpleSelection): # list[SimpleSelection]
        return plot([arg], xAxis=xAxis, title=title, legend=legend, plotfun=plotfun)
    elif isinstance(arg, Expression):
        # arg is an expression
        graphs = []
        for j in range(arg.nGraphs()):
            vars_yd = []
            for i in range(arg.nVariables()):
                s = arg.variables[i][j]
                s.get_data_arrays(xAxis)
                if i == 0:
                    xd = s.xd
                else:
                    assert np.all(s.xd == xd)
                vars_yd.append(s.yd)
            yd = expr.fun(*vars_yd)
            graphs.append((xd,yd))
    else:
        # arg is list[SimpleSelections]
        graphs = []
        for s in arg:
            s.get_data_arrays(xAxis=xAxis)
            graphs.append((s.xd, s.yd))


    if legend is None:
        # no legend provided, create from criteria
        legend = []
        for simpleSelection in arg:
            kv = []
            for k, v in simpleSelection.criteria.items():
                kv.append(f"ncores={v}" if k == 'x' else f"{k}={v}")
            legend.append(' & '.join(kv))

    fig, ax = plt.subplots()
    plt.title(title)

    for ig,g in enumerate(graphs):
        xd, yd = g

        plt.xlabel(xAxis)
        plt.ylabel(f"walltime [s]")

        if legend:
            plotfun(xd,yd,'o--', label=legend[ig])
        else:
            plotfun(xd,yd,'o--')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    figure_name = folder/('png')/f'{title}.png' if title else f'{random.randrange(0,1_000_000):06}.png'
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



if __name__ == '__main__':
    add_napn()

    files = list(folder.glob("**/*.out")) # convert to list so it can be used more than once.

    do_all = True
    if do_all:
        title = 'single node scaling (elpa)'
        criteria = {'nnodes':1, 'x':'*', 'backend':'e'}
        selection = SimpleSelection(criteria, files, order='x').expand()
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
        selection = SimpleSelection(criteria, files,order='x').expand()
        plot(selection, xAxis='na', title=title, legend=legend)

        title = 'single node scaling (scalapack vs elpa)'
        criteria = {'nnodes': 1, 'x': ['4x4','8x8','11x11'], 'backend': '*'}
        selection = SimpleSelection(criteria, files, order=['x', 'backend']).expand()
        legend = ['scalapack, 4x4 cores', 'scalapack, 8x8 cores', 'scalapack, 11x11 cores'
                 ,'elpa, 4x4 cores'     , 'elpa, 8x8 cores',      'elpa, 11x11 cores'
                 ]
        plot(selection, xAxis='na', title=title, legend=legend)

        title = 'single node scaling (scalapack over elpa ratio)'
        selection_s = SimpleSelection({'nnodes': 1, 'x': ['8x8','9x9','10x10','11x11'], 'backend': 's'}, files, order='x').expand()
        selection_e = SimpleSelection({'nnodes': 1, 'x': ['8x8','9x9','10x10','11x11'], 'backend': 'e'}, files, order='x').expand()
        legend = ['8x8','9x9','10x10','11x11']
        expr = Expression([selection_s, selection_e],fun= lambda x,y : x/y)
        plot(expr, xAxis='na', title=title, legend=legend)

        title ='multi-node scaling (elpa)'
        criteria = {'nnodes': [2,3,4,8,18,32,50], 'backend': 'e'}
        selection = SimpleSelection(criteria, files, order='nnodes').expand()
        legend = [ get_max_nprows(int(s.criteria['nnodes']), cpus_per_node=128) for s in selection ]
        legend = [ f'{l}x{l} cores = {l*l} ranks' for l in legend]
        plot(selection, xAxis='na', title=title, plotfun=plt.semilogy, legend=legend)

        title ='multi-node scaling (scalapack)'
        criteria = {'nnodes': [2,3,4,8,18,32,50], 'backend': 's'}
        selection = SimpleSelection(criteria, files, order='nnodes').expand()
        plot(selection, xAxis='na', title=title, plotfun=plt.semilogy, legend=legend)

        title ='multi-node scaling (scalapack vs elpa)'
        criteria = {'nnodes': [2,3,4,8,18,32,50], 'backend': '*'}
        selection = SimpleSelection(criteria, files, order=['nnodes', 'backend']).expand()
        legend = []
        for s in selection:
            backend = 'elpa' if s.criteria['backend'] == 'e' else 'scalapack'
            nprows = get_max_nprows(int(s.criteria['nnodes']), cpus_per_node=128)
            legend.append(f'{backend}, {nprows}x{nprows} cores = {nprows*nprows} ranks')
        # legend = [ , cpus_per_node=128) for s in selection ]
        # legend = [ f'{l}x{l} cores' for l in legend]
        plot(selection, xAxis='na', title=title, plotfun=plt.semilogy, legend=legend)

        title ='multi-node scaling (scalapack over elpa ratio)'
        n = len(selection) // 2
        for s in selection:
            s.xd = s.xd[:8]
            s.yd = s.yd[:8]
        expr = Expression([selection[:n],selection[n:]], fun=(lambda s, e : s/e))
        legend = [ l.split(',')[1].strip() for l in legend]
        plot(expr, xAxis='na', title=title, plotfun=plt.plot, legend=legend)

    title = 'strong scaling (elpa) all'
    criteria = {'napn': '*', 'backend': 'e'}
    selection = SimpleSelection(criteria, files, order='napn').expand()
    plot(selection, xAxis='nranks', title=title, plotfun=plt.semilogy)

    title = 'strong scaling (elpa)'
    criteria = {'napn': ['4096','8192','16384','23170','32768','65536'], 'backend': 'e'}
    selection = SimpleSelection(criteria, files, order='napn').expand()
    plot(selection, xAxis='nranks', title=title, plotfun=plt.semilogy)

