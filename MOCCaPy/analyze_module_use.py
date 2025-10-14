# find all modules and their uses in the files in the `src` dir
# all modules are printed in lower case (Fortran names are case insensitive)
from pathlib import Path
from copy import copy
from collections import Counter
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout

import numpy as np

EXTERNAL_MODULES = ['mpi', 'hdf5', 'iso_fortran_env']

def read_edges_from_src_f90(exclude_external_modules):
    src_dir = (Path(__file__).parent.parent/'src').resolve()
    # print(src_dir)

    module_uses = {}
    for p in src_dir.glob('*.f90'):
        # if p.name == 'fam_testing.f90':
        print(p)
        if p.name == 'run_single.f90':
            print()
        with open(p, 'r') as f90:
            m = None
            lines = f90.readlines()
            for line in lines:
                line = line.strip()
                if line.startswith('module') and not line.startswith('module procedure'):
                    print(line)
                    m = line[6:].strip().split(' ')[0].lower()
                    if not m in module_uses:
                        module_uses[m] = set()
                elif line.startswith('use'):
                    print(line)
                    if not m is None:
                        u = line[3:].strip().split(',')[0].split(' ')[0].lower()
                        if not (exclude_external_modules and u in EXTERNAL_MODULES):
                            module_uses[m].add(u)

            if m is None:
                print('this file does not define a module')
            else:
                print(m, '->', module_uses[m])
            print()
    edges = []
    for m, uses in module_uses.items():
        for u in uses:
            edges.append((m, u))
    return edges

def get_all_right_extensions(connections, edges):
    results = copy(connections)
    done = False
    while not done:
        done = True
        print('*')
        for ic,c in enumerate(results):
            print('+',end='')
            if isinstance(c, tuple):
                c = list(c)
            result = []
            for e in edges:
                # print(connection,e)
                if c[-1] == e[0]:
                    extended = copy(c)
                    extended.append(e[1])
                    result.append(extended)
            if result:
                done = False
                results[ic] = result
            else:
                results[ic] = [c]
        # Flatten results
        flattened = []
        for r in results:
            flattened.extend(r)
        results = flattened
    return results

def get_all_left_extensions(connections, edges):
    results = copy(connections)
    done = False
    while not done:
        done = True
        print('*')
        for ic,c in enumerate(results):
            print('+',end='')
            if isinstance(c, tuple):
                c = list(c)
            result = []
            for e in edges:
                # print(connection,e)
                if c[0] == e[1]:
                    extended = copy(c)
                    extended.insert(0,e[0])
                    result.append(extended)
            if result:
                done = False
                results[ic] = result
            else:
                results[ic] = [c]
        # Flatten results
        flattened = []
        for r in results:
            flattened.extend(r)
        results = flattened
    return results


def unique(results):
    u = set()
    for r in results:
        u.add(tuple(r))
    return u

def get_nodes(edges):
    nodes = set()
    for e in edges:
        nodes.add(e[0])
        nodes.add(e[1])
    return list(nodes)

def encode_edges(edges, nodes):
    e_edges = []
    for e in edges:
        e0 = nodes.index(e[0])
        e1 = nodes.index(e[1])
        e_edges.append([e0,e1])
    return e_edges

def matrix_representation(edges, nodes):


    # empty columns indentify modules at the top of a connection
    # empty rows indentify modules at the bottom of a connection

    n = len(nodes)
    C = np.zeros((n,n),dtype=int)
    for e in edges:
        C[e[0],e[1]] = 1
    return C

def get_connections(C):
    finished = []
    N = C.shape[0]
    # find empty columns
    for icol in range(N):
        column_sum = sum(C[:,icol])
        if column_sum == 0:
            # icol is a starting point
            # find all paths starting with icol
            paths = [[icol]]
            while paths:
                p = paths.pop()
                prev = p[-1]
                next = [j for j in range(N) if C[prev,j] == 1]

                if next:
                    for n in next:
                        pn = p.copy()
                        pn.append(n)
                        paths.append(pn)
                else:
                    finished.append(p)
    return finished


def main(exclude_external_modules=True):
    edges = read_edges_from_src_f90(exclude_external_modules=exclude_external_modules)

    for i,e in enumerate(edges):
        print(i,e)
    edges = list(set(edges))
    for i,e in enumerate(edges):
        print(i,e)
    nodes = get_nodes(edges)
    print(nodes)
    e_edges = encode_edges(edges, nodes)
    for i,e in enumerate(e_edges):
        print(i,e)
    C = matrix_representation(e_edges, nodes)
    N = C.shape[0]
    print(C)
    colsum = np.sum(C,0)
    print(f"{colsum=}")
    for i in range(N):
        name = nodes[i] if colsum[i] == 0 else ''
        print(i,colsum[i],name)
    rowsum = np.sum(C,1)
    print(f"{rowsum=}")
    for i in range(N):
        name = nodes[i] if rowsum[i] == 0 else ''
        print(i,rowsum[i],name)

    print("looking for 'use constants")
    for e in edges:
        if e[1] == 'constants':
            print(e)

    connections = get_connections(C)
    print(connections)
    print(len(connections))
    unique_connections = set()
    for c in connections:
        unique_connections.add(tuple(c))
    print(len(unique_connections))

    histogram = {}
    for c in connections:
        length = len(c)
        if length in histogram:
            histogram[length] += 1
        else:
            histogram[length] = 1
    for k,v in histogram.items():
        print(k,v)
    for c in connections:
        cntr = Counter(c)
        if len(cntr) < len(c):
            print(cntr, c)
    else:
        print('no circular references')

    G = nx.DiGraph()
    G.add_edges_from(edges)
    pos = graphviz_layout(G, prog="dot")
    for k, v in pos.items():
        pos[k] = (-v[1], v[0])

    nx.draw_networkx_nodes(G, pos=pos, node_shape='s', node_size=200,
                           node_color='none', edgecolors='k')
    nx.draw_networkx_edges(G, pos=pos,
                           node_shape='s', width=1, node_size=200)
    nx.draw_networkx_labels(G, pos=pos, font_size=5)

    plt.show()

    # this works but is terribly slow.
    # r = get_all_right_extensions(list(edges), edges)
    # # print(r)
    # r = get_all_left_extensions(r, edges)
    # # print(r)
    # r = unique(r)
    # print(r)

def small():
    # to test get_connections
    edges =[(0,1),(1,2),(2,3),(2,4)]
    nodes = get_nodes(edges)
    C = matrix_representation(edges, nodes)
    print(C)
    connections = get_connections(C)
    print(connections)

if __name__ == '__main__':
    main()
