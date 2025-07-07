import numpy as np
import sys
from pathlib import Path
import json

if __name__ == '__main__':
    if len(sys.argv) == 1:
        slurmdotout = Path('slurm.out')
    else:
        slurmdotout = Path(sys.argv[1])

    timings = {}
    n_doubles = {}
    pathlist = slurmdotout.glob('**/*.out')
    for file in pathlist:
        print(file)
        nranks = None
        ntimings = 0
        with open(file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                # print(f"{file=}:>{line}")
                if line.startswith("Initializing"):
                    words = line.split('x')
                    nprow = int(words[0].split(' ')[1])
                    nranks = nprow*nprow

                elif line[0] == '#':
                    words = line[2:-1].split(':')
                    n = int(words[0].split('x')[0])
                    n += nranks/1000000
                    n_doubles[n] = int(float(words[1].split('=')[1]))

                elif ',' in line:
                    entry = line.split(',')
                    stripped = [ e.strip() for e in entry ]
                    if not nranks:
                        nranks = int(stripped[1])
                    key = ' '.join(stripped[1:4])
                    val = float(stripped[4])
                    if key in timings:
                        timings[key] += val
                    else:
                        timings[key] = val
                    ntimings += 1

    if ntimings != nranks:
        raise RuntimeError(f"{ntimings=} != {nranks=}")
    results = {'timings': timings, 'n_doubles': n_doubles}
    print(results)

    with open(slurmdotout/'results.json', 'w') as f:
        json.dump(results, f, indent=4)