
from pathlib import Path

from mocca.util import title_line
import pandas as pd

n_iter = 7 # 7 iterations with increasing M per case
def df_from(summary):
    lines = summary.split('\n')[1:-1]
    stencil,ur,reduced,h,M,rmse,mean_diff,max_diff =  [],[],[],[],[],[],[],[]
    for line in lines:
        words = line.split()
        ur       .append(      words[0] )
        stencil  .append(  int(words[1]))
        reduced  .append(      words[2] )
        # bc_scale
        h        .append(      words[4] )
        M        .append(  int(words[5]))
        rmse     .append(float(words[6]))
        mean_diff.append(float(words[7]))
        max_diff .append(float(words[8]))

    data = {
        "stencil"  : stencil,
        "u(r)"     : ur,
        "reduced"  : reduced,
        "h"        : h,
        "M"        : M,
        "RMSE"     : rmse,
        "mean_diff": mean_diff,
        "max_diff" : max_diff
    }
    return pd.DataFrame(data=data)


def get_series(i, lines):
    i0 = (i//n_iter)*n_iter
    i1 = i0+n_iter
    print()
    for i in range(i0,i1):
        print(lines[i])
    return lines[i0:i1]

def find(lines:list, cols:dict, **kwargs):
    for il,line in enumerate(lines):
        for k,v in kwargs.items():
            if line[cols[k]] != v:
                break
        else:
            print(il,line)
            return il, line

if __name__ == "__main__":

    with (Path(__file__).parent / "test_generalized_poisson.py.3D").open() as f:
        lines = f.readlines()[21522-1:] # -1 as python starts counting from 0

    lines = [l.replace(', ', ',') for l in lines]
    for il,line in enumerate(lines):
        words  = line.split()
        words[1] =   int(words[1])
        words[3] = float(words[3])
        words[5] =   int(words[5])
        for iw in range(6,11):
            words[iw] = float(words[iw])
        lines[il] = words

    cols = {
        'u'         :  0,
        'stencil'   :  1,
        'reduced'   :  2,
        'bc_scale'  :  3,
        'h'         :  4,
        'M'         :  5,
        'RMSE'      :  6,
        'mean_diff' :  7,
        'max_diff'  :  8,
        'cput_asmbl':  9,
        'cput_solve': 10,
    }

    # print(lines[ 0])
    # print(lines[-1])
    # print(lines[ 0][cols['h']])
    # print(lines[-1][cols['h']])


    issues = []
    for il, line  in enumerate(lines):
        if il%n_iter == 0:
            line0 = None
        if line0 is not None:
            if line[cols['max_diff']] > line0[cols['max_diff']]:
                print(line)
                issues.append(il)
        line0 = line

    print(f"issue 0")
    print(lines[issues[0]])
    il0,line0 = find(lines, cols, u='xy_gauss', reduced='True', h='(1.3333333333333333,1.3366666666666667,1.3366666666666667)', M=12)
    il0,line0 = find(lines, cols, u='xy_gauss', reduced='True', h='(1.3333333333333333,1.3366666666666667,1.3366666666666667)', M=12, stencil=5)
    il0,line0 = find(lines, cols, u='xy_gauss', reduced='True', h='(1.3333333333333333,1.3366666666666667,1.3366666666666667)', M=12, stencil=5, bc_scale=0.0)
    il0,line0 = find(lines, cols, u='xy_gauss', reduced='(True,False,False)', h='(1.3333333333333333,1.3366666666666667,1.3366666666666667)', M=12, stencil=5, bc_scale=1e20)
    """
         ['xy_gauss', 5, '(True,False,False)',   0.0, '(1.3333333333333333,1.3366666666666667,1.3366666666666667)', 12, 0.062071 , 0.00199602, 0.104096 , 0.00311071, 0.00123104 ]
    1856 ['xy_gauss', 5, '(True,False,False)', 1e+20, '(1.3333333333333333,1.3366666666666667,1.3366666666666667)', 12, 0.0280461, 0.00116746, 0.0217683, 0.0020475 , 0.00595267]
      36 ['xy_gauss', 3, 'True'              , 1e+20, '(1.3333333333333333,1.3366666666666667,1.3366666666666667)', 12, 0.0747728, 0.00662339, 0.134259 , 0.001188  , 0.000322292]
    1828 ['xy_gauss', 5, 'True'              , 1e+20, '(1.3333333333333333,1.3366666666666667,1.3366666666666667)', 12, 0.0140231, 0.00116746, 0.0217683, 0.00144442, 0.000694709]
    1842 ['xy_gauss', 5, 'True'              ,   0.0, '(1.3333333333333333,1.3366666666666667,1.3366666666666667)', 12, 0.0140231, 0.00116746, 0.0217683, 0.00160938, 0.000460833]
['xy_gauss', 5, '(True,False,False)', 0.0, '(2.0,2.005,2.005)'                                         ,  8, 0.104318  , 0.00605206 , 0.098741  , 0.00178829, 0.000426]
['xy_gauss', 5, '(True,False,False)', 0.0, '(1.3333333333333333,1.3366666666666667,1.3366666666666667)', 12, 0.062071  , 0.00199602 , 0.104096  , 0.00311071, 0.00123104]
['xy_gauss', 5, '(True,False,False)', 0.0, '(1.0,1.0025,1.0025)'                                       , 16, 0.0289455 , 0.000636171, 0.0510325 , 0.00540104, 0.00655654]
['xy_gauss', 5, '(True,False,False)', 0.0, '(0.8,0.8019999999999999,0.8019999999999999)'               , 20, 0.00926004, 0.000212014, 0.0126265 , 0.00976442, 0.0286603]
['xy_gauss', 5, '(True,False,False)', 0.0, '(0.6666666666666666,0.6683333333333333,0.6683333333333333)', 24, 0.00499166, 0.000100069, 0.00600385, 0.0153395, 0.0871633]
['xy_gauss', 5, '(True,False,False)', 0.0, '(0.5714285714285714,0.5728571428571428,0.5728571428571428)', 28, 0.00311841, 5.40153e-05, 0.00360763, 0.020308, 0.265661]
['xy_gauss', 5, '(True,False,False)', 0.0, '(0.5,0.50125,0.50125)'                                     , 32, 0.00211374, 3.17874e-05, 0.00245869, 0.0270089, 0.730212]
    """
    get_series(issues[0],lines)
    #
    # for i in issues:
    #     print(title_line(text=f"issue {i}", char='=', width=160), end='')
    #     series = get_series(i,summary_5pt)
    #     for line in series:
    #         print(line)
    #     print(title_line(char='-', width=160), end='')
    #     series = get_series(i,summary_3pt)
    #     for line in series:
    #         print(line)
    #     print(title_line(char='-', width=160))
