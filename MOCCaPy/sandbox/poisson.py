import numpy as np
import sympy
from tabulate import tabulate

"""Develop experience with solving poisson equationes
1D: check MOCCaPy/literature/FDM.pdf
2D: check MOCCaPy/literature/electronics-11-02365-v2.pdf
3D: generalize the above?
"""
class PoissonSolver1d:
    """Poisson Solver in 1d
        -u" = f

    Remark: we solve the poisson equation on a standard equidistant grid. Not on the
        gridpoints of the LagrangeMesh which has gridpoints on the center of the cells.
    """
    def __init__(self, M, d, rhs_lambda):
        """
        Args:
            M: number of grid points
            d: grid point spacing
            rhs: right hand side of PE (sympy expression
        """
        self.M = M
        self.d = d
        w = M*d
        N = M//2
        self.x = np.linspace((0.5-N)*d, (N+0.5)*d, M, endpoint=False)
        self.rhs = rhs_lambda(self.x)

        self.A = np.zeros((M,M), dtype=np.float64)
        self.b = np.zeros( M   , dtype=np.float64)
        self.u = np.zeros( M   , dtype=np.float64)

    def assemble(self, bcs):
        """Assemble the matrix A and the column vector b."""
        h = self.d
        A = self.A
        b = self.b
        u = self.u
        # assemble b
        b[:] = self.rhs
        b *= h**2

        # assemble A
        A[0,:2]  = np.array([1.,-1.])
        ai = np.array([-1.,2.,-1.])
        for i in range(1,self.M-1):
            A[i,i-1:i+2] = ai
        A[-1,-2:] = np.array([-1.,1.])
        print(f"{A=}")
        # apply BCs
        i0,i1 = 0,self.M-1

        bc_left = bcs['left']
        if 'D' in bc_left:
            u0 = bc_left['D']
            u[0] = u0
            i0 += 1
            b[i0] += u0
        else:
            raise NotImplementedError

        bc_right = bcs['right']
        if 'D' in bc_right:
            u1 = bc_right['D']
            u[i1] = u1
            i1 -= 1
            b[i1] += u1
        else:
            raise NotImplementedError
        print(f"{b=}")

        self.i0, self.i1 = i0, i1

    def solve(self):
        sl = slice(self.i0, self.i1+1)
        u = np.linalg.solve(self.A[sl,sl], self.b[sl])
        self.u[self.i0:self.i1+1] = u
        return self.u


def case_1D():
    # Use sympy to compute a rhs, given an analytical solution
    x = sympy.symbols("x", real=True)
    sigma = 2.
    ue = sympy.exp(-0.5*(x/sigma)**2)  # a gaussian curve
    ue_lambda = sympy.lambdify(x, ue.simplify(), "numpy")
    fe = -ue.diff(x, 2)
    fe_lambda = sympy.lambdify(x, fe.simplify(), "numpy")

    tbl = []
    M = 6
    d = 1.6

    for i in range(10):
        M *= 2
        d /= 2.
        pe1d = PoissonSolver1d(M=M, d=d, rhs_lambda=fe_lambda)
        bcs = {
            'left' : {'D': ue_lambda(pe1d.x[ 0])},
            'right': {'D': ue_lambda(pe1d.x[-1])}
        }
        pe1d.assemble(bcs)

        for i in range(M):
            print(pe1d.A[i,:], pe1d.b[i])

        u = pe1d.solve()

        ue_x = ue_lambda(pe1d.x)
        diff = np.abs(ue_x - u)
        rmse = np.sqrt(np.sum(np.square(diff))/M)
        mean_diff = float(np.mean(diff))
        max_diff = float(np.max(diff))
        print(f"{M=} : {rmse} {mean_diff=} {max_diff=} ")
        tbl.append([M, rmse, mean_diff, max_diff])
        pass
    print(tabulate(tbl,  tablefmt="simple", headers=["M", "RMSE", "Mean diff", "Max diff"]))

if __name__ == "__main__":
    case_1D()