import numpy as np
import sympy
from contourpy.types import offset_dtype
from tabulate import tabulate
import scipy

from mocca.mesh import LagrangeMesh

"""Develop experience with solving poisson equations
1D: check MOCCaPy/literature/FDM.pdf
2D: check MOCCaPy/literature/electronics-11-02365-v2.pdf
3D: generalize the above?
https://leifh.folk.ntnu.no/teaching/tkt4140/._main055.html
"""
class PoissonSolver1d:
    """Poisson Solver in 1d with Dirichlet boundary conditions
        -u" = f
    """
    def __init__(self, mesh, rhs_lambda, bc_scale=None):
        """P
        Args:
            mesh: Lagrange mesh
            rhs_lamba: right hand side of PE as a numpy `ufunc`, to be applied to the gridpoints.
            bc_scale: if None, boundary conditions are resolved by setting $A_{ij}=\delta_{ij} and
                $b_i=d_i$ where $d_i$ is the value of the dirichlet boundary condition at point $i$.
                This approach destroys the symmetry of $A$. If you want to keep the symmetry of $A$,
                set `bc_scale to a large number $\phi$. The boundary conditions are then resolved by
                setting $A_{ii}$ `+=` $\phi$ and $b_i$ `*=` $d_i \phi$.
                Both methods are described in MOCCaPy/literature/FDM.pdf.
        """
        assert mesh.dim == 1
        self.mesh = mesh
        M = mesh.M[0]
        h = mesh.d[0]
        self.x = np.empty(M+2, dtype=np.float64)
        self.x[1:-1] = mesh.g1D[0]
        self.x[0] = self.x[1] - h
        self.x[-1] = self.x[-2] + h
        self.b = rhs_lambda(self.x)
        if bc_scale is None:
            A = scipy.sparse.linalg.LaplacianNd((M+2,), boundary_conditions='dirichlet').tosparse()
        else:
            A = scipy.sparse.linalg.LaplacianNd((M+2,), boundary_conditions='dirichlet', dtype=np.float64).tosparse()

        A.data *=-1 # because we solve -u" = f
        self.bc_scale = bc_scale
        if bc_scale is None:
            A.data[0,-2] = 0
            A.data[1, 0] = 1
            A.data[1,-1] = 1
            A.data[2, 1] = 0
        else:
            A.data[1, 0] += bc_scale
            A.data[1,-1] += bc_scale
        self.A = A

        # print(self.A.toarray())
        self.u = np.zeros(M+2, dtype=np.float64)

    def assemble(self, bcs):
        """Assemble the matrix A and the column vector b."""
        h = self.mesh.d[0]

        # assemble b
        b = self.b
        b *= h**2

        # apply BCs
        A = self.A
        bc_left = bcs['left']
        if 'D' in bc_left:
            u0 = bc_left['D']
            if self.bc_scale is None:
                b[0] = u0
            else:
                b[0] = u0*self.bc_scale
        else:
            raise NotImplementedError

        bc_right = bcs['right']
        if 'D' in bc_right:
            u1 = bc_right['D']
            if self.bc_scale is None:
                b[-1] = u1
            else:
                b[-1] = u1 * self.bc_scale
        else:
            raise NotImplementedError
        print(b)

    def solve(self):
        self.u = scipy.sparse.linalg.spsolve(self.A, self.b)
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
    for iter in range(10):
        M *= 2
        M -= 2 # to compare with test_generalized_poisson, poisson-bis adds two point to the grid.
        d /= 2.
        mesh = LagrangeMesh(M, d, dim=1, reduced=False)

        for bc_scale in [None, 1e20]:
            pe1d = PoissonSolver1d(mesh, rhs_lambda=fe_lambda, bc_scale=bc_scale)
            bcs = {
                'left' : {'D': ue_lambda(pe1d.x[ 0])},
                'right': {'D': ue_lambda(pe1d.x[-1])}
            }
            pe1d.assemble(bcs)

            if iter == 0:
                # create a table to compare with poisson-bis
                t = []
                A = pe1d.A.toarray()
                for i in range(M+2):
                    ti = [i]
                    for j in range(M+2):
                        ti.append(A[i,j])
                    ti.append(pe1d.b[i])
                    t.append(ti)

            u = pe1d.solve()
            ue_x = ue_lambda(pe1d.x)
            diff = np.abs(ue_x - u)

            if iter == 0:
                for i in range(M+2):
                    t[i].extend([u[i], ue_x[i], diff[i]])
                headers = ['i']
                headers.extend([str(i) for i in range(M+2)])
                headers.extend(['f', 'u', 'ue_x', 'diff'])
                print(tabulate(t, tablefmt="simple", headers=headers))

            rmse = np.sqrt(np.sum(np.square(diff))/M)
            mean_diff = float(np.mean(diff))
            max_diff = float(np.max(diff))
            print(f"{bc_scale=} {M=} : {rmse} {mean_diff=} {max_diff=} ")
            tbl.append([bc_scale, M, rmse, mean_diff, max_diff])
            pass

    print(tabulate(tbl,  tablefmt="simple", headers=["bc_scale", "M", "RMSE", "Mean diff", "Max diff"]))


class PoissonSolverNd:
    """Poisson Solver in 2d with Dirichlet boundary conditions
        -u" = f
    """
    def __init__(self, mesh, rhs_lambda, bc_scale=None):
        """P
        Args:
            mesh: Lagrange mesh
            rhs_lamba: right hand side of PE as a numpy `ufunc`, to be applied to the gridpoints.
            bc_scale: if None, boundary conditions are resolved by setting $A_{ij}=\delta_{ij} and
                $b_i=d_i$ where $d_i$ is the value of the dirichlet boundary condition at point $i$.
                This approach destroys the symmetry of $A$. If you want to keep the symmetry of $A$,
                set `bc_scale to a large number $\phi$. The boundary conditions are then resolved by
                setting $A_{ii}$ `+=` $\phi$ and $b_i$ `*=` $d_i \phi$.
                Both methods are described in MOCCaPy/literature/FDM.pdf.
        """
        if mesh.dim == 3:
            raise NotImplementedError

        self.mesh = mesh
        g1dp = mesh.dim * [0]
        domain = mesh.dim * [0]
        for idim in range(mesh.dim):
            g1dp[idim] = np.empty(mesh.M[0] + 2, dtype=np.float64)
            g1dp[idim][1:-1] = mesh.g1D[0]
            g1dp[idim][ 0] = g1dp[idim][ 1] - mesh.d[idim]
            g1dp[idim][-1] = g1dp[idim][-2] + mesh.d[idim]
            domain[idim] = mesh.M[idim] + 2

        # g1dp = tuple(g1dp)
        domain = tuple(domain)
        self.g1dp = g1dp
        r = np.meshgrid(*g1dp)
        r = [ri.reshape(ri.size) for ri in r]
        self.b = rhs_lambda(*r)
        if mesh.dim == 3:
            raise NotImplementedError
        else:
            invhx2 = 1./(mesh.d[0]**2)
            invhy2 = 1./(mesh.d[1]**2)

            d00 = np.ones(mesh.linear_size, dtype=np.int8)*(-2.*(invhx2 + invhy2))
            dx  = np.ones(mesh.linear_size, dtype=np.int8)*      invhx2
            dy  = np.ones(mesh.linear_size, dtype=np.int8)*               invhy2
            A = scipy.sparse.dia_array([d00,dx,dx,dy,dy], offsets=[0,-1,1,])
        d = mesh.d
        if bc_scale is None:
            A = scipy.sparse.linalg.LaplacianNd(domain, boundary_conditions='dirichlet')
        else:
            A = scipy.sparse.linalg.LaplacianNd(domain, boundary_conditions='dirichlet', dtype=np.float64)
        A = A.tosparse()

        A.data *=-1 # because we solve -u" = f
        self.bc_scale = bc_scale
        if bc_scale is None:
            A.data[0,-2] = 0
            A.data[1, 0] = 1
            A.data[1,-1] = 1
            A.data[2, 1] = 0
        else:
            A.data[1, 0] += bc_scale
            A.data[1,-1] += bc_scale
        self.A = A

        # print(self.A.toarray())
        self.u = np.zeros(M+2, dtype=np.float64)

    def assemble(self, bcs):
        """Assemble the matrix A and the column vector b."""
        h = self.mesh.d[0]

        # assemble b
        b = self.b
        b *= h**2

        # apply BCs
        A = self.A
        bc_left = bcs['left']
        if 'D' in bc_left:
            u0 = bc_left['D']
            if self.bc_scale is None:
                b[0] = u0
            else:
                b[0] = u0*self.bc_scale
        else:
            raise NotImplementedError

        bc_right = bcs['right']
        if 'D' in bc_right:
            u1 = bc_right['D']
            if self.bc_scale is None:
                b[-1] = u1
            else:
                b[-1] = u1 * self.bc_scale
        else:
            raise NotImplementedError
        print(b)

    def solve(self):
        self.u = scipy.sparse.linalg.spsolve(self.A, self.b)
        return self.u


def case_2D():
    # Use sympy to compute a rhs, given an analytical solution
    x, y = sympy.symbols("x y", real=True)
    sigma2 = 2.**2
    ue = sympy.exp(-0.5*(x**2 + y**2)/sigma2)  # a gaussian curve
    ue_lambda = sympy.lambdify([x,y], ue.simplify(), "numpy")
    fe = -(ue.diff(x, 2) + ue.diff(y, 2))
    fe_lambda = sympy.lambdify([x,y], fe.simplify(), "numpy")

    tbl = []
    M = 5
    d = 1.6
    for i in range(10):
        M *= 2
        d /= 2.
        mesh = LagrangeMesh(M, d, dim=2, reduced=False)

        pe1d = PoissonSolverNd(mesh, rhs_lambda=fe_lambda, bc_scale=1e20)
        bcs = {
            'left' : {'D': ue_lambda(pe1d.x[ 0])},
            'right': {'D': ue_lambda(pe1d.x[-1])}
        }
        pe1d.assemble(bcs)
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