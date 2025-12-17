See [feature request 68](https://github.com/IAA-nuclear/tantalus_full/issues/68)

# The Poisson equation
The Poisson equation reads:
[eq 1]
$$

	\mathbf{\Delta} V_c(\mathbf{r}) =- 4 \pi e^2 \rho_c(\mathbf{r}).$$
It relates a charge distribution $\rho_c(\vec{r})$ to its Coulomb potential $V_c(\vec{r})$. $\mathbf{\Delta}$ is the Laplacian $$
\mathbf{\Delta}=\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}+\frac{\partial^2}{\partial z^2}
$$A slightly more general version reads:
[eq 2]
$$
( A \mathbf{\Delta} + b ) F(\mathbf{r}) = f(\mathbf{r})
$$
Where $A$ and $b$ are scalars, $f(\mathbf{r})$ is a scalar function of $\mathbf{r}=(x,y,z)$ and is either symmetric or skew-symmetric wrt the three coordinate axes
The aim is to solve [eq 2] for $F(\mathbf{r})$ on a finite difference grid with the same characteristics (`N`,`d`,`reduce`) as the LagrangeMesh. 

## Boundary conditions
The behavior of the Poisson potential at $\mathbf{r}\rightarrow \infty$ is governed by :
[eq 3]
$$V(\mathbf{r}\rightarrow \infty) \sim \sum_{ \ell m} \frac{Q_{\ell m} Y^*_{\ell m}(r)}{r^\ell} $$
where $\ell \geq 0$ and $- \ell \leq m \leq \ell$ are integers, $Y_{\ell m}$ are spherical harmonics and $Q_{\ell m}$ are multipole moments of the charge distribution. [Eq 3] is used for applying Dirichlet boundary conditions at the outside of the mesh.  

The first term in the rhs is:
$$
e^2 \frac{Z}{r}
$$
$$
Z = \int d\mathbf{r} \rho_c(\mathbf{r})
$$
(For the time being we restrict to only the first term). 
The accuracy of the solver can be improved by increasing the domain, i.e. adding grid points (expensive) or adding higher order terms in [eq 3] (typically cheaper). 

If the mesh is reduced then we have at  $x=0$:
- a Dirichlet bc $V=0$ if $f$ is skew-symmetric wrt $x$  
- a Neumann bc $\frac{\partial V}{\partial x}=0$ if $f$ is symmetric wrt $x$
These bc pose a problem because this introduces a grid spacing of $d_x/2$ between $x=0$ and the first grid point, whereas the spacing between the other grid point is $d_x$. We could solve this by introducing intermediate grid points (expensive). Mocca solves this by not applying bc at $x=0$ and using the symmetry of $f$ in the stencil:
$$
\left. \Delta f(x) \right|_{x = 1/2 d_x} = \frac{f(3/2d_x) - 2f(1/2d_x) + f(-1/2d_x)}{2d_x} = \frac{f(3/2d_x) - 2f(1/2d_x) + s_x f(1/2d_x)}{2d_x}.
$$
where $s_x$ is $+1$ if $f$ is symmetric, and $-1$ if $f$ is skew-symmetric.
> [!Note] Conclusion
> So, actually the FD method uses Dirichlet BC

## Considerations
Obviously, we are not the only ones that need to solve Poisson equations. There are probably efficient and high-quality codes out there. We asked Google

> [!Tip] Question
> geef me wat softwares die een 3D poissonvergelijkeing kunnen oplossen

> [!Warning] Answer
> Er zijn verschillende softwarepakketten en programmeerbibliotheken beschikbaar die de 3D Poisson-vergelijking kunnen oplossen. De keuze hangt af van uw specifieke behoeften, zoals de vereiste nauwkeurigheid, de complexiteit van het domein en of u een commerciële of open-source oplossing prefereert.
>
Hier zijn enkele populaire opties:
>
Commerciële Softwarepakketten
>
Deze pakketten bieden doorgaans een gebruiksvriendelijke grafische interface, geavanceerde pre- en post-processing tools en uitgebreide ondersteuning. Ze zijn ontworpen voor complexe, grootschalige engineering- en fysicasimulaties.
>
>- **COMSOL Multiphysics**: Een krachtig eindige-elementenanalyse (EEM) softwarepakket waarmee u gekoppelde fysische verschijnselen en systemen van partiële differentiaalvergelijkingen (PDE's) kunt modelleren, inclusief de Poisson-vergelijking in 3D.
>- **[ANSYS](https://www.ansys.com/)**: Een toonaangevend commercieel FEA-pakket (Finite Element Analysis) dat bekend staat om zijn uitgebreide analysemogelijkheden in diverse sectoren.
>- **Abaqus**: Een ander prominent FEA-pakket, uitblinkend in complexe simulaties, vooral op het gebied van niet-lineaire analyse. 
>
Open-Source Bibliotheken en Omgevingen
>
>Voor gebruikers die meer flexibiliteit, aanpasbaarheid of kosteneffectieve oplossingen zoeken, zijn er diverse open-source tools. Deze vereisen vaak meer programmeerkennis.
>
>- **FEniCS Project**: Een populair open-source computerplatform voor het oplossen van PDE's met behulp van de eindige-elementenmethode. Het biedt interfaces in Python en C++ en stelt gebruikers in staat om snel wetenschappelijke modellen in efficiënte code om te zetten.
>- **MATLAB**: Hoewel MATLAB zelf een commercieel product is, biedt het een uitgebreide omgeving en ingebouwde functies (`solvepde`, `pdepe`) en de `Partial Differential Equation Toolbox` voor het numeriek oplossen van PDE's, waaronder de 3D Poisson-vergelijking. Veel scripts en voorbeelden zijn online beschikbaar.
>- **Python-bibliotheken (met NumPy, SciPy, etc.)**: Diverse Python-modules en -bibliotheken kunnen worden gebruikt.
 >   - **shenfun**: Een Python-module die spectrale methoden gebruikt om de 3D Poisson-vergelijking in tensorproductdomeinen op te lossen.
 >   - **DssPyLib**: Een open-source Python-software die de eindige-elementen- en randelementenmethoden gebruikt om elektrostatische problemen (die de Poisson-vergelijking betreffen) op te lossen.
>    - **Specifieke solvers op [GitHub](https://github.com/)**: Er zijn diverse gespecialiseerde, lichtgewicht solvers beschikbaar op platforms zoals GitHub, vaak gebaseerd op methoden als multigrid of Fast Fourier Transforms (FFT).
>- **Oomph-lib**: Een open-source, objectgeoriënteerde eindige-elementenbibliotheek geschreven in C++ met Python-bindings, die adaptieve oplossingen voor de 3D Poisson-vergelijking ondersteunt. 
>
>De keuze hangt dus af van uw budget en technische expertise. Voor snelle, gebruiksklare oplossingen zijn commerciële pakketten vaak efficiënter, terwijl voor academisch onderzoek of aangepaste oplossingen open-source bibliotheken en programmeeromgevingen zoals FEniCS of MATLAB uitstekende opties zijn.

Before making a choice, let us write down a wish list to rank candidates:
- we require open source code
- preferably with a Python front-end
- efficient
- not tremendously large
- intuitive expression of the problem to be solved, not just the Poisson equation (Wouter is interested in variants of it e.g. [eq 2])
- flexible handling of boundary conditions, Dirichlet, Neuman, and periodic BC.
- flexible translation of the problem domain and the BC.
- MPI aware.
- well documented, and maintained code
- 
[shenfun](https://github.com/spectralDNS/shenfun) seems to be well suited. I has an intelligently designed API. In addition it has some [demos](https://shenfun.readthedocs.io/en/latest/) featuring the Poisson equations. From the [Introduction](https://shenfun.readthedocs.io/en/latest/introduction.html):

"The spectral Galerkin method solves partial differential equations through a special form of the [method of weighted residuals](https://en.wikiversity.org/wiki/Introduction_to_finite_elements/Weighted_residual_methods) (WRM). As a Galerkin method it is very similar to the [finite element method](https://en.wikipedia.org/wiki/Finite_element_method) (FEM). The most distinguishable feature is that it uses global shape functions, where FEM uses local. This feature leads to highly accurate results with very few shape functions, but the downside is much less flexibility when it comes to computational domain than FEM."

Since we are only interested in rectangular domains, the lack of flexibility on the computational domain is not a problem.

- An additional element on our wish list is that we would prefer to use the grid points of the Lagrange mesh when using a FD method. In a LagrangeMesh the grid points are in the **middle** of the cells, not at the cell boundaries. This has consequences on how the boundary conditions are implemented.

E.g. in 1D with $M=10$, we have 10 cells over a domain $[-5d,5d]$ and grid points 
$$-4.5d, -3.5d, -2.5d, -1.5d, -0.5d, 0.5d, 1.5d, 2.5d, 3.5d, 4.5d$$
Or on a reduced grid 
$$0.5d, 1.5d, 2.5d, 3.5d, 4.5d$$
The boundary conditions, however need to be applied at $x=-5d$, and $x=5d$, and on $x=0, and $x=5d$ for a reduced grid. g which are not grid points.

Assuming a symmetric potential $u(-x)=u(x)$, then $u(-\frac{d}{2}) = u(\frac{d}{2})$. Consequently, 
$$u^"(\frac{d}{2}) = \frac{-1 u(-\frac{d}{2}) + 2u(\frac{d}{2}) -1u(\frac{3d}{2})}{h^2}
$$
$$
= \frac{-1 u(\frac{d}{2}) + 2u(\frac{d}{2}) -1u(\frac{3d}{2})}{h^2}
$$
$$
= \frac{(2-1)u(\frac{d}{2}) - u(\frac{3d}{2})}{h^2}
$$
 For skew-symmetric $u$ we have
$$
u^"(\frac{d}{2}) = \frac{(2+1)u(\frac{d}{2}) - u(\frac{3d}{2})}{h^2}
$$
So we must simply change the $A$ matrix, rather than adding an extra point  or applying an explicit boundary condition.
## Finite difference method in 3D
3 point stencil for the 2nd order derivative:
$$
f″(x) \approx \frac{f(x-h)−2f(x)+f(x+h)}{h^2}
$$
$$
= f″(x) \approx \frac{f_{i-1}−2f_{i}+f_{i+1}}{h^2}
$$
$$
[\frac{1}{h^2} \frac{-2}{h^2} \frac{1}{h^2} ]
$$
For the Laplacian in 3D we have:
$$
\Delta f(x_{ijk}) = (\frac{\partial^2 f}{\partial x^2}+\frac{\partial^2 f}{\partial y^2}+\frac{\partial^2 f}{\partial z^2})(x_{ijk})
$$$$
\approx \frac{f_{i-1,j,k}−2f_{i,j,k}+f_{i+1,j,k}}{h_x^2}
+ \frac{f_{i,j-1,k}−2f_{i,j,k}+f_{i,j+1,k}}{h_y^2}
+ \frac{f_{i,j,k-1}−2f_{i,j,k}+f_{i,j,k+1}}{h_z^2}
$$
$$
= \frac{f_{i-1,j,k}}{h_x^2} + \frac{f_{i,j-1,k}}{h_y^2} + \frac{f_{i,j,k-1}}{h_z^2}
-2(\frac{f_{i,j,k}}{h_x^2} + \frac{f_{i,j,k}}{h_y^2} + \frac{f_{i,j,k}}{h_z^2})
+\frac{f_{i+1,j,k}}{h_x^2} + \frac{f_{i,j+1,k}}{h_y^2} + \frac{f_{i,j,k+1}}{h_z^2}
$$
If $h_x=h_y=h_z=h$ this becomes:
$$
\approx \frac{f_{i-1,j,k}+f_{i,j-1,k}+f_{i,j,k-1}−6f_{i,j,k}+f_{i+1,j,k}+f_{i,j+1,k}+f_{i,j,k+1}}{h^2}
$$
in which case the $h^2$ can be moved to the right hand side and the matrix A is integer with diagonals of $1$ and $-6$ . Otherwise the diagonals are $1/h_x^2$, $1/h_y^2$, $1/h_z^2$, and $-6(1/h_x^2 + 1/h_y^2+ 1/h_z^2)$.

5 point stencil for the 2nd order derivative:
$$
f″(x)=\frac{−f(x-2h)+16f(x-h)−30f(x)+16f(x+h)−f(x+2h)}{12h^2}
$$
$$
\Delta f(x_{ijk}) = (\frac{\partial^2 f}{\partial x^2}+\frac{\partial^2 f}{\partial y^2}+\frac{\partial^2 f}{\partial z^2})(x_{ijk})
$$
$$\approx -\frac{1}{12}
(\frac{f_{i-2,j,k}}{h_x^2} + \frac{f_{i,j-2,k}}{h_y^2} + \frac{f_{i,j,k-2}}{h_z^2}
+\frac{f_{i+2,j,k}}{h_x^2} + \frac{f_{i,j+2,k}}{h_y^2} + \frac{f_{i,j,k+2}}{h_z^2})
$$
$$
+\frac{16}{12}(\frac{f_{i-1,j,k}}{h_x^2} + \frac{f_{i,j-1,k}}{h_y^2} + \frac{f_{i,j,k-1}}{h_z^2}
+\frac{f_{i+1,j,k}}{h_x^2} + \frac{f_{i,j+1,k}}{h_y^2} + \frac{f_{i,j,k+1}}{h_z^2})
$$$$
-\frac{30}{12}(\frac{1}{h_x^2} + \frac{1}{h_y^2} + \frac{1}{h_z^2})f_{i,j,k}
$$$$=
-\frac{1}{12}
(\frac{f_{i-2,j,k}}{h_x^2} + \frac{f_{i,j-2,k}}{h_y^2} + \frac{f_{i,j,k-2}}{h_z^2}
+\frac{f_{i+2,j,k}}{h_x^2} + \frac{f_{i,j+2,k}}{h_y^2} + \frac{f_{i,j,k+2}}{h_z^2})$$
$$
+\frac{4}{3}(\frac{f_{i-1,j,k}}{h_x^2} + \frac{f_{i,j-1,k}}{h_y^2} + \frac{f_{i,j,k-1}}{h_z^2}
+\frac{f_{i+1,j,k}}{h_x^2} + \frac{f_{i,j+1,k}}{h_y^2} + \frac{f_{i,j,k+1}}{h_z^2})
$$$$
-\frac{5}{2}(\frac{1}{h_x^2} + \frac{1}{h_y^2} + \frac{1}{h_z^2})f_{i,j,k}
$$Again, if $h_x=h_y=h_z=h$ all denominators are $12h^2$, which can be moved to the right hand side and the matrix A is integer with diagonals of $-1$, $16$, and $-30$. Otherwise the diagonals are $\frac{-1}{12h_x^2}$, $\frac{-1}{12h_y^2}$, $\frac{-1}{12h_z^2}$, $\frac{4}{3h_x^2}$, $\frac{4}{3h_y^2}$, $\frac{4}{3h_z^2}$ and $-\frac{5}{2}(\frac{1}{h_x^2} + \frac{1}{h_y^2} + \frac{1}{h_z^2})$.

## Construction of the Laplacian Matrix in more dimensions
The Kronecker product of two matrices $A$ ($N_x \times N_x$) and $B$ ($N_y \times N_y$) is defined as:
$$
A \otimes  B = 
\begin{bmatrix}
A_{11}B & A_{12}B & \dots  & A_{1N}B \\
A_{21}B & A_{22}B & \dots  & A_{2N}B \\
\vdots  & \vdots  & \ddots &\vdots   \\
A_{N1}B & A_{N2}B & \dots  & A_{NN}B 
\end{bmatrix}
$$
which is a $N_x N_y \times N_x N_y$ matrix.
Let $D_{2x}^{1D}$ be the derivative matrix ($N_x \times N_x$) for the second derivative wrt $x$, i.e. $\frac{\partial^2f}{\partial x^2} = D_{2x}^{1D}f$.  Then 
$$
D_{2x}^{2D} = D_{2x}^{1D} \otimes I_{N_y}
$$
$$
D_{2y}^{2D} = I_{N_x} \otimes D_{2y}^{1D}
$$
where $D_{2y}^{1D}$ has the same structure as $D_{2x}^{1D}$, but a different shape ($N_y \times N_y$, vs. $N_x \times N_x$).
The 2D Laplacian matrix ($N_x N_y \times N_x N_y$)  then is 
$$
L^{2D} 
= D_{2x }^{1D} \otimes I_{N_y}
+ I_{N_x}      \otimes D_{2y}^{1D}
$$
This formula is for lexicographic order of the indices of the right hand side (i.e. the rightmost index varies the fastest). 
The first term here puts $N_x$ copies of $D_{2y}^{1D}$ along the diagonal, su$ch that the structure of $D_{2y}^{1D}$ is copied onto $L_{2D}$, almost, because at the corners the subdiagonals are interrupted by zeros ($\textcolor{red}{0}$):
$$
\begin{bmatrix}
\ddots & \ddots & \ddots & \ddots & \ddots & 0 & 0 & 0 & 0 & \dots \\
\dots & 0 & 1 & -2 & 1 & 0 & 0 & 0 & 0 & \dots \\
\dots & 0 & 0 & 1 & -2 & \textcolor{red}{0} & 0 & 0 & 0 & \dots \\
\dots & 0 & 0 & 0 & \textcolor{red}{0} & -2 & 1 & 0 & 0 & \dots \\
\dots & 0 & 0 & 0 & 0 & 1 & -2 & 1 & 0 & \dots \\
\dots & 0  & 0  & 0  & 0 & \ddots & \ddots & \ddots & \ddots & \ddots
\end{bmatrix}
$$
Here, the $-2$ above $\textcolor{red}{0}$ is the lower right corner of a D

As demonstrated  in marimo notebook [marimo-kron.py](MOCCaPy/design_docs/marimo-kron.py) the Laplacian matrix in 2D can be constructed from the 1D Laplacian matrix using sums of Kronecker products. For Fortran ordering of the right hand side the formulas read:

In 2D the Kronecker formula reads:
For lexicographic order

$$
L^{2D} = D_{xx}^{1D} \otimes I_{N_y} + I_{N_x} \otimes D_{yy}^{1D}
$$

For Fortran order:
$$
L^{2D} = I_{N_y} \otimes D_{xx}^{1D} + D_{yy}^{1D} \otimes I_{N_x}
$$

In 3D the Kronecker formula reads:
For lexicographic order:
$$
L^{3D} = D_{xx}^{1D} \otimes I_{N_y} \otimes I_{N_z} + I_{N_x} \otimes D_{yy}^{1D} \otimes I_{N_z} + I_{N_x} \otimes I_{N_y} \otimes D_{zz}^{1D}
$$

For Fortran order:

$$
L^{3D} = I_{N_z} \otimes I_{N_y} \otimes D_{xx}^{1D} + I_{N_z} \otimes D_{yy}^{1D} \otimes I_{N_x} + D_{zz}^{1D} \otimes I_{N_y} \otimes I_{N_x}
$$
There exist implementations for the Kronecker product in Numpy ([`np.kron`](https://numpy.org/doc/2.1/reference/generated/numpy.kron.html#numpy.kron)) and SciPy ([`scipy.sparse.kron`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.kron.html)) 

## Reduced axes

3-point stencil
$$
\left. \Delta f(x) \right|_{x = \frac{1}{2} dx} 
= \frac{
f(-\frac{1}{2}dx)  - 2f(\frac{1}{2}dx) + f(\frac{3}{2}dx)
}{dx^2} 
= \frac{
s_x f(\frac{1}{2}dx) - 2f(\frac{1}{2}dx)  + f(\frac{3}{2}dx) 
}{dx^2}.
$$
$$
= \frac{
(s_x - 2)f(\frac{1}{2}dx)  + f(\frac{3}{2}dx) 
}{dx^2}.
$$
that is, $s_x$ is added to the midpoint's matrix entry. This is on the main diagonal.
5-point stencil
$$
\left. \Delta f(x) \right|_{x = \frac{1}{2} dx} 
= \frac{
- f(-\frac{3}{2}dx) + 16f(-\frac{1}{2}dx) - 30f(\frac{1}{2}dx) + 16f(\frac{3}{2}dx) - f(\frac{5}{2}dx)
}{12dx^2} 
$$
$$
= \frac{
- s_xf(\frac{3}{2}dx) + 16s_xf(\frac{1}{2}dx) - 30f(\frac{1}{2}dx) + 16f(\frac{3}{2}dx) - f(\frac{5}{2}dx)
}{12dx^2}
$$
$$
= \frac{
(16s_x-30)f(\frac{1}{2}dx) + (16-s_x)f(\frac{3}{2}dx) - f(\frac{5}{2}dx)
}{12dx^2} 
$$
that is the main diagonal gets an extra $16s_x$ and the +1 diagonal an extra $-s_x$.
Also the second point is influenced:
$$
\left. \Delta f(x) \right|_{x = \frac{3}{2} dx} 
= \frac{
- f(-\frac{1}{2}dx) + 16f(\frac{1}{2}dx) - 30f(\frac{3}{2}dx) + 16f(\frac{5}{2}dx) - f(\frac{7}{2}dx)
}{12dx^2} 
$$
$$
= \frac{
- s_xf(\frac{1}{2}dx) + 16s_xf(\frac{1}{2}dx) - 30f(\frac{3}{2}dx) + 16f(\frac{5}{2}dx) - f(\frac{7}{2}dx)
}{12dx^2}
$$
$$
= \frac{
(16-s_x)f(\frac{1}{2}dx) -30f(\frac{3}{2}dx) + 16f(\frac{5}{2}dx) - f(\frac{7}{2}dx)
}{12dx^2} 
$$
and, consequently, $L$ remains a symmetric matrix, since the -1 diagonal get an extra $-s_x$.
### issues
The 3-point stencil in 1D leads to the linear system:
```
solution: u(r)=gauss, sigma=2, symmetry=1
mesh    : dim=1, M=12, h=0.8, reduced=False
solver  : bc_scale=None, stencil=3, method='direct'
linear system : 12x12
  i\j  point             0    1    2    3    4    5    6    7    8    9    10    11           f
-----  --------------  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ----  ----  ----------
    0  b[0][-4.40000]    1    0    0    0    0    0    0    0    0    0     0     0   0.0889216
    1  i[1][-3.60000]    1   -2    1    0    0    0    0    0    0    0     0     0   0.0709269
    2  i[2][-2.80000]    0    1   -2    1    0    0    0    0    0    0     0     0   0.0576478
    3  i[3][-2.00000]    0    0    1   -2    1    0    0    0    0    0     0     0   0
    4  i[4][-1.20000]    0    0    0    1   -2    1    0    0    0    0     0     0  -0.0855317
    5  i[5][-0.40000]    0    0    0    0    1   -2    1    0    0    0     0     0  -0.150559
    6  i[6][0.40000]     0    0    0    0    0    1   -2    1    0    0     0     0  -0.150559
    7  i[7][1.20000]     0    0    0    0    0    0    1   -2    1    0     0     0  -0.0855317
    8  i[8][2.00000]     0    0    0    0    0    0    0    1   -2    1     0     0   0
    9  i[9][2.80000]     0    0    0    0    0    0    0    0    1   -2     1     0   0.0576478
   10  i[10][3.60000]    0    0    0    0    0    0    0    0    0    1    -2     1   0.0709269
   11  b[11][4.40000]    0    0    0    0    0    0    0    0    0    0     0     1   0.0889216
bc_scale=None M=12 : rmse=np.float64(0.008798406090593368) mean_diff=0.005923604338136682 max_diff=0.01741974052755879 0.00058s 0.02354s
```
the 5-point stencil gives
```
solution: u(r)=gauss, sigma=2, symmetry=1
mesh    : dim=1, M=12, h=0.8, reduced=False
solver  : bc_scale=None, stencil=5, method='direct'
linear system : 12x12
  i\j  point             0    1    2    3    4    5    6    7    8    9    10    11           f
-----  --------------  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ----  ----  ----------
    0  b[0][-4.40000]    1    0    0    0    0    0    0    0    0    0     0     0   0.0889216
    1  i[1][-3.60000]   16  -30   16   -1    0    0    0    0    0    0     0     0   0.851123
    2  i[2][-2.80000]   -1   16  -30   16   -1    0    0    0    0    0     0     0   0.691773
    3  i[3][-2.00000]    0   -1   16  -30   16   -1    0    0    0    0     0     0   0
    4  i[4][-1.20000]    0    0   -1   16  -30   16   -1    0    0    0     0     0  -1.02638
    5  i[5][-0.40000]    0    0    0   -1   16  -30   16   -1    0    0     0     0  -1.8067
    6  i[6][0.40000]     0    0    0    0   -1   16  -30   16   -1    0     0     0  -1.8067
    7  i[7][1.20000]     0    0    0    0    0   -1   16  -30   16   -1     0     0  -1.02638
    8  i[8][2.00000]     0    0    0    0    0    0   -1   16  -30   16    -1     0   0
    9  i[9][2.80000]     0    0    0    0    0    0    0   -1   16  -30    16    -1   0.691773
   10  i[10][3.60000]    0    0    0    0    0    0    0    0   -1   16   -30    16   0.851123
   11  b[11][4.40000]    0    0    0    0    0    0    0    0    0    0     0     1   0.0889216
bc_scale=None M=12 : rmse=np.float64(0.0024529569767277833) mean_diff=0.0022026120883708055 max_diff=0.0034458770525791493 0.00065s 0.01756s
```
We noted above that the 5-point stencil does not converge. The problem is that for point 1 and 23, the boundary points first neighbors, the stencil lacks the $-1$ contribution from the virtual points -1 and 24. Either we should adapt the stencil for these point to account for this symmetry. Or we could provide the extra points, in which we get two extra points, and thus two extra unknowns and two extra equations. This can be accounted for by two extra Dirichlet boundary conditions. This is mathematically sound as their values are obtained from the fact that we know the long range behavior of the solution. Instead of providing extra points we can apply Dirichlet boundary conditions for points 0, 1, 10 and 11 and solve for the unknowns 2-9.
This is simpler (especially in 2D and 3D) than embedding the grid in a larger grid.

This means that for the 5-point stencil we need to collect boundary points and next-to-boundary points. 

The MOCCa solution was to embed the Lagrange grid into a grid with 1 (3-point stencil) or 2 (5-point stencil) points extra in every direction and computing L for the interior points only (that I am not sure of but it seems plausible) and only solve for the unknowns corresponding to the interior nodes (i.e. the LagrangeMesh nodes). This approach replaces the `LagrangeMesh.collect_boundary_nodes(stencil)` by a `LagrangeMesh.add_boundary_nodes(stencil)`

Embedding the `LagrangeMesh` in a larger mesh for solving the Poisson can be done by creating an extra mesh or with a restricted view on the larger mesh. E.g in 1D the Poisson mesh with 20 points for a 5-point stencil would be a numpy array `poisson_mesh` with shape `(20,)` and the corresponding LagrangeMesh would have the points `poisson_mesh[2:-2]`. The difficulty is that the linear system matrix $A^{16 \times 16}$  is no longer the Laplacian matrix $L^{20\times 20}$ over the Poisson mesh. Embedding, however, allows $L$ to be expressed as a linear operator which applies the stencil at every interior point, computing $Au$ (thus it is essentially a matrix-free approach), while reaching into the embedding region when necessary. The system can then be solved iteratively, for the interior points only.  That avoids special treatments of boundary and symmetry boundary conditions all together (except for copying the symmetry points at each iteration). To that end, SciPy provides `LinearOperator` class, whose constructor simply requires a method that computes $Au$ given $u$. The difficulty for this approach is that (except for 1D) is impossible to embed the grid and keep $u$ a contiguous vector. Either we deal with
- a contiguous $u$ and store the boundary regions separately (non-contigously), 
- or with an embedded grid and a non-contiguous u.
This requires the following steps:
0. set the boundary values on the boundary region of the embedded grid,
1. put $u$ on the embedded grid (when using an embedded grid this is a non-contiguous copy operation),  
2. copy the symmetry points in $u$ to the symmetry boundary regions with the correct sign (as $u$  is updated during iterations), 
3. compute $Au$ (when using the embedded grid, this involves applying the stencil to the interior points, which is rather straightforward. If not, this involves applying the stencil on the interior points for which the stencil does not extend beyond the grid, and a whole series of corner cases where the stencil extends into the boundary and symmetry boundary regions) 
4. update (inside iterative solver),
5. back to 1. and repeat until convergence (inside iterative solver).
Embedding the grid seems the less error-prone approach. 

