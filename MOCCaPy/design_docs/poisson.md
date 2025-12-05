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
[shenfun](https://github.com/spectralDNS/shenfun) seems to be well suited. I has an intelligently designed API. In addition it has some [demos](https://shenfun.readthedocs.io/en/latest/) featuring the Poisson equations. From the [Introduction](https://shenfun.readthedocs.io/en/latest/introduction.html):

"The spectral Galerkin method solves partial differential equations through a special form of the [method of weighted residuals](https://en.wikiversity.org/wiki/Introduction_to_finite_elements/Weighted_residual_methods) (WRM). As a Galerkin method it is very similar to the [finite element method](https://en.wikipedia.org/wiki/Finite_element_method) (FEM). The most distinguishable feature is that it uses global shape functions, where FEM uses local. This feature leads to highly accurate results with very few shape functions, but the downside is much less flexibility when it comes to computational domain than FEM."

Since we are only interested in rectangular domains, the lack of flexibility on the computational domain is not a problem.