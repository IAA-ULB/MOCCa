Work on [issue 42 proposal for next Tasks](https://github.com/IAA-nuclear/tantalus_full/issues/42#issuecomment-3543958442):

# Het plan van Wouter
Ik denk dat de beste strategie is om het "ruimtelijke" deel van de spwfs te manipuleren. M.a.w. we laten alle verschillen tussen HF/BCS/Bogoliubov achterwege, want dat betreft eigenlijk alleen maar het "pairing subprobleem" d.w.z. "gegeven een bepaalde set golffuncties op het mesh, hoe bouw ik een mean-field state"? Voor een discussie over de "subproblemen" verwijs ik je graag (opnieuw?) naar [mijn 2019 paper, Eur. Phys. J. A (2019) 55:93](../../literature/Ryssens_et_al_2019_Heavy_ball_dynamics_and_potential_preconditioning.pdf); met name het schema in Fig. 1

Dit betekent om eerst te focussen op "hoe evolueer ik mijn golffuncties op het mesh"? Dat is uiteindelijk hetgene dat tijd kost in kernberekeningen en wat ik denk dat veel beter kan met technieken uit kwantumchemie. We kunnen bovendien nog versimpelen door niet de single-particle hamiltoniaan van de functionaal te gebruiken maar een simpele analytische potentiaal: een Woods-Saxon, i.e. https://en.wikipedia.org/wiki/Woods%E2%80%93Saxon_potential. 

Hieronder een opsomming van de taken, telkens gevolgd door de uitwerking ervan.
## Task 1.  
> [!Tip] Task 1.
> Schrijf een "actionofh" routine die $h \psi$ voor alle (of een individuele) spwf (zie `subroutine actionofh` in [functional.f90](../../../src/functional.f90)) met 
$$h = \hbar^2/2m \Delta + V(r)$$
en $V(r)$ = de Woods-Saxon vorm. De waardes van $\hbar^2/2m$ zijn deel van een .param file.

That is, compute:
$$h\psi$$
where $\psi$ is represented by a MOCCaPy `HFPsi` object.
### The Hamiltonian
$$h = \hbar^2/2m \Delta + V(r)$$
where $\hbar^2/2m=20.73553000$, $\Delta$ represents the Laplacian, and $V(r)$ is the Woods-Saxon potential above. 
### The [Woods-Saxon potential](https://en.wikipedia.org/wiki/Woods–Saxon_potential)
$$V(r)=-\frac{V_0}{1+\exp(\frac{r-R}{a})}$$
with $V_0 = 50$ (MeV), $a=0.5$ (fm), and 
$$R=r_0A^{\frac{1}{3}}$$
$A$ being the mass number (total number of nucleons), $r_0=1.25$ fm
>[!Note]
>MOCCa expresses quantities in MeV and fm (femtometer).
### A bit of notation: correspondence between `HFPsi[:,:,i]` and $|\psi_i\rangle$
$|\psi_i\rangle$ is a 2-component spinor (whatever that may mean)
$$
|\psi_i\rangle
=\psi_i(\vec{r},\sigma)
=\begin{bmatrix}
\psi_i(\vec{r},\sigma=+) \\
\psi_i(\vec{r},\sigma=-)
\end{bmatrix}
$$
$$
=\begin{bmatrix}
\psi_i(\vec{r}, 0) + \mathrm{i}\,\psi_i(\vec{r}, 1) \\
\psi_i(\vec{r}, 2) + \mathrm{i}\,\psi_i(\vec{r}, 3)
\end{bmatrix}
$$
Here, the numbers $0..3$ refer to the four components of `HFPsi[:,0:4,i]`:
$$
|\psi_i\rangle
=\begin{bmatrix}
\mathrm{HFPsi[:,0,}i\mathrm{]} + \mathrm{i}\mathrm{HFPsi[:,1,}i\mathrm{]} \\
\mathrm{HFPsi[:,2,}i\mathrm{]} + \mathrm{i}\mathrm{HFPsi[:,3,}i\mathrm{]}
\end{bmatrix}$$
Applying an operator $O$ to `HFPsi` means applying it to all `4 x n_spwf` wavefunctions and yields a data structure of the same shape as `HFPsi`. 

When computing the expectation value e.g. of the overlap:
$$\langle \psi_i | \psi_j \rangle=\int d^3 r \sum_{\sigma = +-} \psi^*_j (r,\sigma) \psi_i(r,\sigma)$$

$$=\int d^3 r \sum_{\sigma = +-} \mathrm{Re} [\psi^*_j(r,\sigma) \psi_i(r,\sigma)]$$

$$=\int d^3 r \sum_{\sigma = +-} \left ( \mathrm{Re} [\psi_j(r,\sigma)] \mathrm{Re} [\psi_i(r,\sigma)] + \mathrm{Im} [\psi_j(r,\sigma)] \mathrm{Im} [\psi_i(r,\sigma)]\right)$$
 Above we used the symmetry assumption that $\langle \psi_i | \psi_j \rangle$ does not produce an imaginary component. Finally:
 $$\langle \psi_i | \psi_j \rangle=\sum_{r_{ijk}}
\begin{bmatrix}
+\psi_j(\vec{r}_{ijk},0)\psi_i(\vec{r}_{ijk},0) \\
+\psi_j(\vec{r}_{ijk},1)\psi_i(\vec{r}_{ijk},1) \\
+\psi_j(\vec{r}_{ijk},2)\psi_i(\vec{r}_{ijk},2) \\
+\psi_j(\vec{r}_{ijk},3)\psi_i(\vec{r}_{ijk},3) \\
\end{bmatrix}dv
$$
(the big brackets do not represent a matrix but a sum split over 4 lines). Above, $dv$ is the volume of a grid cell and incorporates a symmetry factor 2 for each reduced axis. The matri $\langle \psi_i | \psi_j \rangle$ is in fact the overlap matrix. In the case of an arbitry operator $O$ we get the operator's matrix representation:
$$\langle\psi_i|O|\psi_j\rangle=\sum_{r_{ijk}}
\begin{bmatrix}
+\psi_j(\vec{r}_{ijk},0)O\psi_i(\vec{r}_{ijk},0) \\
+\psi_j(\vec{r}_{ijk},1)O\psi_i(\vec{r}_{ijk},1) \\
+\psi_j(\vec{r}_{ijk},2)O\psi_i(\vec{r}_{ijk},2) \\
+\psi_j(\vec{r}_{ijk},3)O\psi_i(\vec{r}_{ijk},3) \\
\end{bmatrix}dv$$
### Block structure of `HFPsi`
`HFPsi` is a large array with shape `(mesh.linear_size, 4,n_total_wf)`. The single particle wave function are ordered in 8 blocks (originating roughly because we can exploit two spatial symmetries and one non-spatial symmetrie, the proton-neutron symmetry):
- block\[0]: neutrons with positive parity and signature $+\mathrm{i}$
- block\[1]: neutrons with positive parity and signature $-\mathrm{i}$
- block\[2]: neutrons with negative parity and signature $+\mathrm{i}$
- block\[3]: neutrons with negative parity and signature $-\mathrm{i}$
- block\[4]: protons with positive parity and signature $+\mathrm{i}$
- block\[5]: protons with positive parity and signature $-\mathrm{i}$
- block\[6]: protons with negative parity and signature $+\mathrm{i}$
- block\[7]: protons with negative parity and signature $-\mathrm{i}$
This is implemented in `hfblocks[8]` containing the number of single particle wave functions in each block, and `hfblocksrange[8]` containing a tuple giving the range of single particle wave function in the block, e.g. `HFPSi[:, :, hfblockrange[4][0]:hfblockrange[4][1]]` is the 4-th block with protons with positive parity and signature $+\mathrm{i}$. 
#### implementation
`Class Operator` (`mocca/mean_field/states/operator.py`) serves as a base class for operators and provides a standardized way to implement the action of an arbitrary operator and the computation of its matrix representation. Derived classes for the overlap operator (`Overlap`) and a Hamiltionian with a Woods-Saxon (`HamiltionianWoodsSaxon`) potential are provided too. In addition, a `KineticEnergyOperator` is provided from which Hamiltonian operators can derive (e.g. `HamiltionianWoodsSaxon2`).
## Task 2.
>[!Tip] Task 2.
>Schrijf een routine die een "imaginary time evolution" stap doet; Eq. 63 in  [Ryssens et al 2019 EPJA 55:93](../../literature/Ryssens_et_al_2019_Heavy_ball_dynamics_and_potential_preconditioning.pdf). 


## Task 3
>[!Tip] Task 3.
Schrijf een orthornomalisatieroutine die een set golffuncties neemt en die met Gram-Schmidt orthonormaliseert  (zie subroutine GramSchmidt in [wavefunctions.f90](../../../src/wavefunctions.f90)); het is belangrijk dat deze in energie-volgorde gebeurt.
## Task 4
>[!Tip] Task 4.
Schrijf manieren om
$$\langle\psi_i|h|\psi_i\rangle$$
en 
$$\langle\psi_i|h^2|\psi_i\rangle-\langle\psi_i|h|\psi_i\rangle^2$$
uit te rekenen. Vergeet niet dat $h$ hermitisch is, dus voor dat laatste kan je $h|\psi_j\rangle$ en $h|\psi_i\rangle$ veilig sandwichen.
## Task 5
>[!Tip] Task 5.
>Combineer dan alles in het algoritme beschreven in sectie 3.3 van [Ryssens et al 2019 EPJA 55:93](../../literature/Ryssens_et_al_2019_Heavy_ball_dynamics_and_potential_preconditioning.pdf) en verifieer dat je toestandjes convergeren naar eigentoestanden van h, bvb doordat $$\langle\psi_i|h^2|\psi_i\rangle-\langle\psi_i|h|\psi_i\rangle^2$$ heel klein wordt. 
