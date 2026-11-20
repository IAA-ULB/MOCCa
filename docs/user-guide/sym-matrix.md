# The symmetry-block matrix type (`sym_matrix`)

This page is the developer-level reference for `sym_matrix`, the type the
`src_orig/block_matrix_module.f90` module provides for single-particle
matrix elements between the stored Hartree-Fock wavefunctions. It is aimed
at code that computes expectation values, one-body operator matrix
elements or basis transforms — the multipole generator `Qlm_spme_matrix`
in `moments.f90`, the angular-momentum generators of
`angular_momentum.f90` and the gradient generators of
`linear_momentum.f90` are the production examples. End users running MOCCa
never see the type; it appears only in developer output such as the
`print_dense` debug dump.

The authoritative record of conventions, measured facts and design
amendments is the design summary at the repository root
(`mocca-symmetry-block-matrix-design-summary.md`); this page restates
nothing from it — when the two disagree, code and design summary win.

## The symmetry model

The HF basis states carry three discrete good quantum numbers; the sector
word `s = 0..7` encodes them bit by bit:

| bit | quantum number | 0 | 1 |
|-----|----------------|---|---|
| 0 | z-signature $q(R_z)$ | $+i$ | $-i$ |
| 1 | parity $\hat P$ | even | odd |
| 2 | isospin $t_z$ | neutron | proton |

MOCCa's `HFBlocks(B)` bookkeeping uses the same ordering: block $B$
counted the states of sector $s = B - 1$, so block 1 is
neutron/even/$+i$, block 8 is proton/odd/$-i$, and so on.

A one-body operator carries a symmetry word `op` with the same encoding: a
bit is set if the operator flips that quantum number. Selection rules are a
single XOR — an operator of word `op` couples sector $s$ to sector
$\mathrm{IEOR}(s,\mathrm{op})$ and to nothing else. A matrix therefore has
at most 8 non-vanishing dense blocks, and that is exactly what `sym_matrix`
stores: one contiguous rank-2 array per surviving sector pair, at most one
side per pair.

Not every build labels its basis with all three quantum numbers: in a
parity-broken build, for instance, parity is not a good quantum number, the
states of both parities are merged into one storage block, and a
parity-odd operator acts *within* the blocks instead of pairing into the
empty ones. `init` therefore takes the build's good-quantum-number mask
`goodqn` (in MOCCa produced by Hephaestos via `$GOODQN`; default 15 = all
good) and pairs sector $s$ with
$\mathrm{IEOR}\big(s,\,\mathrm{AND}(\text{op},\text{goodqn})\big)$.

## The class system and the storage kind

Every operator the type handles is *classifiable*: it is Hermitian or
anti-Hermitian, and it is even or odd under the y-time-simplex
$ST_y = P \cdot T \cdot R_y$. Since the stored wavefunctions are $ST_y$
invariants, an $ST_y$-even operator has purely **real** matrix elements
(`STY_EVEN`) and an $ST_y$-odd one purely imaginary ones (`STY_ODD`). Both
physical classes are **mandatory** `init` inputs; the storage kind is
derived from them:

$$k = \text{herm} \times \text{sty}$$

Only one half of the (real, imaginary) pair is stored, and the partner
blocks are reconstructed on demand through

$$A(s',s) = k\; A(s,s')^{T}.$$

In a build that conserves the time simplex (`goodqn` bit 3 set — true for
every configuration in the tree), the omitted half is identically zero and
the stored matrix *is* the complete operator.

The TR-folded composites ($J_\mu T$, $\nabla_\mu T$ — the objects whose
ket is the time-reversed state, with MOCCa's fermionic time reversal
$T^2 = -1$) are a special case: their class tags do not follow from the
linear-operator algebra. They are pinned to the kernel-transpose and
mirror oracles of the respective unit arms, which measure the mirror
character directly at every run. **Never re-derive a composite kind
sign**: the unit arms are the arbiter (see the design summary, sections
14.2b and 15.2).

## The type-bound surface

| Procedure | Meaning |
|-----------|---------|
| `init(HFBlocks, symmetry, herm, sty [, goodqn])` | Build the block list; zero-fill. Takes the **global** count array. Idempotent. |
| `fill_blocks(fill, params)` | Populate all stored blocks through a `block_filler` callback, once per block, overwriting the whole array handed out as `intent(out)`. |
| `sector_of(i)` | Global index to sector word; -1 out of range. |
| `dense()` | Partner-reconstructing dense $n\times n$ expansion (freshly allocated) — the shape behind `Qlm_spme_fm`. |
| `print_dense([lun][,title])` | Human-readable dump with sector boundaries; debug tool for small matrices only. |
| `rotate(R)` | Express the same operator in a rotated basis, in place. |

Population is exclusively through the frozen `block_filler` abstract
interface — extend the unlimited-polymorphic `params` object, never the
callback signature. The block's sector words (`srow`, `scol`) are handed
to the callback so producers can apply phase conventions.

## Basis rotations

`sym_rotation` holds one real orthogonal block per populated sector —
the shape MOCCa's own producers make (the HFB canonical basis
diagonalises $\rho$ per MOCCa block, so `cantransfo` is exactly
sector-diagonal). `rotate(R)` transforms every stored block as
$A \leftarrow U_{s_\text{row}}^T A\, U_{s_\text{col}}$ and preserves every
metadata field. Ordinary (real) rotations only: complex unitaries would
break the $ST_y$-invariance of the stored basis, and a rotation that
mixes sector partners would destroy the good quantum numbers the whole
type rests on. Both are refused loudly by the structure guards.

## The generator modules

The three production generators all follow the same pattern: resolve the
operator's symmetry word and physical classes (a small class table each,
arbitered by its unit arm), `init` on `HFBlocks_global` with the build's
`goodqn`, then `fill_blocks` with the legacy kernels — the arithmetic of
the legacy loops, verbatim, so those loops remain the A/B oracle of the
consumer migration.

| Module | Objects |
|--------|---------|
| `moments.f90` — `Qlm_spme_matrix` | Electric multipoles $r^l \mathrm{Re}/\mathrm{Im}\, Y_{lm}$; word $\mathrm{mod}(m,2) + 2\,\mathrm{mod}(l,2)$; `Qlm_spme_fm` returns the dense form |
| `angular_momentum.f90` — `Jmu_spme_matrix`, `Jt_spme_matrix` | Honest $J_x, J_y, J_z$ (words 1/1/0) and the folded $J_{x,y}T$ pair (word 0) |
| `linear_momentum.f90` — `Nabla_spme_matrix`, `NablaT_spme_matrix`, `Nabla_half_fm` | Honest $\nabla_\mu$ (words 3/3/2, all anti-Hermitian) and the folded pair; `Nabla_half_fm` is the dense form used by the centre-of-mass correction |

For the gradient family, `get_nabla_block` implements the *implicit
complete matrix* contract: any block of the doubled basis (stored states
plus their time-reversal images) is retrievable on demand — stored from
the generators, image-image blocks by the tau sign flip, cross blocks by
the honest kernels on materialised images, structural zeros elsewhere.
Only `$TR`-marked call sites may request blocks involving images.

## Testing entry points

The module's own net is standalone — it needs no MOCCa executable:

```
cd testing/unit && bash block_matrix_unit_tests.sh
```

It covers construction invariants, masked pairing, the dense expansion,
the rotation round trip and the structure guards that abort
(`block_matrix_refusal_tests.f90`) — the fastest gate in CI. In the
full executable, the arms of the generator modules (plus the moments
suite) run when `unit_test = .true.` is set in the `&mesh` namelist;
their oracle printouts are the authority the class tags in the modules
have to agree with.
