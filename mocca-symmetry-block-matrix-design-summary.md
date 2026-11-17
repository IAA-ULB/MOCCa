# Symmetry-Block Matrix Type for MOCCa: Design Summary

This document summarizes a design for a dense, block-structured matrix type to be
implemented in the MOCCa repository (github.com/IAA-ULB/MOCCa). It is written so
that an agent (VIBE CLI) or a collaborator can start implementation without the
preceding discussion.

## 1. Purpose and scope

MOCCa calculations are done in a basis of single-particle wavefunctions (spwfs)
that are symmetry-adapted and stored in a fixed order. Operators acting on this
basis respect three discrete symmetries, so an operator matrix decomposes into a
small number of dense blocks. The goal is a derived type that:

- stores only the nonzero blocks dictated by selection rules,
- keeps all block data contiguous and BLAS-friendly,
- performs mat-vec and mat-mat products at dense-BLAS speed,
- supports Hermitian and anti-Hermitian, real and complex operators,
- has zero quantum-number logic at runtime (all resolved at construction).

Constraints agreed with Wouter:

- Basis states are ordered; a simple integer array `HFblocks(0:7)` gives the
  number of states in each symmetry sector, in storage order.
- Isospin is limited to `T_z = ±1/2` (proton / neutron).
- Blocks are dense; matrices are constructed explicitly, not only via mat-vec.
- Total matrix sizes are typically around 1000 x 1000; single-node performance
  matters most. MPI distribution is out of scope for now.

## 2. Symmetry model

The three symmetries are parity (P), z-signature (S), and isospin (T_z). Each
acts on the basis as a Z2 grading:

- parity: even / odd
- z-signature: +1 / -1
- isospin: proton / neutron

Encode a symmetry sector as a 3-bit integer `s = 0..7`. Bit assignment
agreed with Wouter (2026-11-16), chosen so that a sector word coincides
with MOCCa's symmetry-block index minus one:

- bit 0: z-signature (Rz): 0 = q(Rz) = +i, 1 = q(Rz) = -i. Rz is linear
  and antihermitian, so its eigenvalues are +/- i; the bit carries their
  sign (time reversal flips it).
- bit 1: parity (0 = even, 1 = odd)
- bit 2: isospin (0 = neutron, 1 = proton)

An operator carries a symmetry word `op` (3 bits) with the same encoding:
bit 0 = whether the operator flips the Rz q-value sign, bit 1 = the parity
of the operator (odd operators carry 1), bit 2 = isoscalar
(bit clear, preserves proton/neutron) vs isovector-flipping in T_z (bit set).

Selection rule: an operator with symmetry word `op` couples sector `s` to
exactly one sector:

```
s' = IEOR(s, op)
```

Because XOR is an involution, each sector has exactly one partner and an
operator matrix has at most 8 nonzero blocks (fewer when sectors are empty).
This is the only place selection rules appear; everything downstream sees a
flat list of at most 8 dense blocks.

### 2a. Bijection with MOCCa's HFBlocks bookkeeping (confirmed)

Confirmed against the sources: wavefunctions.f90 (block bookkeeping header,
`iniwavefunctions`), gennilsson.f90 and nil8.f90 (all Nilsson seeds have
q(Rz) = +i; `kparz` is the parity), transform.f90 (block pairing/merging when
reading more symmetric solutions) and src_heph/heph_symmetries.py
(`multiply_quantum_numbers`: per-isospin block slots 1,2 carry the positive
parity-q, slots 3,4 the negative; even slots carry q = -i).

MOCCa's `HFBlocks(B)` counts the spwfs of sector `s = B - 1`:

| MOCCa block B | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| sector s = B-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
| isospin | n | n | n | n | p | p | p | p |
| parity | + | + | - | - | + | + | - | - |
| q(Rz) | +i | -i | +i | -i | +i | -i | +i | -i |

Occupation patterns of the basis sectors:

- EV8 symmetries ('Rz,T,P,STy'): blocks 1, 3, 5, 7 -> sectors 0, 2, 4, 6.
- Time-reversal broken: blocks 2, 4, 6, 8 additionally populated ->
  sectors 1, 3, 5, 7.
- Parity broken ('Rz,T,STy' / 'Rz,STy'): all states in blocks 1, 2, 5, 6
  -> sectors 0, 1, 4, 5; blocks 3, 4, 7, 8 remain empty throughout the
  run and the parity bit of a basis sector is vacuous.

Convention note for antiunitary operations: if any of the three symmetry
operations is antiunitary (involves complex conjugation) in the convention
used, the imaginary part of a complex operator carries that bit complemented
(see section 8). The type supports this by letting the real and imaginary
parts carry independent symmetry words; it must be decided once, in the
constructor, based on MOCCa's existing symmetry conventions for q-values.

### 2b. Good quantum numbers of the stored basis (agreed 2026-11-17)

The operator word of section 2 describes the *physical* character of an
operator. The XOR pairing additionally presumes that the sector words
describe definite quantum numbers of the stored basis. The two coincide for
all three bits when the corresponding self-consistent symmetry is conserved,
but a build may break one: then the states of both characters are merged into
a single storage block and the corresponding quantum number is not *good*.
The classic case is parity (section 2a): parity-broken builds store
everything in blocks 1, 2, 5, 6, and the parity bit of a basis sector is
vacuous.

Call the discrete quantum numbers that actually label the basis of a given
build its good quantum numbers, encoded as a 3-bit mask with the same bit
assignment as the sector words:

```
bit 0 : q(Rz) good    - set iff this build conserves Rz
bit 1 : parity good  - set iff this build conserves P
bit 2 : Tz good      - always set for nucleonic bases
```

Rule, applied once, inside `init` (section 4):

```
effective word = AND(physical word, goodqn)  ! bits of quantum numbers that
                                             ! are not good act diagonally
s' = IEOR(s, effective word)
```

- odd-l operator in a parity-broken build: the mask clears the parity bit,
  the blocks are diagonal - exactly what the by-hand `$PBROKEN` guards
  achieve in the current `Qlm_spme`.
- odd-m operator: the signature bit is kept, and populated sectors pair with
  their signature partners - nonzero when both signatures are stored
  (TR-broken builds: the wobbling-relevant operators), and exactly zero when
  only q(Rz) = +i is stored (TR-conserving builds). No special case.
- parity-conserving build whose single-particle space happens to contain no
  odd-parity states: the odd-parity sectors are empty, so an odd-l operator
  yields the (correct) zero matrix through the empty-sector skip.

The mask is a fact of the build scenario and cannot be recovered from the
runtime block counts: `[a,b,0,0,p,q,0,0]` arises both in a parity-broken
build - where an odd-l operator is nonzero within blocks 1, 2, 5, 6 - and in
a parity-conserving, TR-broken build whose basis simply lacks odd-parity
states, where the same operator gives exactly zero. In MOCCa the mask is a
`geninfo` constant injected by Hephaestos at generation time (`$GOODQN`,
computed from the build's symmetry generators by the same identity check the
`$PCONSERVED` substitution uses - bit 0 iff the Rz generator is present, bit
1 iff the P generator is, which avoids e.g. the 'P' inside a 'PT' generator
when parsing SYMSTRING tokens), and handed to the constructor by every
caller. These three bits are the only build-dependent information on this
path: the routine of step 3 (section 13) contains no `$`-marker logic at
all.

## 3. Type definitions

Per-block storage is a rank-2 allocatable, not an offset into one flat buffer.
With at most 8 blocks this costs nothing, keeps each block contiguous with a
natural leading dimension, and makes population callbacks trivial (direct
rank-2 view of the destination).

```fortran
type :: sym_block
  integer :: srow, scol      ! sector indices 0..7 (bookkeeping for construction)
  integer :: is, js          ! 1-based global start of the row/column sectors
  real(8), allocatable :: A(:, :)   ! (m, n), contiguous, LDA = m
  logical :: transposed     ! .true. if this block is reached as +/- transpose
                            ! of a stored partner (Hermitian pairing)
end type

integer, parameter :: KIND_GENERAL      = 0
integer, parameter :: KIND_SYMMETRIC    = 1   ! A(s',s) =   A(s,s')^T
integer, parameter :: KIND_ANTISYMMETRIC = -1 ! A(s',s) = - A(s,s')^T

type :: sym_matrix
  integer :: symmetry       ! operator symmetry word, 0..7 (physical, §2b)
  integer :: goodqn         ! good-quantum-number mask of the basis (§2b);
                            ! 7 = all quantum numbers good
  integer :: kind            ! KIND_GENERAL / KIND_SYMMETRIC / KIND_ANTISYMMETRIC
  integer :: n               ! total basis size, sum(HFblocks)
  integer :: nblocks         ! number of stored blocks, <= 8
  type(sym_block) :: blocks(8)
  integer :: st(0:8)         ! first 1-based global index of each sector,
                             ! st(8) = n+1 (section 3a, scheme 2)
contains
  procedure :: init, fill_blocks, fill_elements, matvec, apply2, matmat, ...
end type
```

(`add_element`, shown in earlier drafts of this snippet, was dropped: the
element-wise path is the `element_filler` callback through fill_elements, see
section 5.)

The `kind` sign generalizes Hermitian/anti-Hermitian for real storage: partner
blocks are related by `A(s',s) = sigma * A(s,s')^T` with `sigma = ±1`. A block
whose sector is its own partner (`s' = s`, which happens per bit where the
operator word is 0) is square and internally symmetric or antisymmetric; store
it dense anyway and use the symmetry only as a consistency check (at 1000 x
1000 scale dense BLAS beats triangle packing).

### 3a. Indexing conventions (Fortran is 1-based; sector words stay 0-based)

Fortran array bounds are a choice. Distinguish three roles:

- Sector words (0..7) are bit patterns, not array indices. Keep them
  0-based so `IEOR(s, op)` applies directly and `s > sp` is meaningful.
  Declare the arrays they index with lower bound 0: `HFblocks(0:7)`,
  `off(0:8)` (or `st(0:8)`), so no `±1` adjustments leak into logic.
  A 1-based sector encoding would require `IEOR(s-1, op) + 1` at every use.
- Global-basis offsets interact with 1-based Fortran vectors. Two coherent
  schemes; pick one and use it everywhere:
  1. 0-based offsets: `off(0) = 0`, slices written `x(j0+1 : j0+n)`.
     Uniform with the sector convention; one `+1` per slice.
  2. 1-based bounds array (recommended): `st(0:8)` with `st(s)` the first
     global index of sector `s` and `st(8) = n+1`, slices
     `x(st(s) : st(s+1)-1)`. The `+1` appears once, in the cumulative sum.

Scheme 2 is adopted (2026-11-16): the matrix component is named `st` and
each block carries the 1-based starts `is`, `js` of its row/column sectors.
- The MOCCa boundary is the only place 1-based input appears: existing
  arrays such as the current `HFblocks` are 1-based with their own block
  ordering. Convert once in the constructor, which also fixes the
  bijection between MOCCa's block ordering and the 3-bit sector word.

## 4. Construction from HFblocks

```fortran
subroutine init(this, HFBlocks, symmetry, kind, goodqn)
  ! HFBlocks(1:8): the GLOBAL MOCCa count array (HFBlocks_global),
  ! in MOCCa block ordering; sector s = block s+1 (section 2a).
  ! goodqn:        optional good-quantum-number mask of the basis
  !                (section 2b); default 7 = all quantum numbers good.
  ! 1. counts(s) = HFBlocks(s+1)   <- the only 1-based-to-0-based step.
  !    st(0) = 1; st(s+1) = st(s) + counts(s).
  ! 2. word = AND(symmetry, goodqn)  <- the only place the mask of section
  !    2b is applied; bits of quantum numbers that are not good act
  !    diagonally.
  !    For s = 0..7: sp = IEOR(s, word).
  !    Skip if counts(s) == 0 or counts(sp) == 0.
  !    If kind /= KIND_GENERAL and s > sp: skip (partner access covers it).
  !    Otherwise append a block with srow=sp, scol=s, is=st(sp), js=st(s),
  !    allocate A(counts(sp), counts(s)) = 0.
end subroutine
```

Properties after init:

- At most 8 blocks; exactly one stored block per surviving sector pair.
- For symmetric/antisymmetric matrices with nonzero symmetry word, blocks
  come in sector pairs and only one per pair is stored; the partner access
  flips `transposed` and applies the sign.
- For a totally symmetric operator (`symmetry = 0`), all blocks are diagonal
  in sector space (square, internally symmetric/antisymmetric).
- A purely real operator can simply leave the imaginary part absent (see
  section 8) without special cases.
- `symmetry` stores the physical word as handed in; the pairing word derived
  in init is `AND(symmetry, goodqn)` (section 2b). With the default
  `goodqn = 7` the pairing reduces to the pure XOR of section 2.

Index-to-sector mapping for element-wise construction is an 8-entry scan:

```fortran
pure integer function sector_of(this, i)
  ! 1 <= i <= n; returns -1 for indices outside the basis
  do s = 7, 0, -1
    if (i >= this%st(s)) return s
end function
```

## 5. Populating the blocks

The clean association mechanism is a callback through a dummy procedure with
an abstract interface; do not store procedure pointers as components.

```fortran
abstract interface
  subroutine block_filler(srow, scol, A, params)
    integer, intent(in) :: srow, scol   ! sectors 0..7 of row/column
    real(8), intent(out) :: A(:, :)     ! block to overwrite entirely
    class(*), intent(in) :: params      ! optional state (fields, couplings, ...)
  end subroutine
end interface

subroutine fill_blocks(this, fill, params)
  ! For each stored block b:
  !   call fill(blocks(b)%srow, blocks(b)%scol, blocks(b)%A, params)
  ! OpenMP parallel over blocks is safe: distinct destinations.
end subroutine
```

Guarantees and rules:

- The filler is called only for stored blocks. Hermitian/anti-Hermitian
  partner consistency is the type's responsibility, not the filler's; the
  type either maintains partners itself or handles them at access time.
- The filler receives `srow`/`scol` so it can apply the correct symmetry
  phase factors (this is where MOCCa's signature/isospin phase conventions
  for matrix elements enter).
- `params` carries state (fields, densities, coupling constants). Keep the
  signature fixed; extend `params`, not the interface.
- Element-wise fallback, confirmed form (2026-11-17): the per-element
  callback `element_filler(srow, scol, i, j, params)` returns `A_ij` at the
  LOCAL row/column indices (i, j) of the block; fill_elements queries it
  for every element of every stored block and OVERWRITES the storage. It
  costs one call per element, so it stays useful for tests and prototypes
  only. An `add_element` addressing global (i, j) was considered and
  dropped: routing non-stored partner elements would either double-count
  (caller adds both sides) or offload partner logic to callers, and partner
  consistency is anyway resolved at access time, never by writing.

Constructor sugar: `new_sym_matrix(HFblocks, symmetry, kind, fill, params)`
builds and fills in one statement (implemented 2026-11-17). A named
function, not a generic `sym_matrix(...)`: in Fortran a name cannot serve
as both a derived-type name and a function name. Equivalent to init
followed by fill_blocks with the same arguments. Step 3 (section 13)
extends the sugar's signature with the optional `goodqn` argument,
exactly like `init`.

## 6. Runtime operations

All runtime code loops over `blocks(1:nblocks)` only. No quantum numbers, no
XOR, no selection-rule checks.

Mat-vec (`y = A x`, full vectors of length n):

- For each block with `.not. transposed`: `dgemv('N', ...)` accumulating into
  `y(i0+1 : i0+m)`.
- For each block with `transposed`: `dgemv('T', ...)` with the sign from
  `kind`, scattering into the partner row range.
- OpenMP: parallelize over non-transposed blocks (independent output
  slices); handle transposed scatter separately (serial pass or atomics).

Mat-mat (`C = A B`): if `A` has symmetry word `a` and `B` has `b`, then
`C` has word `IEOR(a, b)` and one block per sector:

```
C(s, s⊕c) = A(s, s⊕a) · B(s⊕a, s⊕a⊕b),   c = a⊕b
```

implemented as at most 8 `dgemm` calls, with block lookup by scanning the
two 8-entry block lists once at init of the product.

Also needed eventually: transpose/conjugate (same symmetry word, blocks
transposed, sign rules from `kind`), addition (assert equal symmetry and
kind), scalar and complex scaling, trace.

## 7. Performance notes

- Link against an optimized real BLAS (OpenBLAS/MKL); all kernel calls are
  `dgemv`/`dgemm`.
- Each block is a separate contiguous rank-2 array: good locality, natural
  LDA, no stride games.
- Parallelize fill and mat-vec over the (at most 8) blocks.
- The element-wise path (`fill_elements`) is for correctness and tooling;
  hot paths should fill whole blocks via the callback.
- At the intended scale (total n around 1000, blocks around 100-500 per
  side) per-block dense storage is the right tradeoff. Triangle packing
  inside blocks is deliberately not done.

## 8. Complex operators: real/imaginary split

Represent a complex operator as two real matrices, never as complex blocks:

```fortran
type :: complex_op
  type(sym_matrix) :: re, im   ! each may be independently absent
end type
```

The algebra, with `O = R + iI`, `R` and `I` real:

- O Hermitian       : R symmetric,      I antisymmetric
- O anti-Hermitian  : R antisymmetric,  I symmetric

Each half is a `sym_matrix` with its own `kind`. The mapping from the complex
operator's Hermiticity class to the pair of real kinds is fixed:

- Hermitian:      (re: KIND_SYMMETRIC,     im: KIND_ANTISYMMETRIC)
- Anti-Hermitian: (re: KIND_ANTISYMMETRIC, im: KIND_SYMMETRIC)

Symmetry words: for unitary symmetries both halves share the operator word.
If a symmetry is antiunitary or involves complex conjugation in the adopted
convention, the imaginary part carries that bit complemented:
`im%symmetry = IEOR(op, antiunit_mask)`. Decide the mask once from MOCCa's
q-value conventions; this is the mechanism behind FAM facts such as fields
related by complex conjugation connecting complementary signature sectors.

Mat-vec with packed two-column vectors. With `x = x_R + i x_I`, `y = O x`:

```
y_R = R x_R - I x_I
y_I = I x_R + R x_I
```

Pack `X = [x_R | x_I]` as `n x 2`; each half of the operator is one `dgemm`
pass per block (`apply2`), with the sign/swap folding the multiplication by
`i`. Two real passes over an `n x 2` array stay in real BLAS and beat one
complex pass on memory traffic. The same packed layout serves iterative
solvers over the complex field (mat-vec closure over length-2n vectors).

Mat-mat: `C_R = A_R B_R - A_I B_I`, `C_I = A_R B_I + A_I B_R`, each term a
real sym_matrix mat-mat with XOR'd symmetry words as in section 6.

A purely real operator leaves `im` unallocated (`nblocks = 0`); no special
cases anywhere else.

## 9. Testing and validation

- Selection rules: for every operator word 0..7, verify the block list
  matches the XOR pairing of sectors surviving in `HFblocks`.
- Symmetry checks: for KIND_SYMMETRIC/ANTISYMMETRIC, verify partner-block
  consistency including the sign and the transposed access path.
- Mat-vec correctness: compare against a dense reference mat-vec built by
  expanding all blocks into a full `n x n` matrix, for random `HFblocks`
  (including empty sectors), all symmetry words, all kinds, real and
  complex.
- Mat-mat: compare against dense multiplication; verify symmetry word of the
  product is `IEOR(a, b)`.
- Complex split: verify that reassembling `R + iI` and applying as complex
  `zgemv` reproduces the packed two-column result, and that Hermitian /
  anti-Hermitian complex operators produce the promised real kind pairs.
- Benchmarks: time mat-vec against a dense `dgemv` on the same data;
  profile with likwid or gprof; check OpenMP scaling over blocks.

## 10. MOCCa integration notes

- Never edit generated solver code. The type belongs in `src_orig/` (solver
  source), and any build-system or preprocessing adjustments in
  `Hephaestos.py` / `src_heph/` / `configs/`. The module should be plain
  Fortran 2003+ with no preprocessor dependencies so it can be included
  wherever needed.
- Quantum-number conventions (bit assignments for parity, signature,
  isospin; phase factors under signature; the antiunitary mask of section 8)
  must be taken from MOCCa's existing symmetry machinery for q-values, not
  re-invented. The sector bits and the XOR pairing are internal to the
  module; only the constructor converts MOCCa conventions to sector words.
- The natural consumers are the linear-response (FAM) machinery, where
  response matrices and external fields with various symmetry words appear,
  and any place where dense operator matrices on the single-particle basis
  are built explicitly.
- Module name: `block_matrix_module` in `src_orig/block_matrix_module.f90`,
  following MOCCa naming style; `HFblocks` is the existing concept aligned
  with in the constructor. The file must be listed in `no_process_files`
  in `src_heph/preprocess.py` (and in `FORTRANFILES` in `Hephaestos.py` for
  the dry-run report) so Hephaestos copies it verbatim into the build tree;
  an unlisted file is silently skipped by the Makefile-per-file invocation.

## 11. Suggested implementation order

1. Core module: sector encoding, `sym_block`, `sym_matrix`, `init` from
   `HFblocks`, block-list construction with XOR pairing and kind-based
   partner skipping.
2. `fill_blocks` with the abstract interface, plus the element-wise
   fallback.
3. Generalised multipole matrix elements (detailed plan in section 13):
   the first production consumer. Move the `Qlm_spme`-based machinery from
   `fission_MOI.f90` to `moments.f90` and rebuild it on `sym_matrix` for
   all (l, m, Re/Im) electric multipoles in natural units `fm^l`, using the
   good-quantum-number mask of section 2b. The legacy routine stays in
   place as the testing oracle. Also lands the mask in `init`/`geninfo`
   and the dense-expansion helper in the module.
4. `matvec` (including transposed access with sign), tested against dense
   reference.
5. Symmetric/antisymmetric partner maintenance and consistency checks.
6. `matmat`, add, scale, trace.
7. `complex_op` wrapper with the packed two-column `matvec` and the
   antiunitary-mask constructor logic.
8. OpenMP over blocks for fill and mat-vec; benchmarks.
9. Integration into a FAM-relevant code path, behind a Hephaestos config
   option if appropriate.

## 12. Design decisions to keep

- At most 8 blocks; runtime code never touches quantum numbers.
- Selection rules are one XOR, resolved entirely in `init`.
- Per-block rank-2 allocatables, not one flat buffer.
- Population via dummy-procedure callback with abstract interface; no
   procedure-pointer components.
- Complex operators are pairs of real matrices; re/im carry independent
   symmetry words and kinds; packed two-column vectors for mat-vec.
- Dense blocks, no triangle packing; single node only for now.
- MOCCa boundary: `init` takes the global `HFBlocks_global(1:8)` directly;
  sector word = MOCCa block index - 1 (section 2a). One constructor, one
  place where the bijection appears.
- Offset scheme 2 throughout: 1-based bounds `st(0:8)` with `st(8) = n+1`,
  block starts `is`, `js`; no 0-based global offsets anywhere.
- Operator words handed to `init` are physical (section 2); the
  good-quantum-number mask of section 2b is applied once, inside the
  constructor. MOCCa's mask lives in `geninfo`, derived from the build's
  `SYMSTRING`; producers and consumers of these matrices carry no
  `$`-marker logic.

## 13. Step 3 detailed: generalised multipole matrix elements

Target (agreed 2026-11-17): move the single-particle matrix elements of the
electric multipole operators from `fission_MOI.f90` into `moments.f90`,
rebuild them on `sym_matrix`, generalise to all multipoles in natural
units, and keep the old routine in place as a testing oracle. Consumers
found in the tree: `fission_MOI` (collective inertia tensor), `functional`
(`CompCOMCorrection`, Q_10), `fam` (`get_qpme_op`: multipole fields, z_com),
`fam_testing` ([R, P] unit tests).

### 13.1 Semantics of the new generator

Operator: multiplication by `r^l Re Y_lm` or `r^l Im Y_lm`, evaluated on
the COM-shifted mesh (`meshx_shifted/...`, exactly as the current
`Qlm_spme`), in natural units `fm^l` (the `b^{l/2}` rescaling is dropped).
Loud input checks (`l <= MaxMoment`, `0 <= m <= l`);
symmetry-forbidden combinations return zero matrices silently - guarding
user input against forbidden moments is the business of Hephaestos' moment
list at the input layer (e.g. `read_inertia`), not of this routine.

Output: the full complex ME matrix in the HF basis, split per section 8
into two `sym_matrix` halves (routine names provisional, final naming at
implementation):

```fortran
! <i|Q_lm|j> = me_re + i*me_im : me_re KIND_SYMMETRIC, me_im KIND_ANTISYMMETRIC
! me_im is optional (skip where the build kills it; one extra block fill)
subroutine Qlm_spme_matrix(l, m, Im, me_re, me_im)
! Legacy dense shape of the real symmetric half: drop-in replacement,
! real(nwt,nwt) as before but in fm^l
function Qlm_spme_fm(l, m, Imaginary) result(me)
```

- physical symmetry word (section 2): `op = mod(m,2) + 2*mod(l,2)` - bit 0
  flips q(Rz) for odd m, bit 1 is the parity of the operator, bit 2 is
  always clear (electric multipoles are isoscalar). Both operator halves
  (Re/Im Y_lm) carry the same word: MOCCa's Rz and P are linear, no
  antiunitary complement (section 8) arises for electric multipoles.
- constructor call: `init(HFBlocks_global, op, kind, goodqn = ...)`, with
  the mask of section 2b taken from `geninfo`. Every build then takes the
  identical, marker-free code path; empty partner sectors give zero blocks
  through init's occupancy skip.
- all matrix elements are direct mesh integrals of the stored basis
  functions - no time-reversal/signature-partner shortcuts, hence no
  partner-phase conventions to get wrong. Block fill: one dense product
  per spin component per block (`A = sum_k X_k^T (Q dv) X_k`, a `dgemm`
  equivalent), with `dv` folded into the weight mesh; the antisymmetric
  imaginary half crosses the Re/Im wavefunction components the same way.
- helper added to `block_matrix_module`: a partner-reconstructing dense
  expansion of a `sym_matrix` (refactored out of the `print_dense` core).
  It is the return shape behind `Qlm_spme_fm` and later the dense
  reference builder for the mat-vec/mat-mat tests of steps 4 and 6.

The antisymmetric half is computed in all builds when asked:
TR-conserving builds yield roundoff-level zeros there, TR-broken builds
genuinely need it. `$TR`/`$NTR` markers are therefore absent from the
numerics; a short-circuit for TR-conserving builds remains an optional
optimisation for later. The old routine's restrictions (even m only, real
matrix part only, `b^{l/2}` units, marker-guarded selection) all
disappear.

### 13.2 The old routine stays, as oracle

The current `Qlm_spme` moves `verbatim` from `fission_MOI.f90` to
`moments.f90` (a pure move; behaviour preserved bit-for-bit) and remains
until step 4 (mat-vec) has landed and survived a shakeout period. It is
the oracle for the agreements below; the conversion identity is
`old [b^{l/2}] = new [fm^l] * 100^(-l/2)`.

### 13.3 Sequencing (one commit per row)

| # | Commit | Contents | Regression net |
|---|---|---|---|
| 0 | park FAM anchor | build `exec/fam.NLO.reference.exe` from the un-refactored tree and park it | none needed |
| 1 | goodqn mask | optional `goodqn` in `init` (+ module header), `$GOODQN` injected by Hephaestos into `geninfo`, masked-init unit tests | `block_matrix_unit_tests`; compile NLO, NLO-T, NLO-P, NLO-TP |
| 2 | pure move | `Qlm_spme` verbatim to `moments.f90` (which gains explicit `use wavefunctions`), deprecation banner | `fission.sh`, `minimal.sh` golden: outputs bit-identical |
| 3 | generator | dense-expansion helper (with its unit tests), `Qlm_spme_matrix`, `Qlm_spme_fm`, and an off-by-default oracle hook in `fission_MOI`: compute both paths for every requested moment, print the max deviation, run symmetry self-checks, compare `me_im` against a brute-force conjugated double-loop reference | unit tests; oracle runs NLO / NLO-T / NLO-P; all four configs compile |
| 4 | inertia switch | `calc_collective_inertia` takes `Qsp` from `Qlm_spme_fm * 100^(-l/2)`: the boundary rescale holds every presented number in `b` units | `fission.sh`, `minimal.sh` golden within tolerances |
| 5 | units flip | drop the boundary rescale: inertia chain in `fm` units end to end; `print_collective_inertia` convention strings updated; `fission.sh` reference inertias scaled by `100^(-(l_i+l_j)/2)` (I_2020 x 1e-4, I_2030 x 1e-5, I_3030 x 1e-6) and its header text | `fission.sh` |
| 6 | consumer migration | `functional` COM correction (drop `*10`), `fam` multipole/z_com operators (drop `*(100**(l/2))`, `*sqrt(100)`), `fam_testing` [R, P] (drop `*sqrt(100)`) | `fission.sh`, `minimal.sh` golden, `cranking.sh` (NLO-T), `fam_unit_tests`, local runs vs the parked FAM reference |

`Ksum_Mij_HF/BCS/HFB`, `calc_Q20` and the LAPACK inversion are
deliberately untouched in this step beyond the source of `Qsp`.

### 13.4 Validation and risks

- Oracle agreement to machine precision modulo summation order (BLAS vs
  the old scalar loops - the golden tolerances of `golden_diff.sh` exist
  for exactly this reassociation drift).
- The antisymmetric half has no legacy oracle (the old code computes the
  real part only): covered by the brute-force conjugated reference inside
  the hook, antisymmetry and Hermiticity identities, the [R, P] checks in
  `fam_testing`, and vanish-to-roundoff in TR-conserving runs.
- Parity-mask correctness: masked-init unit tests plus the NLO-P oracle
  run - the old routine's `$PBROKEN` diagonal-blocks path is the reference
  there.
- No available config breaks Rz today (all four SYMSTRINGs contain 'Rz'),
  yet the mask derivation is generic: an Rz-broken build clears bit 0 and
  odd-m operators act within the merged signature blocks - same rule, no
  new code.
- The units-flip arithmetic of the `fission.sh` references is verified by
  computation before editing the file; `make.inc` stays frozen (the trust
  anchor is untouched).
- MPI out of scope (the sprint is serial-only): construction takes
  `HFBlocks_global` and the dense wrappers assume the whole basis on one
  rank, matching `fission_MOI`'s current inert-under-MPI behaviour
  (`N_inertia = 0`).
- Park the FAM reference BEFORE any `fam.f90` edit: commit 0 precedes
  commit 6.

### 13.5 Deferred (not this step)

- consuming the re/im halves inside FAM response assembly (needs
  `complex_op`, step 7),
- odd-m and Im-operator FAM fields end to end,
- inertia for Im-operator moments (`Ksum`/`calc_Q20` generalisation),
- the optional `$TR` short-circuit for the antisymmetric half,
- deleting the old `Qlm_spme`.