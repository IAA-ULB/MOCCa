# Symmetry-Block Matrix Type for MOCCa: Design Summary

This document summarizes a design for a dense, block-structured matrix type to be
implemented in the MOCCa repository (github.com/IAA-ULB/MOCCa). It is written so
that an agent (VIBE CLI) or a collaborator can start implementation without the
preceding discussion.

## 0. Implementation status

This document records design intent, agreed conventions, amendments and
measured records - each amendment stays inline in the section it amends,
dated and commit-referenced. It deliberately carries no running progress
narrative: live step status (one row per landed step, commit hashes and
executed gates) is in `SPRINT.md`, open user decisions in `DECISIONS.md`,
and out-of-scope findings in `FOUND_BUGS.md` - all at the repo root.

## 1. Purpose and scope

MOCCa calculations are done in a basis of single-particle wavefunctions (spwfs)
that are symmetry-adapted and stored in a fixed order. Operators acting on this
basis respect three discrete symmetries, so an operator matrix decomposes into a
small number of dense blocks. The goal is a derived type that:

- stores only the nonzero blocks dictated by selection rules,
- keeps all block data contiguous and BLAS-friendly,
- hands the stored blocks to consumers that contract them with
  densities, amplitudes and other operator matrices at dense-BLAS speed,
- supports Hermitian and anti-Hermitian operators,
- has zero quantum-number logic at runtime (all resolved at construction).

Mat-vec and mat-mat operations were part of the original goal and are
dropped by the consumer audit of 2026-11-19 (section 6).

Constraints agreed with Wouter:

- Basis states are ordered; a simple integer array `HFblocks(0:7)` gives the
  number of states in each symmetry sector, in storage order.
- Isospin is limited to `T_z = ±1/2` (proton / neutron).
- Blocks are dense; matrices are constructed explicitly, in full.
- Total matrix sizes are typically around 1000 x 1000; single-node performance
  matters most. MPI distribution is out of scope for now.
- STy conservation can be assumed in all cases.

## 2. Symmetry model

The three symmetries that induced block structure in all matrices 
are parity (P), z-signature (Rz), and isospin (T_z). Each acts on 
the basis as a Z2 grading:

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

### 2a. The y-time-simplex characteristic (agreed 2026-11-17)

The y-time-simplex is STy = P x T x Ry, parity times time-reversal times
y-signature: a hermitian, antilinear involution, STy**2 = +1. The spwfs
MOCCa stores are its invariants,

  STy \psi = \psi^*(x,-y,z,sigma) = \psi .

For an operator O of definite characteristic STy O STy^-1 = sigma x O,
sigma = +/- 1, the matrix elements over the stored basis therefore obey

  < i | O | j > = < STy \psi_i | O | STy \psi_j > = sigma x < i | O | j >* ,

so STY_EVEN operators (sigma = +1) have purely real and STY_ODD operators
(sigma = -1) purely imaginary matrix elements. For even (odd) operators the
one-body matrix elements are the REAL (IMAGINARY) part of the spinor inner
product of the spwfs.

STy introduces NO block structure: every stored state is an invariant, so
there is no grading to pair, no XOR bit, and the pairing machinery of section
2 is untouched. The characteristic is a phase constraint on the stored
numbers and joins Hermiticity on the structure side of the type:

- Classifiability (agreed 2026-11-18): every operator handled by the type is
  classifiable - Hermitian or anti-Hermitian, and STy-even or STy-odd. Both
  classes are MANDATORY `init` inputs, `herm` (HERMITIAN / ANTIHERMITIAN)
  and `sty` (STY_EVEN / STY_ODD); there is no STY_NONE and no general,
  non-(anti-)Hermitian kind input.
- The storage kind is DERIVED from the classes as

      kind = herm x sty

  i.e. the stored numbers - the real part (STY_EVEN) or the imaginary part
  (STY_ODD) of the matrix elements - form a symmetric matrix for a
  Hermitian STY_EVEN or an anti-Hermitian STY_ODD operator, and an
  antisymmetric one for a Hermitian STY_ODD or an anti-Hermitian STY_EVEN
  operator. Hermiticity makes the real part symmetric and the imaginary
  part antisymmetric; an anti-Hermitian operator flips both.
- Completeness: in a build that conserves the time simplex (goodqn bit 3,
  the case for every SYMSTRING in the tree) the omitted half of the
  complex matrix elements is identically zero and the single stored matrix
  IS the complete operator. With bit 3 clear the same inputs still define
  a consistent single matrix - the selected half of the (re, im) pair
  model of section 12, the other half a separate matrix - but no
  completeness claim. The mask is stored verbatim by init and read by
  consumers as the completeness criterion; init does not gate on it.
- Producers: fill callbacks compute the real or the imaginary part of the
  spinor product according to the characteristic; the callback interface
  is frozen (section 5), so the producer hands the same value to its
  callback through `params`.
- Consumers (audited 2026-11-19): every consumer in the tree contracts
  the matrices with other matrices or reads them elementwise - traces
  against densities (moments, CompCOMCorrection), elementwise products
  of two operator matrices (the `Ksum_Mij` ATDHFB sums), elementwise
  contractions with the FAM amplitude matrices (`X = -(F + dH)`, the
  strength sums), and within-block matmuls against consumer-owned
  Bogoliubov/canonical transformations. No consumer applies an operator
  to a basis vector, so the planned mat-vec combination of packed
  two-column complex vectors is dropped (section 6).
- Caveat for composed objects: the classifiability assumption covers
  physical one-body operators. It does NOT cover products of two
  non-commuting classifiable operators (AB of Hermitian A and B is
  neither Hermitian nor anti-Hermitian: AB = {A,B}/2 + [A,B]/2) - one
  more reason the mat-mat product is out of scope entirely (dropped
  2026-11-19, section 6): no consumer forms operator products, and
  one-body Q^2 through sp matrices is impossible in any MOCCa basis
  (the truncated basis does not resolve Q|i>).
- Character table, electric multipoles (agreed 2026-11-17): r^l Re Y_lm is
  STY_EVEN and r^l Im Y_lm is STY_ODD for every m - the characteristic is
  m-independent, unlike the symmetry word, which flips the signature bit
  for odd m. J_y is STY_ODD. Character tables are producer-side knowledge;
  for operators acting on spin they depend on MOCCa's spin convention for
  STy and must be taken from the existing symmetry machinery, not
  re-derived.

Two convention safeguards. The design presumes the +1 invariance above: a
basis holding STy \psi = -psi states would need per-state signs theta_i
theta_j in the reality constraint, which MOCCa's conventions do not
produce. And a single-particle state cannot simultaneously be an Rz
eigenstate and a sigma_z = Rz.T invariant - MOCCa stores Rz eigenstates,
so no second, z-simplex-based reality constraint exists and STy is the one
antiunitary characteristic the type carries.


### 2b. Bijection with MOCCa's HFBlocks bookkeeping (confirmed)

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
- Parity broken ('Rz,T,STy' / 'Rz,STy'): the eight-slot grid PERMITS
  blocks 1, 2, 5, 6 (sectors 0, 1, 4, 5), but at runtime only blocks
  1, 5 are populated: the parity blocks MERGE into them, every stored
  state is q(Rz) = +i (user-confirmed 2026-11-19), and blocks 2, 6
  carry zero count for as long as T is conserved; blocks 3, 4, 7, 8
  remain empty throughout the run and the parity bit of a basis sector
  is vacuous. The s = block - 1 encoding is thereby untouched in every
  build.

Convention note for antiunitary operations: if any of the three symmetry
operations is antiunitary (involves complex conjugation) in the convention
used, the imaginary part of a complex operator carries that bit complemented
(see section 8). The type supports this by letting the real and imaginary
parts carry independent symmetry words; it must be decided once, in the
constructor, based on MOCCa's existing symmetry conventions for q-values.

### 2c. Good quantum numbers of the stored basis (agreed 2026-11-17)

The operator word of section 2 describes the *physical* character of an
operator. The XOR pairing additionally presumes that the sector words
describe definite quantum numbers of the stored basis. The two coincide for
all three bits when the corresponding self-consistent symmetry is conserved,
but a build may break one: then the states of both characters are merged into
a single storage block and the corresponding quantum number is not *good*.
The classic case is parity (section 2): parity-broken builds store
everything in blocks 1, 2, 5, 6, and the parity bit of a basis sector is
vacuous.

Call the discrete good quantum numbers that actually label the basis of a
given build its good quantum numbers, encoded as a 4-bit mask: bits 0-2
with the same bit assignment as the sector words, bit 3 for the
y-time-simplex characteristic of section 2a:

```
bit 0 : q(Rz) good    - set iff this build conserves Rz
bit 1 : parity good  - set iff this build conserves P
bit 2 : Tz good      - always set for nucleonic bases
bit 3 : STy good     - set iff this build conserves the y-time-simplex,
                       whose invariants the stored spwfs then are
```

Bits 0-2 and bit 3 play different roles, applied exactly once,
inside `init` (section 4):

```
effective word = AND(physical word, goodqn)  ! bits of quantum numbers that
                                             ! are not good act diagonally
s' = IEOR(s, effective word)

goodqn itself is stored verbatim; bit 3 is NOT a pairing bit - it is the
completeness criterion for the STy characteristic (section 2a), read by
consumers to know whether a tagged matrix is the entire operator
```

The physical classes `herm` and `sty` are mandatory init inputs and are
never masked: a classifiable operator keeps its classes in every build.
Bit 3 only determines whether the stored half is the complete (anti-)
Hermitian operator (bit set) or one half of the complex pair model
(section 12).

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
1 iff the P generator is, bit 3 iff the STy generator is, which avoids e.g.
the 'P' inside a 'PT' generator when parsing SYMSTRING tokens), and handed
to the constructor by every caller. These four bits are the only
build-dependent information on this path: the routine of step 3 (section 13)
contains no `$`-marker logic at all. Bit 3 does not mask the pairing (a
conserved time simplex neither merges nor separates blocks) and is not a
gate; it is the consumer's completeness criterion for the characteristic
of section 2a, and every SYMSTRING in the tree today declares
STy, so it is set in all existing configs.

### 2d. TR-folded composites (agreed 2026-11-19, design stage)

In a build that stores only the q(Rz) = +i sectors - every
parity-conserving config, blocks 1, 3, 5, 7 - the stored basis is one
half of the doubled single-particle basis, and the omitted -i partners
are the time reverses of the stored states (MOCCa's T is fermionic,
T**2 = -1, section 14). A `$TR`-marked consumer that needs the
cross-signature content of a signature-flipping operator can then not
read it off the honest operator matrix - the partner sectors are empty
and the honest matrix is structurally zero (the odd-m multipole
precedent). What it needs is the composite

    < i , O T , j >  =  < psi_i | O | T psi_j >  ,

the honest matrix element of O between a stored row state and the
time-reversal image of a stored column state: a piece of the honest
operator matrix over the doubled basis, folded back onto the stored half
- the matrix elements "of nabla with the time reverses of the states", in
the words of the agreement. Rather than handling this case per family,
every generator family exposes the folded twin through ONE uniform
convention:

- Construction: the honest kernel with the time-reversed image of the
  ket, exactly as the wavefunctions.f90 machinery already implements it
  - the `t` kernels for the angmom family, the TR flag of
  `Transfer_psi`/`Transfer_derpsi` for the gradient family. Producers
  wrap the existing kernels verbatim, values per stored pair, as
  everywhere in this sprint; no re-derivation of the fold conventions.
- Word: T conserves parity and isospin and flips the Rz grading, so the
  twin's physical word is `IEOR(word, 1)` - observed in both existing
  families (J_x word 1 -> word 0; \nabla_x word 3 -> word 2).
- Herm: the composite is HERMITIAN as an operator for the angmom family
  (measured 2026-11-19, section 14.1): with MOCCa's FERMIONIC time
  reversal, T**2 = -1, the user-accepted chain
  ` (J.T)^dagger = -T**2 (J.T) = +(J.T)` holds. The blanket ANTIHERMITIAN
  rule asserted earlier is REFUTED by the oracle: both folded angmom
  twins measure SYMMETRIC storage (kind +1). CAVEAT: the section-2a
  kind = herm * sty machinery does not transfer to antilinear composites
  - the sty-odd twin (K_y) needs herm = ANTIHERMITIAN as the
  kind-completing LABEL (herm = kind * sty), not as a physical class.
  Classes of future families (nabla section 15, spin, ...) remain
  per-family assertions, each settled by its own kernel-transpose oracle:
  kind signs are test assertions, never assumptions.
- Sty: carried over from the honest operator, per component. Kind is
  settled by measurement, not derived (the caveat above): the angmom
  twins measure kind +1 BOTH (symmetric storage, section 14.1); the
  \nabla twins of section 15.2 are pending their own oracle.
- Existence and marker policy: the twins EXIST in every build and the
  generator modules stay marker-free; their callers gate - only `$TR`
  -marked call sites construct them, dead code in TR-broken builds by
  construction, not by preprocessing.

The convention absorbs its own special cases. The twin of an honest
operator whose word has bit 0 clear (J_z word 0, \nabla_z word 2) has
bit 0 set: in a `$TR` build its partner sectors are empty and it is
structurally zero through the ordinary empty-sector skip - "needs
no folded twin" is a corollary, not a rule. In a `$TR` build exactly one
of {honest, twin} carries content per family, decided by bit 0 of the
honest word; in a `$NTR` build the honest objects carry everything and
no twin is ever constructed.

Every TR-conserving config stores the +i half only, parity-broken
included (user-corrected 2026-11-19): when parity breaks, the parity
blocks MERGE - the basis never gains -i states while T is conserved -
so the collapse above covers NLO-P just like NLO. The +i-only reading
of the composite - "the honest matrix elements with the time reverses
of the stored states" - is what every `$TR` call site consumes, and the
doubled basis materialising T-images alongside the stored states is the
reference shape for the oracle proposed in 15.6. A bookkeeping caution
(user-confirmed 2026-11-19): 2b's occupation row lists the grid-PERMITTED
slots - under merged parity only blocks 1 and 5 are populated at
runtime and every stored state is q(Rz) = +i, so the constructor's
`s = block - 1` encoding stands in every build.

Consumers address objects by name (honest vs folded), never through
build-dependent array slots - dissolving the fill-slot semantic overload
that CompNablaMelements' array carries today (section 15.1).

Scope: a PRODUCER-SIDE convention. The type changes not at all - a
folded composite is an ordinary classifiable matrix over the stored
basis (one stored half, its classes, its masked word), which is why
block_matrix_module must not learn about time reversal; T-knowledge
stays with the kernels of the generator modules. Families covered so
far: the angular-momentum generators of section 14 (K_x, K_y) and the
gradient generators of section 15 (\nabla_x.T, \nabla_y.T). Multipoles
are deliberately NOT folded: section 13.5 remains the gate for their
odd-m / imaginary-operator content end to end.

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
end type

integer, parameter :: HERMITIAN           =  1   ! physical class under O†
integer, parameter :: ANTIHERMITIAN     = -1   ! physical class under O†
integer, parameter :: STY_EVEN           =  1   ! <i|O|j> purely real (§2a)
integer, parameter :: STY_ODD            = -1  ! <i|O|j> purely imaginary (§2a)
integer, parameter :: KIND_SYMMETRIC    =  1    ! A(s',s) =   A(s,s')^T
integer, parameter :: KIND_ANTISYMMETRIC = -1   ! A(s',s) = - A(s,s')^T

type :: sym_matrix
  integer :: symmetry       ! operator symmetry word, 0..7 (physical, §2b)
  integer :: goodqn         ! good-quantum-number mask of the basis (§2c);
                            ! 0..15, bits 0-2 = sector gradings, bit 3 = STy
                            ! good (§2a); 15 = every quantum number good
  integer :: kind            ! storage kind DERIVED in init as kind =
                            ! herm * sty (§2a): KIND_SYMMETRIC or
                            ! KIND_ANTISYMMETRIC, never zero
  integer :: sty            ! y-time-simplex characteristic (§2a), the input
                            ! class verbatim: STY_EVEN / STY_ODD. No block
                            ! structure: tells producers and consumers which
                            ! part of the complex matrix elements the stored
                            ! real numbers are (real part / imaginary part).
                            ! herm is recoverable as kind * sty
  integer :: n               ! total basis size, sum(HFblocks)
  integer :: nblocks         ! number of stored blocks, <= 8
  type(sym_block) :: blocks(8)
  integer :: st(0:8)         ! first 1-based global index of each sector,
                             ! st(8) = n+1 (section 3a, scheme 2)
contains
  procedure :: init, sector_of, fill_blocks, dense, print_dense
end type
```

(`add_element`, shown in earlier drafts of this snippet, was dropped, and
with it (2026-11-20, DECISIONS.md D-6) the whole element-wise path: the
call audit found no caller of `fill_elements`/`element_filler` or the
`new_sym_matrix` sugar outside the block-matrix unit net, and every real
producer populates whole blocks. partner access is `dense()`'s kind
algebra alone. KIND_GENERAL existed in earlier drafts and was dropped with the
classifiability agreement of section 2a: the kind is derived from the two
mandatory physical classes and can never be zero.)

The (anti-)Hermitian operator relates partner blocks of the stored real
numbers by `A(s',s) = kind * A(s,s')^T` with `kind = ±1`. A block
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
subroutine init(this, HFBlocks, symmetry, herm, sty, goodqn)
  ! HFBlocks(1:8): the GLOBAL MOCCa count array (HFBlocks_global),
  ! in MOCCa block ordering; sector s = block s+1 (section 2a).
  ! herm:          mandatory physical class: HERMITIAN or ANTIHERMITIAN.
  ! sty:           mandatory y-time-simplex characteristic (section 2a):
  !                STY_EVEN (stored numbers = real part of the matrix
  !                elements) or STY_ODD (= imaginary part).
  ! goodqn:        optional good-quantum-number mask of the basis
  !                (section 2c); 0..15, bits 0-2 = sector gradings,
  !                bit 3 = STy good (section 2a), default 15; stored
  !                verbatim - bits 0-2 enter the pairing word, bit 3 is
  !                the consumers' completeness criterion for sty.
  ! 0. kind = herm * sty  <- the derived storage kind; stored in
  !    this%kind together with the verbatim characteristic in this%sty.
  ! 1. counts(s) = HFBlocks(s+1)   <- the only 1-based-to-0-based step.
  !    st(0) = 1; st(s+1) = st(s) + counts(s).
  ! 2. word = AND(symmetry, goodqn)  <- the only place the grading bits of
  !    the mask are applied; bits of quantum numbers that are not good act
  !    diagonally.
  !    For s = 0..7: sp = IEOR(s, word).
  !    Skip if counts(s) == 0 or counts(sp) == 0.
  !    If s > sp: skip (partner access covers it; the derived kind is
  !    always +/- 1).
  !    Otherwise append a block with srow=sp, scol=s, is=st(sp), js=st(s),
  !    allocate A(counts(sp), counts(s)) = 0.
end subroutine
```

Properties after init:

- At most 8 blocks; exactly one stored block per surviving sector pair.
- Blocks come in sector pairs and only one per pair is stored; the partner
  is reconstructed at access time by the derived sign and transposition
  (dense()).
- For a totally symmetric operator (`symmetry = 0`), all blocks are diagonal
  in sector space (square, internally symmetric/antisymmetric).
- A purely real operator can simply leave the imaginary part absent (see
  section 8) without special cases.
- `symmetry` stores the physical word as handed in; the pairing word derived
  in init is `AND(symmetry, goodqn)` on the grading bits (section 2c). With
  the default `goodqn = 15` the pairing reduces to the pure XOR of section 2.
- `kind` stores the DERIVED storage kind `herm * sty` of section 2a and is
  always +/- 1; `sty` stores the characteristic verbatim. The block list is
  independent of both. Producers fill the real (STY_EVEN) or imaginary
  (STY_ODD) part of the spinor product accordingly; consumers read the
  components to interpret the stored numbers and to combine packed complex
  vectors.

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
- Amendment (2026-11-20, DECISIONS.md D-6): the per-element fallback
  (`element_filler` through `fill_elements`) and the `new_sym_matrix`
  constructor sugar are DROPPED - the call audit found no caller outside
  the block-matrix unit net; every producer fills whole blocks through
  `block_filler`, so population is the single path above. `add_element`,
  addressing global (i, j), had already been considered and rejected:
  routing non-stored partner elements would either double-count or
  offload partner logic to callers, and partner consistency is anyway
  resolved at access time, never by writing.

## 6. Runtime operations

All runtime code loops over `blocks(1:nblocks)` only. No quantum numbers,
no XOR, no selection-rule checks.

Amendment (agreed 2026-11-19): the mat-vec and mat-mat steps of the
original design are DROPPED. A consumer audit over the tree found no
place where MOCCa applies a one-body operator matrix to a basis vector:

- moments, CompCOMCorrection, the moments unit net: traces Tr(rho Q)
  against the density, elementwise;
- fission_MOI `Ksum_Mij_*`: elementwise products Qa(i,j) Qb(i,j) summed
  over particle-hole pairs - the ATDHFB sums contract two operator
  matrices, they are never matrix products;
- FAM: `X = -(F + dH)` elementwise on the n x n amplitude matrices,
  strength sums `S = sum conjg(F(i,j)) X(i,j)` elementwise, Bogoliubov
  qp<->sp transforms as within-block matmuls against consumer-owned U/V
  blocks, the free response diagonal in the sp/qp energies, elementwise;
- momentsofinertia: within-block canonical/HF transforms plus
  Belyaev-type trace sums.

Root cause: MOCCa acts with fields on wavefunctions on the mesh through
the density route, never in HF coefficient space. Single-particle ME
matrices exist only to be contracted with other matrices or read
elementwise; `y = A x` over global basis vectors has no caller. Products
of operators have no caller either, and a product of two non-commuting
classifiable operators is itself not classifiable (section 2a), so
mat-mat would not even fit the type.

What remains of the runtime surface:

- `dense()` / `print_dense`: the partner-reconstructing full n x n
  expansion - the consumer shape behind `Qlm_spme_fm`. It is the only
  place partner blocks are synthesised at access time, applying the
  derived kind sign.
- partner-consistency CHECKS: a unit-test operation, not a runtime path -
  verify stored blocks obey `A(s',s) = kind * A(s,s')^T`.
- traces and elementwise contractions stay at the call sites, as today.

The transposed-partner access machinery that mat-vec would have needed
is not built; nothing else requires it. The `sym_block%transposed`
component drafted for it was removed with the element-wise path
(2026-11-20, DECISIONS.md D-6) - it was never set and every reader
ignored it; `dense()` reconstructs partners through the kind relation.

### 6a. Basis rotations (agreed 2026-11-19)

The type gains exactly one runtime operation beside `dense()`: the change
of basis of a stored operator by a REAL ORTHOGONAL transformation that is
sector-diagonal - one square orthonormal block U_s per populated sector s.
This is the shape MOCCa's rotation producers actually build: `Canonical`
(HFB.f90) diagonalises rho per MOCCa block (`do B=1,8`) and fills
`cantransfo` exactly block-diagonal sector-wise; the merged
signature-pair loops of `basis_transform`'s `transform_mat` are
conservative about this but mix nothing. In the classifiable,
STy-invariant setting (section 2a) the general unitary group collapses
to this orthogonal sector-diagonal shape: no other unitary preserves
the type's invariants. Out of scope, deliberately (agreed 2026-11-19):

- complex unitaries: a complex combination of STy-invariant spwfs is not
  an STy invariant (STy is antilinear), so the sty classes - and with
  them the whole stored-half model - would lose their meaning; there is
  no caller today.
- rotations mixing sector partners, e.g. the q(Rz) = +i/-i pairs: the
  rotated basis has no definite signature, so no sym_matrix can represent
  an operator over it unless the signature bit stops being a good
  quantum number (the goodqn-mask machinery of section 2c would cover
  that case if ever needed); no config in the tree breaks Rz and no
  producer makes such rotations.

Why this does not reopen the dropped mat-mat (section 6): operator
products AB remain out - the product of two non-commuting classifiable
operators is not classifiable. A similarity transform is a different
animal: it is the SAME operator expressed in a rotated basis, and for a
real sector-diagonal orthogonal C the transformation preserves
everything the type stores:

- the rotated basis states are real linear combinations of STy
  invariants, hence STy invariants themselves: the sty classes of all
  operators keep their meaning, and the stored numbers remain the
  selected (real/imaginary) half of the complex matrix elements;
- every stored block transforms independently,
      A'(s', s) = U_{s'}^T A(s', s) U_s  ,
  and the partner relation survives:
      A'(s, s') = U_s^T A(s, s') U_{s'}
                = U_s^T (kind * A(s', s)^T) U_{s'}
                = kind * (U_{s'}^T A(s', s) U_s)^T
                = kind * A'(s', s)^T  .
  herm, sty, kind, the symmetry word, goodqn, the block list and the
  st bounds are ALL unchanged: rotate touches numbers, never metadata.
  Two dense products per stored block (intrinsic matmul: the module
  stays free of explicit BLAS calls, so the standalone unit suite links
  no library); the unstored partner blocks are never
  touched, exactly as everywhere else in section 6.

Type and interface:

```fortran
! One square orthonormal block per sector; empty sectors stay unallocated.
type :: sector_rotation
  real(8), allocatable :: U(:, :)
end type

type :: sym_rotation
  integer :: n                    ! total basis size, sum(HFBlocks)
  integer :: st(0:8)              ! 1-based sector bounds, scheme 2 (3a)
  type(sector_rotation) :: U(0:7) ! populated sectors only
contains
  procedure :: init              ! (HFBlocks): identity per populated sector
  procedure :: init_from_dense   ! (HFBlocks, C, tol): extract the per-sector
                                 !  blocks from a full n x n rotation array
                                 !  (basis_transform convention); refuse a C
                                 !  with nonvanishing off-sector submatrices
                                 !  or non-orthonormal blocks
  procedure :: set_sector        ! hand one sector's block directly (tests);
                                 !  orthonormality checked on the same tol
end type
```

and one addition to sym_matrix's `contains` list:

```fortran
procedure :: rotate  ! in place: A_b <- U_{srow}^T A_b U_{scol} for every
                     !  stored block; refuses a rotation whose n and st do
                     !  not match the operand's exactly (same HFBlocks)
```

Conventions:

- C convention of `basis_transform`: each COLUMN of C is a new basis
  vector expressed in the old basis (`psi' = C^T psi`), and the matrix
  of the same operator in the rotated basis is `C^T A C`. `rotate`
  implements exactly this, per block. Consumer-side direction flips
  exist in the tree (e.g. functional.f90 applies `C P C^T` to some of its
  matrices - the inverse direction): callers wanting that hand `rotate`
  the sector-wise transpose of their C. The inverse rotation is always
  the sector-wise transpose (U_s orthogonal: U_s^T = U_s^{-1}); no
  separate interface is built.
- Composition: rotating by R1 and then by R2 is a single rotation by
  the sector-wise composite U_s = U1_s . U2_s; no compose helper is
  built - callers matmul the per-sector blocks, or simply rotate twice.
- goodqn plays no role: rotate acts on storage sectors, never on quantum
  numbers, and inherits every piece of metadata from the operand.

Checks (loud, the module's print-and-stop style), all inside the module:

- `init_from_dense`: every off-sector submatrix of C must vanish within
  the tolerance (proposed default 1e-12 absolute, a named argument;
  MOCCa's producers write exact zeros there today) and every extracted
  block must be orthonormal within it. A C that mixes sectors is
  REFUSED: the rotated matrix would not be representable, and accepting
  it would silently corrupt physics.
- `set_sector`: orthonormality of the handed block, same tolerance.
- `rotate`: the rotation's n and st must match the operand's exactly.

Consumers, now and later (agreed 2026-11-19: this is an infrastructure
step - no consumer migrates with it; the block_matrix unit net is the
net). The production candidates, in the natural order: the `U^T O U`
factor of the FAM qp-space operator in `get_qpme_op` (the full
Bogoliubov image `U^T O U - V^T O^T V +/- ...` is NOT a rotation - the
V terms carry transposes - and stays consumer-side unless its output is
ever made a sym_matrix over a qp bookkeeping), the canonical-basis
operator chains of the momentsofinertia / Ksum families, and the
whole-matrix `transform_mat` call sites once their callers migrate.

## 7. Performance notes

- Link against an optimized real BLAS (OpenBLAS/MKL); the kernel calls
  are the dense products inside the fill callbacks and at consumer call
  sites (`dgemm`, and `dgemv` only where a consumer owns a genuine vector).
- Each block is a separate contiguous rank-2 array: good locality, natural
  LDA, no stride games.
- Parallelize fill over the (at most 8) blocks.
- At the intended scale (total n around 1000, blocks around 100-500 per
  side) per-block dense storage is the right tradeoff. Triangle packing
  inside blocks is deliberately not done.

## 8. Testing and validation

- Selection rules: for every operator word 0..7, verify the block list
  matches the XOR pairing of sectors surviving in `HFblocks`.
- Symmetry checks: for KIND_SYMMETRIC/ANTISYMMETRIC, verify stored blocks
  obey `A(s',s) = kind * A(s,s')^T` (unit-level check; the runtime
  partner access path is not built, section 6).
- Mat-vec/mat-mat checks are dropped with the operations (section 6).
- Basis rotations (section 6a): identity is an exact no-op; the round trip
  (rotate by R, then by the sector-wise transposed rotation) recovers the
  blocks at the reassociation level of two dense products per block; the dense reference
  `C^T dense(A) C` agrees with the rotated blocks at their global
  positions; the partner relation `A'(s', s) = kind * A'(s, s')^T` holds
  for the rotated matrix of BOTH kinds (a KIND_ANTISYMMETRIC operand
  stays antisymmetric); the diagonal blocks' per-sector traces are
  invariant under rotation; both sty classes are exercised; the refusal
  paths (non-orthonormal U, off-sector-nonzero C, HFBlocks-mismatched
  rotate) stop loudly.
- Imaginary half of `Qlm_spme_matrix`: the full-box reconstruction oracle
  (section 13.1, validation) - selected off-diagonal elements recomputed
  by direct mesh quadrature of the chamber-extended spinors, judged
  against physics: purely imaginary oracle values for the STY_ODD
  operator, `<j|O|i> = <i|O|j>*` with the stored antisymmetry.
- Benchmarks: profile the fill callbacks' dense products; check OpenMP
  scaling over blocks.

## 9. MOCCa integration notes

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

## 10. Suggested implementation order

Live status of the steps below, with commit hashes and executed gates:
the step table in `SPRINT.md` (repo root) - not this section.

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
4. (amended 2026-11-19, replaces mat-vec) full-box reconstruction oracle
   for the imaginary half of `Qlm_spme_matrix`: selected off-diagonal
   elements recomputed by direct mesh quadrature of the chamber-extended
   spinors, as a development hook in `moments_testing` beside the moments
   unit net. Together with the partner-consistency checks this closes the
   validation gap the mat-vec step had been assigned.
5. FAM-side migration: `get_qpme_op` off the legacy b-units `Qlm_spme`
   (design step 9); regressions against the parked `fam.NLO.reference.exe`.
   Landed 2026-11-19 (commit 20fca2e5): both call sites ('multipole' l>0
   and 'zcom') feed from `Qlm_spme_fm` in fm^l. A/B against the parked
   anchor on a shared O18/SLy4/HFB state: Q20 and Q22 strength
   famfiles bit-identical, GMRES residuals identical; fam unit tests
   (SLy4) green; fam builds of NLO/NLO-T/NLO-P compile. The zcom call
   site has no working runtime gate - every parity-odd FAM response
   NaNs identically in the parked pre-refactor anchor (l=1 multipole,
   including the committed fam_spurious_mode_subtraction scenario;
   operator_type 'zcom') - so it is source-validated (moments full-box
   oracle covers `Qlm_spme_fm(1,0)` on all builds) with prefactors
   preserved verbatim. The imaginary half stays deferred (section 13.5).
6. (dropped 2026-11-19) mat-mat: no consumer, and products of two
   non-commuting classifiable operators are not classifiable (section 2a).
   The remaining small operations (add, scale, trace) land if and when a
   consumer needs them.
7. (agreed 2026-11-19, intermediate infrastructure step - precedes the
   deferred section 13.5 items) basis rotations: the `sym_rotation` type
   and the `rotate` method of sym_matrix (section 6a), plus the
   extension of the block_matrix unit net (section 8). No consumer
   migrates with this step; the first production caller is decided
   later (candidates listed in section 6a).
   Landed 2026-11-19: `sector_rotation`/`sym_rotation` (identity `init`,
   `init_from_dense` with the structure checks, `set_sector`), the
   in-place `rotate`, the orthonormality helper; the unit net extended
   to 53562 checks, plus the 8-scenario refusal program
   `block_matrix_refusal_tests.f90` - the structure guards stop the
   process, so the script drives them per exit status. Gates: unit net
   green (0 failures); golden minimal NLO 4849/4849 numbers
   bit-identical; the moments suite green on NLO / NLO-T / NLO-P
   including the cross-build comparison; mf builds of all four configs
   compile with the frozen make.inc.

## 12. Design decisions to keep

- At most 8 blocks; runtime code never touches quantum numbers.
- Selection rules are one XOR, resolved entirely in `init`.
- Per-block rank-2 allocatables, not one flat buffer.
- Population via dummy-procedure callback with abstract interface; no
   procedure-pointer components.
- Dense blocks, no triangle packing; single node only for now.
- MOCCa boundary: `init` takes the global `HFBlocks_global(1:8)` directly;
  sector word = MOCCa block index - 1 (section 2a). One constructor, one
  place where the bijection appears.
- Offset scheme 2 throughout: 1-based bounds `st(0:8)` with `st(8) = n+1`,
  block starts `is`, `js`; no 0-based global offsets anywhere.
- Operator words handed to `init` are physical (section 2); the
  good-quantum-number mask of section 2c is applied once, inside the
  constructor. MOCCa's mask lives in `geninfo`, derived from the build's
  `SYMSTRING`; producers and consumers of these matrices carry no
  `$`-marker logic.
- Every operator is classifiable (section 2a): `init` takes the mandatory
  physical classes `herm` (HERMITIAN / ANTIHERMITIAN) and `sty`
  (STY_EVEN / STY_ODD) and DERIVES the storage kind as kind = herm * sty.
  No general kind input exists. STy is the one antiunitary constraint in
  play (an Rz eigenstate cannot be a sigma_z = Rz.T invariant, and MOCCa
  stores Rz eigenstates); it carries no block structure - the classes
  select the stored half of the complex matrix elements, never a block.
- goodqn bit 3 (STy good) is not applied by the constructor: the mask is
  stored verbatim and read by consumers as the completeness criterion -
  with the bit set, an STy-tagged matrix is the entire (anti-)Hermitian
  operator; with the bit clear, it is one half of the (re, im) pair below
  and both halves must be built.
- Complex operators are represented as one stored real half plus the
  characteristic telling which part it is; in STy-conserving builds (all
  current configs) that half is the entire operator. The (re, im) pair
  model is retained for any future build that does not declare the time
  simplex: both halves are then separate classifiable-classes-verbatim
  matrices, sharing the symmetry word.
- No mat-vec and no mat-mat in the type (dropped 2026-11-19): no consumer
  applies an operator to a basis vector - MOCCa acts with fields on the
  mesh through the density route - and operator products are both
  consumer-less and unclassifiable. Partner machinery exists only as
  `dense()` and the consistency checks.
- Basis rotations (section 6a) are real orthogonal and sector-diagonal
  only; one bound method `rotate`, in place, metadata-free, touching the
  stored blocks alone. The structure checks (off-sector vanishing of the
  dense source array, orthonormality of every handed block, exact n/st
  match at rotation time) are constructor/rotate-time, loud. Operator
  products remain out of the type.
- TR-folded composites (section 2d) are a producer-side convention, not
  type machinery: the twin kernel is the honest kernel with T's image of
  the ket, its physical word is `IEOR(word, 1)`, it is ANTIHERMITIAN
  (asserted from MOCCa's fold conventions per family, the
  kernel-transpose oracle as arbiter), sty carries over per component,
  and only `$TR`-marked call sites construct it. The type never learns
  about time reversal; folded composites are ordinary classifiable
  matrices over the stored basis.

## 13. Step 3 detailed: generalised multipole matrix elements

Target (agreed 2026-11-17): move the single-particle matrix elements of the
electric multipole operators from `fission_MOI.f90` into `moments.f90`,
rebuild them on `sym_matrix`, generalise to all multipoles in natural
units, and keep the old routine in place as a testing oracle. Consumers
found in the tree: `fission_MOI` (collective inertia tensor), `functional`
(`verify_COM_motion`, Q_10 - a pushing-model diagnostic identified at
migration time as DEAD CODE, no call site anywhere; the earlier
"CompCOMCorrection" label was wrong, that routine never used `Qlm_spme`),
`fam` (`get_qpme_op`: multipole fields, z_com), `fam_testing`
([R, P] unit tests, uncalled). All migrated to the sym_matrix generator
by 2026-11-19 (step 5 plus the cleanup commit after it); the only
remaining `Qlm_spme` reference is the off-by-default oracle hook in
`fission_MOI`.

### 13.1 Semantics of the new generator

Operator: multiplication by `r^l Re Y_lm` or `r^l Im Y_lm`, evaluated on
the COM-shifted mesh (`meshx_shifted/...`, exactly as the current
`Qlm_spme`), in natural units `fm^l` (the `b^{l/2}` rescaling is dropped).
Loud input checks (`l <= MaxMoment`, `0 <= m <= l`);
symmetry-forbidden combinations return zero matrices silently - guarding
user input against forbidden moments is the business of Hephaestos' moment
list at the input layer (e.g. `read_inertia`), not of this routine.

Output: the stored half of the ME matrix in the HF basis, per the
operator's STy class, as one `sym_matrix`:

```fortran
! Im = .false.: r^l Re Y_lm, Hermitian, STY_EVEN - stored numbers are
!   Re <i|Q_lm|j>, init derives kind = KIND_SYMMETRIC
! Im = .true. : r^l Im Y_lm, Hermitian, STY_ODD - stored numbers are
!   Im <i|Q_lm|j>, init derives kind = KIND_ANTISYMMETRIC; in fm^l
subroutine Qlm_spme_matrix(l, m, Im, me_half)
! Dense legacy shape of the stored half, real(nwt,nwt), in fm^l
function Qlm_spme_fm(l, m, Imaginary) result(me)
```

- both operator halves are Hermitian and differ in their STy class:
  `init(HFBlocks_global, op, HERMITIAN, STY_EVEN/STY_ODD, goodqn = ...)`
  per the character table of section 2a (agreed 2026-11-19): STY_EVEN for
  r^l Re Y_lm, STY_ODD for r^l Im Y_lm, for every m. This supersedes the
  step-3 revert note of 2026-11-17: the within-block imaginary matrix
  elements are not artifacts - the STY_ODD class makes the imaginary part
  the physical content of the operator, and with MOCCa's linear
  antihermitian Rz the only antiunitary constraint on the stored basis is
  STy itself, which constrains the matrix elements to be purely imaginary
  but does not annul them.
- physical symmetry word (section 2): `op = mod(m,2) + 2*mod(l,2)` - bit 0
  flips q(Rz) for odd m, bit 1 is the parity of the operator, bit 2 is
  always clear (electric multipoles are isoscalar). Both operator halves
  (Re/Im Y_lm) carry the same word: MOCCa's Rz and P are linear, no
  antiunitary complement (section 8) arises for electric multipoles.
- chamber weighting on the reduced mesh (agreed 2026-11-19, the fill
  callback `Qlm_spme_fill_block` consumes the weights from its `params`):
  the stored mesh covers 2^NUMSYM chambers with `dv` carrying the
  multiplicity; the spwf components extend across the fold planes with
  the per-component reflection signs `sx/sy/sz` of `wavefunctions.f90`
  (the convention the derivative machinery uses), and the operator mesh
  function is of definite reflection character `O(-axis) = oc(axis) O`
  with `oc = +/- 1` estimated numerically from the mesh values. Every
  stored-point component-pair product of a row and a column state
  carries the exact 0-or-1 weight
  `w(k,k') = prod_(folded axis) (1 + s(k,i) s(k',j) oc) / 2`,
  which depends on nothing but the sector pair, the component pair and
  the operator: `Qlm_spme_build_weights` checks the sector-wise
  uniformity of the `sx/sy/sz` signs once, precomputes the entire table
  into the fill state `Qlm_spme_params` once per matrix, and the fill
  callback never touches a reflection character.
  For every block pair allowed by the XOR pairing all weights are 1 and
  the fill reduces to the plain dense products of the legacy routine -
  bit-identical numbers on the real half, so the oracle comparison
  against `Qlm_spme` survives unchanged; a weight 0 zeroes a component
  pair whose chamber contributions cancel through the reflection
  characters, the per-amplitude selection that the sector word cannot
  express.
- all matrix elements are direct mesh integrals of the stored basis
  functions - no time-reversal/signature-partner shortcuts, hence no
  partner-phase conventions to get wrong. Block fill: one dense product
  per spinor-component term per block (real half: `A = sum_k X_k^T
  (Q dv) X_k`; imaginary half: crossed products per spin,
  `A = sum_s [Re_s^T (Q dv) Im_s - Im_s^T (Q dv) Re_s]`), with `dv`
  folded into the weight mesh.
- helper in `block_matrix_module`: the partner-reconstructing dense
  expansion of a `sym_matrix` (refactored out of the `print_dense`
  core). It is the return shape behind `Qlm_spme_fm` and the reference
  shape for the elementwise oracle checks of the amended step 4
  (section 10).
- validation status: the real half is oracle-checked against the legacy
  `Qlm_spme` (agreement ~1e-15 relative on all requested moments, also
  with the weighted fill active). The imaginary half is covered by the
  full-box reconstruction oracle (amended step 4): selected off-diagonal
  elements recomputed by extending the spinors across the fold planes
  with the `sx/sy/sz` per-component reflection signs (the folding
  convention of the derivative machinery) and integrating
  `r^l Im Y_lm`, evaluated fresh at the full-box coordinates with the
  same COM shift, with plain cell volumes. The oracle shares only the
  trusted raw inputs with production - the chamber-weight table and the
  numeric `oc` character estimation are bypassed, so the weighting
  machinery itself is under test. Element selection: the
  largest-magnitude element plus a few random off-diagonal ones per
  stored block, in both orderings (i,j) and (j,i). Physics checks:
  purely imaginary oracle values for the STY_ODD operator,
  `<j|O|i> = <i|O|j>*`, and agreement with the stored blocks to ~1e-12
  relative (scalar summation vs the block dgemm products reorders the
  sum). Even-m imaginary halves carry within-sector antisymmetric
  content in every build, so NLO already exercises the type; NLO-T
  populates the odd-m signature-partner blocks additionally and NLO-P
  the parity-masked diagonal case. In the parity-broken build only the
  folded axes are reconstructed (the z-axis is stored unfolded).
  Landed status (2026-11-19): all three builds green - identity at
  reassociation level (max 5.3e-13 / 1.8e-13 / 6.8e-14 absolute on
  NLO / NLO-T / NLO-P), the opposite halves vanish to <= 6.1e-16 in
  every build and every (l, m, half), and the odd-m signature-partner
  imaginary blocks of the TR-broken build agree to ~1e-14.

### 13.2 The old routine stays, as oracle

The current `Qlm_spme` moves `verbatim` from `fission_MOI.f90` to
`moments.f90` (a pure move; behaviour preserved bit-for-bit) and remains
until the amended step 4 (the full-box oracle) validates the imaginary
half AND the FAM consumer `get_qpme_op` has migrated (step 5 / design
step 9) - the deletion gate was redefined 2026-11-19 with the dropping
of the mat-vec step. Amendment 2026-11-19, after step 5 landed: two
further legacy consumers were found in the tree and block deletion -
`functional.f90`'s `verify_COM_motion` (Q10; dead code, no call site
anywhere) and `fam_testing.f90`'s test_RP_commutator (uncalled test
code). Both were migrated in the cleanup commit after step 5, so the
gate is now purely a decision: either `Qlm_spme` and the `fission_MOI`
spme oracle go together (the oracle is its only remaining reference and
uses it as its reference matrix), or both stay. The self-checks in the moments unit net cover the
generator either way. It is the oracle for the agreements below; the
conversion identity is `old [b^{l/2}] = new [fm^l] * 100^(-l/2)`.

### 13.3 Sequencing (one commit per row)

| # | Commit | Contents | Regression net |
|---|---|---|---|
| 0 | park FAM anchor | build `exec/fam.NLO.reference.exe` from the un-refactored tree and park it | none needed |
| 1 | goodqn mask | optional `goodqn` in `init` (+ module header), `$GOODQN` injected by Hephaestos into `geninfo`, masked-init unit tests | `block_matrix_unit_tests`; compile NLO, NLO-T, NLO-P, NLO-TP |
| 2 | pure move | `Qlm_spme` verbatim to `moments.f90` (which gains explicit `use wavefunctions`), deprecation banner | `fission.sh`, `minimal.sh` golden: outputs bit-identical |
| 3 | generator | dense-expansion helper (with its unit tests), `Qlm_spme_matrix` (both operator halves as of 2026-11-19: real `STY_EVEN`/imaginary `STY_ODD`, chamber-weighted fill - formerly the reverted `me_im` of 2026-11-17), `Qlm_spme_fm`, and an off-by-default oracle hook in `fission_MOI`: compute both paths for every requested moment, print the max deviation, run symmetry self-checks. The imaginary half ships unvalidated (no reference path covers it yet) | unit tests; oracle runs NLO / NLO-T / NLO-P (real halves); all four configs compile |
| 4 | inertia switch | `calc_collective_inertia` takes `Qsp` from `Qlm_spme_fm * 100^(-l/2)`: the boundary rescale holds every presented number in `b` units | `fission.sh`, `minimal.sh` golden within tolerances |
| 5 | units flip | drop the internal rescale: the inertia chain (`Qsp`, `Ksum`, `calc_Q20`, the LAPACK inversion, the stored `collective_inertia`/`M1`/`M2`/`M3` arrays) runs in `fm` units end to end. The OUTPUT boundary restores the historical numbers for downstream users (amended 2026-11-19): `print_collective_inertia` and the inertia file (`write_inertias` in IO.f90) present every number in `b^{l/2}` units - the tensor element `I(i,j)` multiplied by `100^{(l_i+l_j)/2}`, the intermediate matrices `M_n(i,j)` divided by the same factor (they scale with `Q_i Q_j`, the inertia with `1/(Q_i Q_j)`). Printed and file values, their convention strings and the `fission.sh` references are therefore unchanged | `fission.sh` with unchanged references; `minimal.sh` golden |

`Ksum_Mij_HF/BCS/HFB`, `calc_Q20` and the LAPACK inversion are
deliberately untouched in this step beyond the source of `Qsp`.

### 13.4 Validation and risks

- Oracle agreement to machine precision modulo summation order (BLAS vs
  the old scalar loops - the golden tolerances of `golden_diff.sh` exist
  for exactly this reassociation drift).
- Parity-mask correctness: masked-init unit tests plus the NLO-P oracle
  run - the old routine's `$PBROKEN` diagonal-blocks path is the reference
  there.
- No available config breaks Rz today (all four SYMSTRINGs contain 'Rz'),
  yet the mask derivation is generic: an Rz-broken build clears bit 0 and
  odd-m operators act within the merged signature blocks - same rule, no
  new code.
- The units flip leaves every presented number unchanged (row 5, amended
  2026-11-19): the fm-based internal chain and the b-units output boundary
  agree with the legacy chain up to floating-point reassociation only, well
  inside the `fission.sh` tolerances; the `fission.sh` reference arithmetic
  of the original row 5 is no longer needed. The one internal consumer of
  the stored arrays, `calc_vibrational_correction` (`Tr M3^{-1} M2`), is
  scale-invariant and needs no rescale. `make.inc` stays frozen (the
  trust anchor is untouched).
- MPI out of scope (the sprint is serial-only): construction takes
  `HFBlocks_global` and the dense wrappers assume the whole basis on one
  rank, matching `fission_MOI`'s current inert-under-MPI behaviour
  (`N_inertia = 0`).
- Park the FAM reference BEFORE any `fam.f90` edit: commit 0 precedes
  commit 6.

### 13.5 Deferred (not this step)

- odd-m and Im-operator FAM fields end to end,
- inertia for Im-operator moments (`Ksum`/`calc_Q20` generalisation),
- the optional `$TR` short-circuit for the antisymmetric half,
- deleting the old `Qlm_spme`.

## 14. Angular momentum operators: the angular_momentum module

Target (agreed 2026-11-19): a new module `src_orig/angular_momentum.f90`
hosting generators of the angular-momentum one-body matrices on
sym_matrix - the first producer family since the multipoles of section
13, and the supplier for the rotation consumers of section 6a (the
canonical-basis J-matrix chains of the momentsofinertia/Belyaev
machinery are the first production consumer).

Base conventions recorded on MOCCa's side (user-confirmed 2026-11-19):
- MOCCa's time reversal is FERMIONIC: T**2 = -1 on the spwfs.
- The angmom kernels of wavefunctions.f90 are the authoritative ME
  implementations; the generators below wrap them verbatim as values
  per stored pair, so the legacy dense loops of
  `calcJ2andBelyaev_HFB` are themselves the A/B oracle (the step-3
  precedent: generator vs moved legacy path).

### 14.1 The five operators

| object | word | herm | sty | kind | stored numbers | kernel | measured mirror (14.3) |
|---|---|---|---|---|---|---|---|
| `J_x` | 1 | HERMITIAN | STY_EVEN | +1 symmetric | `Re <i,J_x,j>` | `angmom_x_real` | symmetric, 7.8e-15 (NLO-T) |
| `J_y` | 1 | HERMITIAN | STY_ODD | -1 antisym. | `Im <i,J_y,j>` | `angmom_y_imag` | antisymmetric, 2.1e-15 (NLO-T) |
| `J_z` | 0 | HERMITIAN | STY_EVEN | +1 symmetric | `Re <i,J_z,j>` | `angmom_z_real` | symmetric, 2.5e-14 (NLO) |
| `J_x.T` | 0 | HERMITIAN | STY_EVEN | +1 symmetric | `Re <i,J_x T,j>` | `angmom_xt_real` | symmetric, 7.5e-15 (NLO-T) |
| `J_y.T` | 0 | ANTIHERMITIAN* | STY_ODD | +1 symmetric | `Im <i,J_y T,j>` | `angmom_yt_imag` | symmetric, 2.1e-15 (NLO-T) |

MEASURED (2026-11-19, the kernel-transpose oracle of 14.3, O16/SLy4
unit state, maxiter=0): every tagged kind above agrees with the
kernel-measured mirror in every build with stored pairs; the table's
kinds are test assertions that PASSED, not assumptions. In TR-conserving
builds the honest word-1 objects are structural zeros (the odd-m
precedent): 0 stored blocks on NLO for J_x/J_y, exactly as predicted.
NLO-P (frozen -O3): J_z and both folds symmetric at 8.9e-16 / 1.1e-15 /
4.8e-16, on the merged-parity diagonal blocks - the parity-masked case
covered; the word-1 honest objects are structural zeros there too.

- Honest trio: `J_x`/`J_y` flip q(Rz) (word 1), `J_z` conserves it
  (word 0); all three are parity-even and isoscalar (bits 1, 2 clear);
  J_y is the STY_ODD member per the character table of section 2a. In
  TR-conserving builds the word-1 partners of every stored sector are
  empty, so honest `J_x`/`J_y` come out structurally zero there - the
  odd-m multipole precedent, no special case.
- Folded pair: with only q(Rz) = +i stored, the TR-conserving consumer
  (`calcJ2andBelyaev_HFB`, `$TR` branch) needs the time-reversal-folded
  objects `<i, J_mu T, j>` the `t` kernels implement. Both are word 0
  (T and J_mu each flip the signature, the composite preserves it) and
  live within the stored +i sectors. The composite is HERMITIAN as an
  operator (user-accepted 2026-11-19): with MOCCa's FERMIONIC time
  reversal, T**2 = -1, the chain
  ` (J.T)^dagger = -T**2 (J.T) = +(J.T)` holds - and the ORACLE
  confirms it operationally: both folded kernels measure SYMMETRIC
  storage (kind +1). The section-2a (herm,sty)->kind machinery does NOT
  transfer to the antilinear composite - the sty-odd member `J_y.T`
  combines a symmetric kind with sty = STY_ODD, so its `herm` input is
  the kind-completing sign `herm = kind * sty` (marked * in the table):
  a storage label, not the physical adjoint class. `J_x.T` carries the
  plain HERMITIAN label. The empirical arbiter of every kind sign is
  the kernel-transpose oracle of 14.3.
  NOTE (open, user to resolve): the labelling convention itself -
  expressing the measured folded kinds through the overloaded `herm`
  input - is PROVISIONAL. Whether it stands permanently, or the type
  grows an honest folded-composite constructor that does not reuse
  `herm` for a non-adjoint sign, is Wouter's decision; either way the
  oracle remains the arbiter of every kind sign, never re-derivation.
- `J_z` needs no folded twin: the `$TR` branch consumes honest
  `angmom_z_real` within the half (J_z preserves the grading).

### 14.2 Convention audit (outcome 2026-11-19: the fill is NOT bugged)

The earlier suspicion that the `$TR` fill `jx(jj,ii) = jx(ii,jj)`
(momentsofinertia.f90:832) was a sign bug (the asserted ANTIHERMITIAN x
STY_EVEN classes of K_x demanding an antisymmetric fill) was REFUTED by
measurement and A/B:

- The kernel-transpose oracle measures BOTH folded kernels SYMMETRIC
  (14.1 table) - the symmetric fills of `calcJ2andBelyaev_HFB` are
  correct, and the composite is Hermitian as an operator (14.1, the
  user-accepted chain with fermionic T).
- An attempted antisymmetric fill was A/B-tested with `symmetries.sh`
  (Mg24/SLy4/HFB constrained, the four NLO-family executables): in the
  pre-change run the routine's two internal definitions of
  `<J^2> - <J>^2` (pairing cutoff vs collective) agree on every
  component and the constrained Belyaev run converges to
  I_B(X) ~ 0.9295; with the antisymmetric fill alone EXACTLY the
  x-consuming rows move - the two `<J^2>-<J>^2` definitions split on X
  (pre-change 1.4880 both ways; broken 6.8473 vs 2.2486) and the Belyaev
  X rows move 0.9295 -> 3.2213, with the Y/Z rows
  bit-identical. The fill and the `$TR` consumer algebra are a coherent
  pair; the source was reverted, no production line changed.

The measuring oracle printout in the unit arm is the permanent arbiter:
any future claim about folded kinds re-runs it, never re-derives.

### 14.2b The angmom opposite-half purity prints (RESOLVED 2026-11-19)

The claim in the code, verbatim - the block comment above the
j-matrix fill of `calcJ2andBelyaev_HFB` (momentsofinertia.f90, two
branches, `$NTR` lines 774-775 and `$TR` line 792):

    jz_ij is a always a real number, and only non-zero if i and j are
    in the same 'half' of the basis (...)
    jz_ij and jx_ij are always real, and jy_ij is always imaginary.

This is not a free-standing comment: it is the COMPLETENESS premise of
the storage - section 2a's corollary for J_z (HERMITIAN x STY_EVEN):
with goodqn bit 3 set (STy conserved, every config in the tree) every
matrix element of J_z over the stored basis is purely real, the stored
`Re <i,J_z,j>` half is the complete operator, and `Im <i,J_z,j>` is
identically zero.

The measurement (angular_momentum_testing purity print, 2026-11-19,
O16/SLy4 unit state, maxiter=0): `angmom_z_imag(psi_i, psi_j, dwf_j)`
under the SAME calling convention that the arm validates for
`angmom_z_real` (stored spinors, HFdpsi or fresh Derive_* gradients)
returns values of the SAME ORDER as the real half, in every build:

    NLO   : opposite-half max 1.4146   vs stored-half max 1.9146
    NLO-T : opposite-half max 0.8594   vs stored-half max 1.6500
    NLO-P : opposite-half max 0.8064   vs stored-half max 1.0000

The tension is sharp: if `angmom_z_imag` is what its name says - the
imaginary half of `< wf2 | J_z | wf1 >`, and a structural reading of
its spin and orbital terms against the complex product supports that
reading - then the measured O(1) contradicts the identically-zero
claim and the section 2a completeness argument behind the whole
stored-half model. x and y are in the same position (`angmom_x_imag`,
`angmom_y_real`, also O(1), also production-dead: the arm is their
first caller). NOTHING in production is affected whatever the answer:
the stored real halves are mirror-oracle-clean, and the dead kernels
have no consumer. The arm PRINTS the magnitudes, it does not judge
them.

Two candidate resolutions, deliberately not adjudicated here:

(a) REDUCED-MESH ARTIFACT (suspected, UNVERIFIED): the purity
    evaluation integrates one chamber of the folded mesh with the dv
    multiplicity, exactly like the generator fills - valid only for
    component-pair products UNIFORM across chambers. The
    opposite-half products are precisely the cross-parity terms the
    multipole chamber-weight table (13.1) assigns weight ZERO: they
    cancel across the fold planes in the full box, so the true
    Im <i|J_z|j> is ~0 as the comment claims, and the naive
    single-chamber integral produces the O(1) residue. On this
    reading the premise stands, the kernel is worthless as called on
    the reduced mesh, and no bug exists - only a trap.

(b) GENUINE KERNEL BUG or packing-convention mismatch in the dead
    kernel: the premise then needs re-examination wherever the
    comment feeds consumer algebra.

The same-halves reading holds for the opposite-half kernels of the
honest x/y components (`angmom_x_imag`, `angmom_y_real`): also
production-unreferenced, also O(1) on the reduced mesh, same two
candidates, same deciding test.

RESOLUTION (2026-11-19, audit + deciding test both executed):

A term-by-term rederivation of all three dead kernels against the
packing conventions pinned by the validated kernels themselves
(spin_z_real, angmom_z_real, angmom_z_quad's action convention,
angmom_x_real, and angmom_y_imag - the production-validated imaginary
half) found:

- `angmom_z_imag`: CORRECT. Every spin and orbital term is exactly
  Im < wf2 | J_z | wf1 >. No bug.
- `angmom_y_real`: CORRECT. The honest Re half of the STY_ODD J_y.
- `angmom_x_imag` (and its clone `spin_x_imag`): the spin part's
  down-bra/up-ket block was mistranscribed - it carried the sigma_y
  pairings (wf2(3)wf1(1), wf2(4)wf1(2): the pairing angmom_y_imag
  correctly uses for y) dressed with x-style antisymmetric signs,
  instead of the honest
      Im[ conj(d_down) c_up ] = wf2(3) wf1(2) - wf2(4) wf1(1),
  thereby also breaking the pointwise mirror-antisymmetry of the
  imaginary half. FIXED; the up-block, the entire orbital part and
  angmom_x_real were and are correct. Both kernels are now
  banner-annotated as caller-less (the arm's purity print is
  angmom_x_imag's only reference in the tree; spin_x_imag has none).

The deciding test IS BUILT: the check_fullbox / fullbox_angmom_me
pair in angular_momentum_testing - the 13.1 step-4 full-box
reconstruction extended to the derivative-and-spin family. Spinors
are extended across the fold planes with the sx/sy/sz per-component
signs, gradient fields with the additional chain-rule sign (-1)^flip
per differentiated axis, the operator's coordinates are the mirrored
COM-shifted mesh values per chamber, the spin part follows the Pauli
structure, cell volume dv/2^nch. Per sampled pair of every stored
block it JUDGES two facts: identity (the main half is
chamber-uniform: the full-box quadrature equals the kernel evaluation
on the single stored chamber) and the STy zero (the opposite half
vanishes over the full box). The oracle's chamber algebra reproduces
the validated kernels term-by-term on the unflipped chamber, so the
identity judgment validates the oracle itself - the trap of this
section cannot be re-set silently.

VERDICT (2026-11-19, O16/SLy4 unit state, maxiter=0, the moments-suite
scenario; x/y honest are structural zeros in NLO and NLO-P - the
odd-m precedent - so their pairs live in NLO-T only):

    build  comp  pairs  identity max  opposite-half max
    NLO    z     260    1.6e-14       2.1e-17
    NLO-T  x      68    2.3e-14       1.4e-17
    NLO-T  y      68    1.7e-14       2.1e-17
    NLO-T  z     136    1.6e-14       1.0e-17
    NLO-P  z      72    7.3e-15       2.1e-17   (merged-parity blocks)

Candidate (a) CONFIRMED on every measurable case: the full-box
Im <i|J_z|j> vanishes at the 1e-17 level in every build while the
same kernel returns O(1) on the single stored chamber - the purity
prints measured chamber residues all along. The mechanism is
per-term and explains both halves at once: on the stored +i sign
patterns every opposite-half term is fold-odd (the spin term
wf2(1)wf1(2) carries x-character sx(1)sx(2) = -1; the orbital term
wf2(1) y d_x wf1(1) carries sx(1) x (-sx(1)) = -1 - a derivative
flips the fold sign) while every main-half term is fold-even - which
is exactly why the main half is valid on the reduced mesh without a
weight table and the opposite half is not (the 13.1 chamber-weight
zero, term by term). The premise "jz_ij is always a real number" is
JUDGED GREEN, never assumed; the purity print is relabeled as
residues and the arm fails loudly if either judgment breaks. The
nabla analogue of section 15.6 remains open - the derivative family
has no full-box oracle yet, and this machinery is its ready template.
No production line changed at any point: the fixed kernels remain
caller-less.

### 14.3 Module shape, marker policy, validation

- `src_orig/angular_momentum.f90`: plain Fortran, like
  block_matrix_module free of `$`-markers; takes the angmom kernels
  from wavefunctions; generators built as
  `new_sym_matrix`-style constructors
  (`init(HFBlocks_global, word, herm, sty, goodqn)` per the 14.1 table;
  the sugar itself has since been dropped, section 5)
  with block_filler callbacks looping the kernels over the stored
  pairs of each block (the kernels are per-pair pure functions; no
  dgemm restructure is attempted - performance equals the legacy
  loops). Wiring: MF_SRC/FAM_SRC in the Makefile and
  `no_process_files` in src_heph/preprocess.py.
- Marker policy (agreed 2026-11-19): the folded generators EXIST in
  every build and the module stays marker-free; their CALLERS gate -
  only `$TR`-marked consumers construct the folded pair. Dead code in
  the TR-broken builds by construction, not by preprocessing.
- Validation:
  *) kernel-transpose oracle: per component, the kinds are settled
     numerically - `K(wf_i, wf_j)` vs `K(wf_j, wf_i)` measured by the
     kernels themselves on real states fixes the relation A_ij =
     (+/-) A_ji per component. The kind signs of 14.1 become test
     assertions, not assumptions. Runs in every build (kernel
     property); the folded objects' CONSUMER chain is TR-only, but
     the sign checks do not need the TR data model.
  *) legacy A/B: the existing dense loops of
     `calcJ2andBelyaev_HFB` remain in place as the oracle until the
     consumer migration; the moved path must be bit-identical for
     `J_z` within any block and for the `$NTR` honest fills.
  *) partner-consistency checks: `A(s',s) = kind * A(s,s')^T` per
     stored block, both kinds, all three builds; XOR block-coverage
     against the words of 14.1.
  *) full-box reconstruction oracle (2026-11-19, the deciding test of
     14.2b, now part of the arm): check_fullbox judges, per sampled
     pair of every stored block of the honest trio, the chamber-
     uniformity identity of the main half and the STy zero of the
     opposite half; its chamber algebra reproduces the kernels on
     the unflipped chamber, so the identity judgment self-validates
     the oracle.
  *) unit arm placement: beside the moments net (moments_testing or
     an angular_momentum_testing module), armed by the `unit_test` flag
     of &mesh/ per the §13 test infrastructure.

### 14.4 Sequencing (one commit per row)

| # | Commit | Contents | Regression net |
|---|---|---|---|
| 1 | ~~K_x fill fix~~ | CANCELLED by the 14.2 audit: the symmetric fill is correct (measured); no source change, no commit |
| 2 | generators | LANDED 2026-11-19: angular_momentum.f90 (the five generators with the measured 14.1 classes), build wiring (MF_SRC/FAM_SRC, no_process_files), the unit arm `angular_momentum_testing.f90` (kernel-transpose oracle, partner checks, structural-zero and purity prints, legacy A/B check) armed from Run_MOCCa beside run_moments_tests | arm green on NLO / NLO-T / NLO-P on the frozen make.inc (2026-11-19; debug-profile NLO-P is unusable per the known pre-existing coulomb-folding OOB, AGENTS.md open item); golden minimal NLO 4849/4849 bit-identical (the arm runs only when unit_test is armed); block_matrix net untouched; mf+fam of NLO/NLO-T/NLO-P and mf of NLO-TP compile clean |
| 3 | consumer migration | `calcJ2andBelyaev_HFB` builds its j-matrices from the generators (`$NTR`: honest trio; `$TR`: honest J_z plus the folded pair from `$TR`-marked call sites); the cantransfo rotations go through sym_rotation (section 6a) in a following row | fam unit net unaffected here; inertia smoke + minimal golden; fission.sh as the heavy gate at the end of the rotation chain |

Consumers further downstream (Bogoliubov J20/J11 images, qp machinery)
belong to the qp-bookkeeping discussion deferred from section 6a and
stay out of this module's scope.

## 15. Gradient operators: the linear_momentum module

Target (agreed 2026-11-19, design stage): the single-particle matrix
elements of the gradient operator family migrate from densities.f90's
`CompNablaMelements` to sym_matrix generators in a new marker-free
module `src_orig/linear_momentum.f90` - separate from the
angular_momentum module for now, to be merged with it once both exist
(user decision; both wrap wavefunctions.f90 kernels and both carry the
folded pair of section 2d). The legacy routine stays in place, verbatim,
as the A/B oracle until a deletion-gate decision, the section 13.2
precedent. MPI is DROPPED in the generators (user decision): the sprint
is serial-only, the type is single-node by design, and the routine's own
TODO already listed "removing the MPI, as this is functionality that is
not currently used". The legacy routine keeps its MPI machinery for as
long as it exists, so the deletion gate is the only place MPI users lose
anything.

### 15.1 What CompNablaMelements actually computes (audit 2026-11-19)

The function advertises `NablaMElements(3,2,nwt,nwt)`, but only three of
the six `(mu, half)` slots are ever filled, each with the partner
triangle written by a sign:

    (1,1) : Re <i|\nabla_x|j>  , partner fill MINUS (antisymmetric)
    (2,2) : Im <i|\nabla_y|j>  , partner fill PLUS  (symmetric)
    (3,1) : Re <i|\nabla_z|j>  , partner fill MINUS (antisymmetric)

one stored half per component, exactly as the classes of section 2a
demand; the other three slots stay zero in every build (goodqn bit 3 is
set in all current configs). The `$TR`/`$NTR`/`$PBROKEN` branch maze is
the masked-XOR pairing of section 2c hand-rolled:

- honest `\nabla_z`, word 2 (P-odd, Rz-even): block pairs 1-3 / 2-4 in
  every build; the `$PBROKEN` within-block loops are the parity bit
  cleared in the mask, the diagonal case of the same word.
- honest `\nabla_x`/`\nabla_y`, word 3 (P-odd, Rz-flipping): `$NTR`
  loops over pairs 1-4 / 2-3; the `$PBROKEN $NTR` 1-2 loop is word 3
  masked to 1 by the parity-broken mask.
- the `$TR` hidden_TR loops are NOT honest operators at all: they
  compute the section 2d folded pair `\nabla_x.T`/`\nabla_y.T` (word 2)
  over the same block pairs 1-3 / 2-4 as `\nabla_z`, plus the `$PBROKEN`
  diagonals.

The load-bearing wart: the SAME array slots hold honest matrix elements
in `$NTR` builds and folded ones in `$TR` builds. The migration
dissolves this semantic overload into named generators, and consumers
request honest-vs-folded by construction, per section 2d. No fill-bug
suspicion arises here (and the 14.2 audit dissolved the angmom one
anyway): the sign of every partner fill matches the kind demanded by the
classes below - but per the amended section 2d the folded kinds are
MEASURED by their own kernel-transpose oracle before they arm, not read
off the class rule (the angmom precedent: the asserted classes were
wrong while the fills were right).

### 15.2 The five generators

| object | word | herm | sty | kind | tau (15.8) | stored numbers |
|---|---|---|---|---|---|---|
| `\nabla_x` | 3 | ANTIHERMITIAN | STY_EVEN | -1 antisym. | +1 (T-even) | `Re <i,\nabla_x,j>` |
| `\nabla_y` | 3 | ANTIHERMITIAN | STY_ODD | +1 sym. | +1 (T-even) | `Im <i,\nabla_y,j>` |
| `\nabla_z` | 2 | ANTIHERMITIAN | STY_EVEN | -1 antisym. | +1 (T-even) | `Re <i,\nabla_z,j>` |
| `\nabla_x.T` | 2 | HERMITIAN | STY_EVEN | +1 sym. | +1 (inherited) | `Re <i,\nabla_x T,j>` |
| `\nabla_y.T` | 2 | ANTIHERMITIAN | STY_ODD | +1 sym. | +1 (inherited) | `Im <i,\nabla_y T,j>` |

AMENDED 2026-11-19 by the mirror oracle at the first unit-arm run: the
folded x row drafted here as ANTIHERMITIAN x STY_EVEN (kind -1, after
the legacy in-loop partner fill) was measured SYMMETRIC - max
|a_ij - a_ji| ~ 5e-15 against max |a_ij + a_ji| ~ 0.7 - and the tags
follow the arbiter (section 2d's own rule): folded `\nabla_x.T` is
HERMITIAN x STY_EVEN, kind +1. Exactly the J family's story: the
legacy's transposed-position fill carries a sign bug for the
composite (the 14.2-style wart; class tags drafted from it inherit
the bug), invisible to every consumer because CompOMCorrection
contracts squares. The folded y row was confirmed as drafted. Both
folded objects now carry the same per-component (herm, sty) shape as
the landed angular-momentum family's K_x / K_y.

- The tau column is the time parity of section 15.8's mechanism 1:
  the gradient as written is time-EVEN (T antilinear means the i of
  p = -i(nabla) carries the momentum's time-oddness), user-confirmed
  2026-11-19; the composites inherit it via
  T(OT)T^{-1} = tau(O)(OT). The oracle check (one image-image block
  computed both ways) remains the standing arbiter before the algebra
  path arms - first measurement 2026-11-19: agrees at 4e-15.

- Classes user-confirmed 2026-11-19: all five antihermitian; the
  per-component characteristics follow the vector pattern of section
  14.1 with the parity bit set (`\nabla` is a P-odd, isoscalar vector;
  `\nabla_x`/`\nabla_y` additionally flip q(Rz), `\nabla_z` conserves
  it). The twins obey the section 2d convention; the kind signs read
  off the legacy fills (antisymmetric Re slots, symmetric Im slot) and
  are settled as test assertions by the kernel-transpose oracle.
- `\nabla_z` has no folded twin with content: honest word 2 -> twin word
  3 -> partner sectors empty in a `$TR` build (every TR-conserving
  build stores the +i half, parity-broken included), the section 2d
  collapse. The honest object serves every such build, like `J_z`.
- Suggested interface (one function per operator, the 14.3 shape):
  `init(HFBlocks_global, word, ANTIHERMITIAN, sty, goodqn)` with a
  `block_filler` callback; kernels are `Transfer_psi`/`Transfer_derpsi
  (direction, TR flag, basis)` - the folded generators pass the TR flag
  the legacy `$TR` loops pass, TimeReverse applied to the transferred
  derivative, verbatim. Like section 14 no dgemm restructure of the
  kernel arithmetic is required, but the block fill MAY materialise
  per-block component arrays and use dense products, the section 13.1
  pattern: the `\nabla_z`/`\nabla_x` Re halves are
  `A = sum_c Psi_c^T (dv) Der_c`; the `\nabla_y` Im half is the crossed
  four-product-per-spin shape of the 13.1 imaginary half,
  `A = Psi_1^T D Der_2 - Psi_2^T D Der_1 - Psi_3^T D Der_4 +
  Psi_4^T D Der_3` with D the dv weight. No chamber-weight table: the
  derivative machinery owns the fold planes and the legacy fill is
  plain `dv*sum(...)`, preserved exactly.
- The generators take the same optional basis argument as the legacy
  routine (`'HF'`/`'CAN'`/`'DEN'`, default `'DEN'`) and fill from
  `Transfer_psi(basis=...)/Transfer_derpsi(basis=...)`. This is
  deliberately NOT a sym_rotation consumer yet - the basis strategy is
  section 15.4.
- Wiring the module per section 14.3: MF_SRC/FAM_SRC in the Makefile,
  `no_process_files` in src_heph/preprocess.py; unit arm beside the
  moments net (moments_testing or a linear_momentum_testing module),
  armed by the `unit_test` flag of &mesh/.

### 15.3 Consumers (audit 2026-11-19)

Four call sites for `CompNablaMelements`, two live:

1. `fam.f90` `get_qpme_op`, `operator_type='zmomentum'` (live): wants
   the full complex `\nabla_z` sp matrix in the HF basis, assembled as
   `f_spme = dcmplx(NablaMElements(3,2,:,:), -NablaMElements(3,1,:,:))`.
   With `\nabla_z` STY_EVEN and goodqn bit 3 set, the Imaginary slot is
   identically zero in every build - so
   `f_spme = -i * dense(\nabla_z stored half)`, complex with zero real
   part. This is the FIRST live use of the section 2a/2c completeness
   claim outside the multipole family: the -i and the hbar convention of
   P_z stay at the call site, the step-5 precedent. Runtime caveat,
   inherited from the operator's word 2: a P-odd response field is not
   representable in parity-conserving FAM builds and is expected to NaN
   exactly as the zcom/l=1 fields do (the parked-anchor criterion); the
   working gate is the NLO-P fam build, per the zcom precedent.
2. `functional.f90` `CompCOMCorrection`, the `COM2body = 1` path
   (live): contracts squared nabla matrix elements against
   `rho_can`/`kappa_can` in the DEN basis. `$TR`: the folded pair plus
   honest `\nabla_z` (engineering factors 1/4 and the factor 2 for
   the TR half-sum kept verbatim); `$NTR`: the honest trio plus the
   `conjugp` double sums (nabla_{i j} nabla_{i bar j bar} against
   `kappa_can`), with the existing stop guard for TR-broken HF.
   All contractions are elementwise and remain consumer-side (section
   6); the migration replaces the SOURCE of the numbers, mainly so
   that the `$TR` branch's request reads as the folded pair by
   construction instead of "slot (1,1) semantics".
3&4. `fam_testing.f90` `test_RP_commutator` and `functional.f90`
   `verify_COM_motion`: dead code (no call sites), migrated alongside
   as in the 4ca0f181 cleanup. verify_COM_motion already does a
   by-hand `cantransfo` rotation of a dense nabla matrix in its HFB
   branch - its migration can foreshadow the section 15.4 audit.

### 15.4 Basis strategy (audit, then rotate - user decision 2026-11-19)

The generator phase fills per requested basis, verbatim legacy
semantics, so every consumer migrates safely without a semantics bet.
Separately, the audit question: is the DEN basis verbatim the HF basis
rotated by `cantransfo` (sector-diagonal per section 6a, `Canonical`
diagonalises per MOCCa block)? If yes, CompCOMCorrection switches to
the HF generators plus sym_rotation(cantransfo) - the FIRST production
consumer of the section 6a rotations, the natural showcase of `rotate`
's metadata preservation. Until the audit says yes, per-basis filling
is the shipped behaviour; `Rsq_spme`-style direct oracles decide the
identity, not an assumption.

### 15.5 Sequencing (one commit per row)

Updated 2026-11-19 for the implicit-retrieval contract of 15.8 and the
tau table of 15.2.

| # | Commit | Contents | Regression net |
|---|---|---|---|
| 0 | scenarios | define `COM2body = 1` scratch states: one NLO (the folded path), one NLO-T HFB state (the honest + conjugp path); run each with the pre-migration binary and keep the logs as reference numbers | none |
| 1 | generators + accessor + unit arm | `linear_momentum.f90`: the three honest generators, the retrieval accessor of 15.8 (stored blocks served as stored, TT blocks by the `TAU_NABLAFAM = +1` sign flip, cross blocks by the TR-flag kernels; the folded twins reachable as its eager packaging), wiring (Makefile + `no_process_files`), the unit arm (kernel-transpose oracle, partner-consistency checks, XOR block-coverage vs the 15.2 words, the tau both-ways image-image check, doubled-basis full traversal against legacy `CompNablaMelements`); serial only (15.6) | the new unit arm on NLO / NLO-T / NLO-P; block_matrix net untouched; mf builds of all four configs compile (fam builds compile the module too - functional.f90 links it); legacy routine untouched |
| 2 | COM correction consumer | `CompCOMCorrection` from the generators: `$TR` = fold content retrieved via the accessor (twin packaging) + honest `\nabla_z`, the fac/4 and 2x factors kept verbatim at first (15.8's connection note governs any later re-derivation); `$NTR` = honest trio + conjugp sums; basis 'DEN' per-basis fill | `COM2body = 1` scratch A/B vs the row-0 reference logs (reassociation-level drift expected); minimal golden bit-identical (`COM2body` defaults to 0) |
| 3 | dead-code consumer | `verify_COM_motion` (mf side) migrated off `CompNablaMelements` | compiles; source-level |
| 4 | deletion gate | the CompNablaMelements deletion decision (with its MPI machinery - they go together) | decision recorded |

Deferred by user direction (2026-11-19): everything FAM-side. The
zmomentum consumer of `get_qpme_op` keeps its legacy
`CompNablaMelements('HF')` call; no fam binaries are parked for this
step; `test_RP_commutator` (fam_testing) stays on the legacy routine
with the deferred fam consumers of 13.5. The zmomentum row as
originally drafted (binary parking, `-i * dense(\nabla_z)`, the
parity-odd NaN criterion) re-attaches to the fam phase when it comes.

Not in this table, deliberately: the DEN-vs-cantransfo audit and the
possible sym_rotation switch of CompCOMCorrection (15.4) - a separate
follow-up step, not a consumer-migration row; and the merge of
linear_momentum with angular_momentum (user direction), which waits
until both modules exist and brings the J family's own tau entries
with it.

STATUS (2026-11-19, evening): rows 1-3 all LANDED. Row 2's A/B
verdicts and the fac/4-2x open question are recorded in 15.6. Row 3
(verify_COM_motion, dead code) also migrated its one
compNablaMElements call to Nabla_half_fm(AX_Z, .false., 'DEN') -
identical content to the legacy slot (3,1) it consumed; gate is
source-level plus the NLO/NLO-T mf compiles. The remaining
CompNablaMelements references are exactly the deferred FAM pair
(fam.f90 zmomentum, fam_testing.f90 test_RP_commutator), the
definition itself, and the unit arm's A/B oracle - i.e. row 4 is now
purely the deletion decision, pending the FAM phase.

### 15.6 Validation and risks

- No golden path exercises this operator family today: `COM2body`
  defaults to 0, no integration script sets `zmomentum`. Net is built
  deliberately (row 0); the legacy routine is the oracle
  (`CompNablaMelements` stays in place verbatim, like `Qlm_spme`).
- Reassociation: legacy fills are scalar per-pair sums, the generator's
  block fill is dense products - expect ~1e-13-level drift of the
  golden-tolerance philosophy, judged per the oracle hook's max
  deviation, not by bit-identity.
- Structural coverage: NLO exercises the folded path and honest
  `\nabla_z`; NLO-T the honest trio with `$PBROKEN` masked lists only
  via config differences; NLO-P the parity-bit-masked diagonal case.
  The honest path needs pairing (HF + `$NTR` stopped by construction),
  so the `$NTR` A/B is an HFB/BCS scratch state.
- Main rewire risk: the `$TR` branch of CompCOMCorrection must receive
  the FOLDED generators, not honest word-3 objects (which are
  structurally zero there anyway - a wrong choice fails loudly as
  zeros, not silently as wrong numbers).
- The DEN-basis audit precedes any sym_rotation consumer (section
  15.4); per-basis filling stays the shipped safety until then.
- MPI behaviour is unchanged for as long as the legacy routine exists;
  dropping it at the deletion gate is a user decision, recorded as such.
- NLO-P bookkeeping (user-confirmed 2026-11-19): a parity-broken build
  does NOT store both signature sectors - the parity blocks MERGE into
  runtime-populated blocks 1, 5 (slots 2, 6 are permitted by the grid
  but carry zero count while T is conserved, so the legacy `$PBROKEN`
  block-2 loops are no-ops). Every stored state is q(Rz) = +i: the
  generators see the same +i-only basis as in NLO, and the
  constructor's `s = block - 1` encoding stands unchanged. The
  occupation row of 2b is amended accordingly.
- OPEN FINDING, measured 2026-11-19 at the first unit-arm run,
  deliberately NOT judged: the "purity" identity of the honest
  generators - opposite half (eta-stripped) blocks vanishing, i.e.
  <STy-|nabla|STy+> ~ 0 - FAILS for the derivative family. Magnitudes
  O(0.2) against a O(0.35) element scale, on both the maxiter=0 and
  a converged state, on the crossed-component pairs. The basis is not
  at fault (the multipole full-box oracle proves STy cleanliness at
  1e-16 on this very machinery), so the suspicion points at how the
  derivative kernels treat the crossed component pairs - the section
  13.1 chamber-weight territory. The angular-momentum family's
  analogue was the same signature and is now ADJUDICATED (14.2b:
  reduced-mesh artifact confirmed, one dead-kernel pairing bug
  fixed); by mechanism the artifact reading is favoured here too -
  the crossed-component products are the chamber-cancelling ones -
  but the derivative finding stays deliberately UNJUDGED until the
  derivative family gets its own full-box reconstruction: the
  angmom arm's check_fullbox is the ready template. The unit arm
  MEASURES and reports both halves; judgment (and any fix) is
  deferred until the derivative-kernel chamber audit - which is a
  prerequisite for
  nothing in this section, since every consumer contracts real
  densities against these elements and the legacy agrees to 3e-15
  exactly as filled.
- MIGRATION A/B VERDICTS (2026-11-19, row 2 of the 15.5 sequencing):
  the CompCOMCorrection consumer fed from the linear_momentum
  generators (dense twin packaging of 15.8) passes every gate against
  the pre-migration parked binary - the fac/4 and 2x legacy factors
  survive VERBATIM with bit-similar results: $TR (NLO, HF and BCS,
  ph and pp sums) - energies bit-identical, only last-digit DE/DR
  drift at the 1e-16 level; $NTR (cranked Mg24 HFB, honest trio +
  conjugp sums) - body of the logs BIT-IDENTICAL. This A/B validates
  the migration's numerics, NOT yet the tau bookkeeping of 15.8's
  connection note: whether the fac/4 and 2x factors are mechanistically
  the mechanism-1 T-pair mirroring remains unadjudicated (no
  dedicated unit-arm test of the factors was built; any later
  re-derivation is governed by that note).
- Doubled-basis oracle (proposed 2026-11-19 with the section 2d discussion, user reading: the folded composites are the honest matrix
  elements of O between stored rows and the time reverses of the stored
  columns): materialise the T-images as explicit basis states, hand
  `init` a doubled HFBlocks (+i entries copied to their -i partner
  slots; pure XOR, no new type machinery), build the honest O over the
  doubled basis by its honest kernels, and check its stored-stored
  blocks against the honest generator and its cross blocks against the
  twin. This settles the twin class table empirically - the strongest
  arbiter available, subsuming the kernel-transpose checks - and makes
  the identity between the legacy hidden fold
  (`TimeReverse(derivative)`) and the literal matrix element with the
  T-image state (`derivative(TimeReverse)`) a test rather than an
  interpretation.

### 15.7 Deferred (not this step)

- the sym_rotation switch of CompCOMCorrection (after the 15.4 audit),
- folding the multipole family (a section 13.5 decision, if ever),
- deleting `CompNablaMelements` (and its MPI machinery with it).

### 15.8 The complete matrix, implicit (user direction 2026-11-19)

Complete-fill (eager: the generators always materialise every block of
the doubled-basis matrix - was considered as the previous wording of
this section) is rejected in favour of the intermediate the user chose:
the additional matrix elements - those involving the T-image states -
are carried IMPLICITLY, not computed and stored at construction but
RETRIEVABLE by consumers on demand. This extends the philosophy the
type already applies to the Hermitian partner blocks (section 6: one
block per sector pair stored, partners reconstructed at access time to
dense() by the kind sign): the complete doubled-basis matrix of an
operator is a documented CONTRACT backed by three retrieval mechanisms
of increasing cost:

1. Algebra from stored numbers - the computational-gain route
   (reinstated 2026-11-19: the per-family ruling of section 2d bans
   ASSUMING the relation for untested families, not USING it where the
   unit arm confirms it):
   (a) the Hermitian partner blocks, kind sign - the existing algebra
   of dense();
   (b) the image-image blocks of bit-0-clear operators:
       <Tu|O|Tv> = tau(O) * conj(<u|O|v>) - they cost a sign,
       not a fill. tau(O), the operator's time parity
       (T O T^{-1} = tau O), is recorded per operator as a tau column
       in the producer's class table (15.2) - in code, a named integer
       parameter per family in the generator module (the accessor is
       per-operator, so its own code path is the metadata; no runtime
       lookup, no sym_matrix component). The value is asserted per
       family and ARBITERED EMPIRICALLY: the unit arm computes one
       image-image block both ways (honest kernels and sign flip)
       before the algebra path is armed for that operator. Untested
       families default to the honest recompute of mechanism 2.
   Arithmetic caution, learned the hard way (2026-11-19):
   tau is a property of the operator as written, not of its physics
   family - T is antilinear, so the i of p = -i(nabla) is what makes
   the momentum time-odd; nabla itself is time-EVEN
   (tau = +1, user-confirmed 2026-11-19; the gradient family's
   relation is TT = +SS in single-half storage), NOT the -SS a
   first draft of this section claimed, borrowing the momentum's
   parity. The composite OT inherits the parent's tau:
   T(OT)T^{-1} = tau(O)(OT). The tau table is therefore per-operator
   and oracle-checked, never intuition-fed.
2. On-demand kernels: the cross blocks of bit-0-set operators - the
   fold content of section 2d - are irreducible, nothing stored
   determines them. They are computed when requested, by the honest
   kernels on the TimeReverse images (the TR flags of
   Transfer_psi/Transfer_derpsi). The cost equals the eager twin fill
   of 15.2, paid lazily, once per consumer request (cacheable by the
   caller within a consumer call).
3. Structural zeros: there are no other blocks - the XOR pairing
   annihilates them (no stored x image blocks for bit-0-clear
   operators, no stored x stored blocks for bit-0-set ones).

Placement: the retrieval contract lives in the GENERATOR MODULES -
block_matrix_module stays a plain data container (no stored procedure
components, section 5; no T-knowledge, section 2d), and the T-status
gate stays with the caller: only `$TR`-marked call sites may ask for
image-involving blocks, the standing marker policy. One accessor per
family is the suggested shape:

    ! Any block of the complete doubled-basis matrix of \nabla_mu:
    ! stored and algebraically reconstructible blocks served from the
    ! stored sym_matrix, image-involving blocks materialised by the
    ! kernels on demand; the $TR call-site gate supplies the images.
    function nabla_block(mu, srow, scol, basis) result(A)

Packaging note: for consumers contracting against stored-basis
quantities (rho_can, kappa_can over nwt), the natural packaging of
the cross blocks is the section 2d twin - the image columns relabelled
by their stored partners (the stored-k <-> image-k correspondence) -
so that the `$TR` consumer code keeps the shape of the legacy loops.
The five generators of 15.2 are then the EAGER materialisation of
what the accessor yields lazily (identical numbers, the legacy A/B
unchanged), and implementors may keep them behind it. Consumers under
the implicit design:

- zmomentum: stored-stored blocks of `\nabla_z` only - as today;
- `CompCOMCorrection` `$TR`: requests the fold content via the
  accessor, packaged as the twins, once per call - numerically
  identical to the eager design, with no family-wide doubled storage;
- `$NTR` consumers: the full basis already is the object; the
  accessor degenerates to the stored blocks.

What the implicit design buys over the five-generator eager one:
numerically nothing changes, the consumers keep their code shape; the
win is conceptual - three physical operators with one retrieval
contract instead of five named objects - and economical, no eager cost
for content nobody consumes. The doubled-basis oracle of 15.6 remains
the arbiter: the unit arm exercises the contract by requesting every
block of the complete matrix and judging it against the legacy
routine's honest half and fold - and, per mechanism 1, the tau table
by one image-image block per operator computed both ways.

Connection to the legacy bookkeeping (to be audited at migration
time, not assumed): CompCOMCorrection's `$TR` factors - fac/4 with
"rho_can is twice too large if T is conserved" and the doubling
"2 * COM2_*_debug" for "we summed over only half the states" - look
like the mechanism-1 algebra applied by hand inside the consumer: the
assumption that the omitted T-pairs' contributions mirror the summed
ones. The tau table makes those factors verifiable by the unit arm
instead of maintained as magic numbers; whether they are pure
T-pair mirroring or also carry the fermionic-T pairing phase is an
open empirical question for the migration's A/B.
