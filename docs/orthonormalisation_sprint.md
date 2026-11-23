# Sprint: energy-ordered Cholesky orthonormalisation (successor of `optimize_sphamil`)

## Scope and motivation

`GramSchmidt` (`src_orig/wavefunctions.f90`) is the serial default
orthonormalisation. It executes a classical modified Gram-Schmidt over
each symmetry block with one Fortran whole-array inner product and one
axpy per spwf pair:

- O(N^2) sequential passes over the full `HFPsi` array - bandwidth
  bound, no BLAS, no threading potential.
- Its one structural virtue: it orthonormalises **in order of ascending
  single-particle energy** (via `OrderSpwfsSym`, when `diagsphamil`),
  and it respects the **locked/active** partition of the spwfs. Both
  properties feed the convergence machinery (locking by index, Fermi
  filling by index, max-spe estimate), so a replacement that loses
  either one changes the physics trajectory of a calculation, not just
  its speed.

The dense-algebra alternative already exists:
`Cholesky_orthonormalisation` reshapes each block to a contiguous
`(4mv, N)` matrix through a pointer remap and performs one DGEMM
(overlap matrix), one `dpotrf` (Cholesky factor L) and one DTRSM
(`X = Psi L^-T`). It is the default under MPI, unused serially.

## Measured state of affairs (from the sphamil sprint, Ca40 20^3 NLO)

At 800 spwfs (nwn=nwp=400, 10 iterations, gfortran, serial):

| | Gram-Schmidt | Cholesky |
|---|---|---|
| Orthonormalization | 24.54 s | 0.84 s (29x) |
| whole run | 153.7 s | 123.9 s (-19%) |
| energy trajectory (10 it) | identical to 10 digits | |

All other blocks (sphamil, densities, HF derivatives) match to ~1%.
So the routine-level win is large and real, and in this 10-iteration,
no-locking, well-separated-guess test the trajectories were
indistinguishable. Two obstructions stand between this and a serial
default switch:

1. **Ordering.** Cholesky mixes the block's spwfs into a triangular
   combination. The span is right; the individual columns and hence the
   slot semantics (which column is "the lowest state") are lost. In
   longer runs with locking and mid-run Fermi filling this is suspected
   to matter for convergence.
2. **Locked spwfs.** The serial GS routine orthonormalises active
   states against locked ones and leaves locked states untouched.
   The block-level DTRSM in the Cholesky routine has no such concept;
   it reprocesses locked columns like any other.

There is also an unexplained side observation to resolve before any
default change: at 50 spwfs, CHOLESKY was reproducibly ~1 s *slower*
overall (HF-basis Derivatives +20%, Density: ph +40%) with identical
call counts, while orthonormalisation itself was 3x faster. Suspicion:
denormals in the cancellation-laden spwf tails after a single triangular
solve. Verify by scanning `HFPsi` for magnitudes below ~1e-300 after a
Cholesky orthonormalisation and, if confirmed, flushing tails to zero
(or checking FTZ behaviour) - the fixed-loop blocks cannot otherwise
differ in cost.

## (a) GS-equivalent Cholesky

### The identity to exploit

Classical (not modified) Gram-Schmidt applied to the energy-ordered
columns of a block computes exactly the triangular factor that
`dpotrf` computes on the overlap matrix of those same ordered
columns: if `S = Psi^T Psi = L L^T`, then `X = Psi L^-T` is the
GS basis (up to the signs of the normalisation, made positive by
`dpotrf`'s convention). The ordering is therefore recovered by a
**column permutation only** - the dense algebra itself needs no
change.

### Implementation plan

1. **Order.** Obtain `indices = OrderSpwfsSym(b)` per block (already
   used by `GramSchmidt` when `diagsphamil`; needs a fallback
   definition when `diagsphamil` is false, where GS falls back to
   storage order - Cholesky should inherit exactly that fallback).
   Check how the overlap DGEMM cost compares between orderings: none,
   `S` is permutation-similar (P^T S P), cost identical.
2. **Permute.** Bring the block's `HFPsi` columns into the sorted
   order before the reshape. Options:
   - physically permute the `HFPsi(:,:,:)` columns in place
     (~one pass over 4mv x N, 25 MB at 800 spwfs - milliseconds,
     negligible against the 24 s GS cost), or
   - fold the permutation into the reshape pointer (gathers of
     non-contiguous columns are copy-unfriendly; the explicit permute
     is simpler and honest).
   If the callers require the original slot ordering afterwards,
   scatter back after the DTRSM - decide whether re-scattering or
   re-numbering (rho array, symmetry tables, `spwf_map`) is cheaper;
   re-numbering avoids the second copy but touches the mapping tables.
3. **Solve.** The existing `dpotrf` + DTRSM then yield the
   energy-ordered GS basis to rounding. Note the difference classical
   vs modified GS only affects numerical stability, not the result:
   both sprint routines already coexist on this (the Cholesky header
   is explicitly "happily ignorant" of stability issues).
4. **Numerical acceptance test.** Frozen-orbitals off, compare
   `HFPsi` after `GramSchmidt` vs ordered-Cholesky on a converged
   state: max deviation should sit at rounding level relative to the
   GS result, and `X^T X = I` to ~1e-14. Acceptance: Ca40/O16 minimal
   goldens bit-identical energies, `fam_unit_tests.sh` passes,
   long-run convergence trajectories (E vs iteration, and the locking
   log) indistinguishable from GS.
5. **Measure.** Repeat the 800-spwf Ca40 benchmark and a *long*
   (100+ iteration, locking enabled) run; record the ortho timer and
   the iteration count to convergence. The 10-iteration window of the
   sphamil sprint could not exercise the convergence question.

### Validation matrix

NLO and N2LO mf builds, `minimal.sh` goldens, locking on and off,
`diagsphamil` on and off, and one MPI build (the ordered variant also
fixes the MPI path, which currently lives with unordered Cholesky).

## (b) Locked-state handling

### What GS does today

For each spwf `nw` in the block, in energy order:
- if `nw` is locked: untouched (used as a fixed direction);
- if `nw` is active: project out **all locked** spwfs of the block
  (inner products + axpys against each locked column), then
  normalise, then subtract itself from all later active spwfs.

Locked states thus pin the basis; the active set is orthonormalised
against and within that pinned subspace.

### Dense formulation

Partition the block's columns as `Psi = (Q_l, Q_a)` (locked, active).
The GS result is: `Q_l` unchanged, `Q_a` transformed by
`Y = Q_a - Q_l (Q_l^T Q_a) + ...` and mutual orthonormalisation of
the active block. In matrix form with the block overlap matrix,

```
S = Psi^T Psi = ( Q_l^T Q_l   Q_l^T Q_a )
                ( Q_a^T Q_l   Q_a^T Q_a )
```

the GS-with-locks outcome is obtained by a block-bordered
construction (the `Phi` = projected-active overlaps have their own
Cholesky factor) - i.e. it remains all BLAS + one dpotrf per active
sub-block plus a DGEMM for the projections of active against locked.

Concretely:
1. `Qi_a = Q_l^T Q_a` (DGEMM), `Ca = Q_a^T Q_a - Qi_a^T Qi_a`
   (DSYRK/DGEMM),
2. `dpotrf` on `Ca`,
3. `X_a = (Q_a - Q_l Qi_a) Ca^-T` (DGEMM + DTRSM) - locked columns
   never modified.

This is the same asymptotic cost as the unpartitioned Cholesky
(leading term still the overlap DGEMM) and it drops the per-pair
axpys against locked columns entirely.

### Plan and acceptance

1. Implement as `Cholesky_orthonormalisation` option or a sibling;
   keep the public `Orthonormalize` interface point untouched.
2. Equivalence test against `GramSchmidt` with locks in play: locked
   columns byte-identical, active columns equal to GS to rounding.
   A locked-heavy case (converged heavy nucleus, locking_treshold
   active) is the right test vehicle.
3. Acceptance: minimal + fam goldens; a locking production run
   comparing iteration-by-iteration energies and the lock journal.
4. Only after ((a) and (b) land, and the denormal question is
   resolved): propose `ortho_strategy` default 'CHOLESKY' for serial.
   The LOEWDIN stub (`stp('Implementation error')`) can stay stubbed.

## Secondary: copy and allocation cleanup in the subspace rotation

Code-level analysis of the serial `calc_sphamil` +
`apply_subspace_rotation` (evolution.f90, `T_calc_sph` / `T_subrot_*`
timers). Not measured yet - the Ca40 benchmarks of the sphamil sprint
ran with subspace rotation INACTIVE. The algebra (matrix-element
DGEMM, DSYEV, rotation DGEMM) is already well-shaped; the wins are in
the Fortran shuffling around it. Expected size: single-digit percent
of the subspace step, which itself is dominated by the ~800
`apply_sphamil_prec` evaluations it launches (0.28 s per call-set at
800 spwfs) - tidy-up, not a 29x.

- **`calc_sphamil`: remove the `temp_wfs` copies.** The per-spwf
  `temp_wfs(:,:,i) = hfpsi(:,:,wave)` exists only to feed
  `apply_sphamil_prec` (which does not modify its input; the slice
  `hfpsi(:,:,wave)` is contiguous and can be passed directly) and to
  build the contiguous matrix-element DGEMM operand. `HFPsi` is stored
  `(mv,4,N)`, which is memory-identical to `(4mv, N)` - so with no
  locked spwfs in the block, the DGEMM operand is a free pointer
  remap of `hfpsi(:,:,si+1:si+N)`. Gathers remain only for
  locked-containing blocks, over the active columns only.
- **`apply_subspace_rotation`: five passes over the wavefunction set
  to two or three.** Today: per-(k,wave) gather into `temp_wfs` (1),
  `temp_wfs_in = temp_wfs` full copy (2, with `temp_wfs_in` allocated
  by implicit assignment per block), DGEMM read/write (3, 4),
  per-(k,wave) scatter back (5). With no locks, gather and scatter are
  identity remaps: rotate straight out of `HFPsi` into a scratch
  buffer, then one copy back. With locks, keep the gather of active
  columns but drop the `temp_wfs_in` copy by rotating into scratch and
  scatting only the active columns. `Momentum_Updates` gets the same
  fix.
- **Allocation churn**: `hpsi`, `temp_wfs`, `active_sph`,
  `active_indices` per block per call in `calc_sphamil`; the DSYEV
  workspace query + reallocation per block in
  `apply_subspace_rotation`. Hoist all of it to maximum-block-size
  workspaces.
- **DSYEV vs DSYEVD/DSYEVR**: at per-block N around 100 the
  divide-and-conquer routine may be faster; benchmark, do not assume.
- Cosmetic: the `active_sph`-into-`sph` scatter and the `findloc` per
  spwf in the eigenvalue bookkeeping are trivial next to the above.

Acceptance for all of the above: bit-comparison of the rotated
wavefunctions and eigenvalues against the current implementation on a
rotation-heavy run, plus the minimal/fam goldens.

## Non-goals

- No change to the GS routine itself (it remains the reference
  implementation and the fallback).
- No MPI-side redesign beyond reusing the ordered/locked-aware block
  routines where ScaLAPACK allows; testing validity with MPI remains
  the developer's explicit choice.
