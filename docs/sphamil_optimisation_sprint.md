# Sprint: apply_sphamil optimisation (branch `optimize_sphamil`)

## Scope and motivation

`apply_sphamil` (in `src_orig/functional.f90`, body emitted by Hephaestos for
the `$SKYRMEACTION` marker) applies the single-particle Hamiltonian to one
spwf. It is called per spwf from `apply_sphamil_block`
(`src_orig/evolution.f90`) and dominates the wall time of most MOCCa
calculations. Its loop structure is generated code, so optimising it means
touching `src_heph/heph_fields.py::GenerateAction` and/or the derivative
kernels of `src_orig/derivatives.f90`.

## Cost structure of the generated action (per apply_sphamil call)

Generated bodies inspected in `build/<CONFIG>/src_mf/functional.f90` for
NLO, N2LO, BXL and BXL-N2LO (regenerate with
`python3 Hephaestos.py <CONFIG> mf 0 build/<CONFIG>/src_mf src_orig/functional.f90`):

- The bodies of N2LO and BXL-N2LO are **byte-identical**: the BXL difference
  lives in the pairing fields, not in the action of the s.p. Hamiltonian.
- BXL (vs NLO) unfolds `F_Nm_Nm` into the full Cartesian tensor
  `F_N_N(i,a,b,it)`: three times the axpy sweeps; the off-diagonal
  `(a,b)`/`(b,a)` sweeps are redundant if `F_N_N` is (it is — confirmed) but
  exploiting the symmetry systematically in Hephaestos is deemed too hard
  for this sprint.

Per-call cost items, in generated order:

1. `F = Fin` deep copy of the whole `PotentialVector` plus
   `combine_potentials` (Coulomb + constraints + cranking), executed per spwf.
2. On-the-fly branch: `Derive_tot` per spinor column (k = 1..4).
3. Kinetic-energy loop over `ddpsi` diagonal terms.
4. `$SKYRMEACTION`: per field term a sequence of
   - whole-array `temp = 0.0` wipes (13 in N2LO),
   - one separate `do i=1,mv` axpy sweep per field component
     (up to 27 per block, e.g. `G_NkNk_NS` in N2LO),
   - batched `Derive_X/Y/Z(temp(:,:), sym, dtemp(:,:,.))` on (mv,4) — the
     s.p. GEMM dimension is only 4 — or, for N2LO second-order terms,
     four per-column `Derive_lap` calls (twice per `F_NmNm_NkNk` term),
   - whole-array `hpsi(:,k) = hpsi(:,k) +/- c*temp(:,k)` accumulation passes
     per block.
5. The per-spwf callers pass trailing-dimension slices (`x(:,:,wave)`,
   `dx(:,:,:,wave)`); these are contiguous (leading dimensions complete),
   so no argument copies occur — an earlier suspicion that gfortran copies
   them was wrong. Batching's payoff is elsewhere: wider GEMM s.p.
   dimensions and amortized `F` loads, not copy removal.

## Agreed candidate list and discussion outcome

1. **Hoist `F = Fin` + `combine_potentials` out of the per-spwf loop.**
   Keep the public interface of `apply_sphamil` unchanged;
   `apply_sphamil_block` combines once into a private work vector and calls
   an internal routine that assumes pre-combined potentials. The public
   `apply_sphamil` keeps copy+combine for its remaining callers. No dirty
   flags, no invalidation logic (rejected: silent-wrong-answers risk).
2. **Batch across spwfs.** `x(mv,4,m)` reshapes contiguously to
   `(mv, 4*m_block)`, so the DGEMM s.p. dimension widens from 4 to
   `4*m_block`. Requires wider batched derivative kernels in
   `src_orig/derivatives.f90` first, then a batched sibling of the generated
   action. The shared `F` arrays amortize across spwfs in the batched loops:
   this is where the memory-bandwidth win lives.
3. **Fuse the generated axpy sweeps.** Honest estimate: saves only the
   `temp` read-modify-write traffic and loop overhead, ~10-20% if `temp`
   stays cache-resident, up to ~2x if it spills. Not transformative alone,
   but compounds with 2 (the `F` loads amortize over the batch). Do it while
   building the batched variant, not as a standalone.
4. **Reduce `hpsi` accumulation passes.** Non-derivative blocks can
   accumulate directly into `hpsi` (no `temp` round-trip); derivative
   kernels with GEMM-style `alpha/beta` accumulation
   (`hpsi += c * D temp`) absorb the final pass into the kernel call.
5. **Batched `Derive_lap`.** New sibling with per-column signs, mirroring the
   `Derive_X` batched interface; mainly worth it as part of the kernel
   widening of 2.
6. **`F_N_N(a,b) = F_N_N(b,a)` symmetry**: confirmed symmetric, but dropped
   from scope (hard to exploit systematically in the generator).

## Staging

- **Step 0 (this sprint, now): instrument and time.** Split `T_sphamil`
   into sub-timers (see below) and get numbers for NLO and N2LO before any
  optimisation lands; establish golden baselines first.
- **Stage A (low risk)**: candidates 1, 4, 3.
- **Stage B (structural, main prize)**: candidate 2 with widened
  batched derivative kernels; removes the non-contiguous slice copies.

## Instrumentation (step 0)

The timer table is a nested context tree (`print_all_timers_aux`), so
timers started inside a running timer appear as its subtimers — the kernel
time attributed to `apply_sphamil` comes out of the box, without touching
the generator. `golden_diff.sh` strips the timers table, so goldens are
unaffected.

New timers:

- `T_sph_copy` (Sphamil: pot copy+combine): wraps `F = Fin` +
  `combine_potentials` inside `apply_sphamil` (`src_orig/functional.f90`).
- `T_sph_dtot` (Sphamil: on-the-fly derive): wraps the `if(OnTheFly)`
  block, i.e. the `Derive_tot` calls.
- `T_sph_action` (Sphamil: potential action): wraps the kinetic part +
  `$SKYRMEACTION` block, including the `Derive_X/Y/Z/lap` calls it makes.
- `T_der_kernel` (Derive kernel time): wraps the batched spwf kernels
  `Derive_X_spwf`, `Derive_Y_spwf`, `Derive_Z_spwf` and the
  `Derive_lap_{3D,1d}{,_complex}` bodies (`src_orig/derivatives.f90`).
  Nested contexts separate sphamil's kernel share from other callers.

Attribution: within the printed tree, `Sphamil: potential action` minus
its `Derive kernel time` subtimer = candidate-3 sweep time;
`Sphamil: pot copy+combine` and on-the-fly derive are direct read-outs.

## Implementation plan (proposed, agreed scope: staging below)

Corrected cost model after step-0 timing (share of `T_sphamil`):
pot copy+combine 21-23%, DGEMM kernels ~36% of the action,
axpy sweeps + accumulation ~64% of the action, on-the-fly derive 2-4%.
Note: trailing-dimension spwf slices are contiguous, so the per-spwf calls
make no hidden argument copies; batching wins via wider GEMMs and
amortized `F` loads.

### Stage A — incremental, each step one commit

**A1. Hoist the potential copy+combine out of the per-spwf loop.**
- `src_orig/functional.f90`: split `apply_sphamil` into the public
  function (unchanged signature and semantics: copy `Fin`, combine, then
  dispatch) and an internal/public `apply_sphamil_prec(psi,...,F)` that
  assumes a combined `F` and skips the `T_sph_copy` timer (keep timing
  the public one).
- `src_orig/evolution.f90`: in `apply_sphamil_block` and
  `apply_sphamil_block_no_derivative_storage`, do `Fwork = Fin` +
  `combine_potentials` once before the wave loop, then call the
  pre-combined variant. Same for the FAM-side per-wave loops
  (`fam_testing.f90`) after checking they do not modify `F` per wave.
- Expected: removes ~95% of the 21-23% copy share (per block instead of
  per spwf). Risk: low. Separate `PotentialVector` work copies must be
  per-block; no caching across iterations (explicitly rejected).

**A2. Generator: fused sweeps, no temp round-trip (candidates 3+4).**
- `src_heph/fortran_templates/GenerateAction_templates.py` +
  `heph_fields.py::GenerateAction`: per term block emit ONE
  `do i=1,mv` loop that accumulates all component terms into
  `temp(i,1:4)`, initializing on the first statement (`=` not `+`),
  eliminating the whole-array `temp = 0.0` wipes.
- For blocks whose temp is *not* subsequently derived (pure axpy terms
  such as `G_I_NS` sym=+1): emit the accumulation directly into `hpsi`
  with the block coefficient folded in — no `temp` array and no separate
  `hpsi = hpsi + c*temp` pass.
- Expected: bounded by the sweep share (~46% of action measured as
  "internal"); realistic reduction 20-40% of the action. Risk: medium —
  this is generated code deduplication, must be validated across several
  configs since every `$SKYRMEACTION` differs. One mechanism, many
  expansions: verify with dry-run + goldens for NLO/N2LO/BXL/BXL-N2LO and
  the -T variant.

**A3. Accumulating derivative kernels (candidate 4, kernel side).**
- `src_orig/derivatives.f90`: give the leaf DGEMM kernels an accumulate
  form (GEMM `beta=1` on the output) and add `Derive_X/Y/Z_spwf_acc`
  siblings with signature `(f, px, fx, coef)` implementing
  `fx = fx + coef * dX f`. Thread through the derive interfaces as
  distinct module procedures (no optional args — keep the existing
  signatures call-compatible).
- `GenerateAction`: emit the accumulating call for the final action of
  each derivative block, folding the `0.5d0`-type block coefficient into
  `coef` (removes one full `mv x 4` pass per derivative block: 6 in NLO,
  more in N2LO).
- Expected: part of the 64%-share; ~6-9 whole-array passes per spwf
  removed.

Stage A total expectation: T_sphamil ~0.25 (A1) + part of the action;
to be re-measured per commit. **Each step: build NLO + N2LO + BXL +
NLO-T, run goldens, re-time the Ca40 20^3 run against the step-0
baseline; record before/after in the commit body.**

### Stage B — structural (only after Stage A is landed and goldens green)

**B1. Wide batched derivative kernels.**
- New `Derive_X/Y/Z_block(f(mv,4,mb), px(4), fx(mv,4,mb))` (whole block
  shares `sx/sy/sz` by construction of the block callers).
- Batched DGEMM by contiguous remaps of `(nx,ny,nz,4,mb)`:
  X as `(nx, ny*nz*4*mb)` in one product; Z as `(nx*ny, nz*4*mb)` in
  one product; Y per z-plane with widened `rest` dimension
  (plane loop unavoidable because `y` is not the slowest index —
  acceptable, Y is already the cheapest to widen partially).
- Keep the unbatched kernels; dispatch by block callers only.

**B2. Batched generated action.**
- `GenerateAction` gains a wave dimension in the emitted code
  (`psi(i,k,w)`, `do i` fused across `w=1,mb`), so `F%...` loads
  amortize over the whole block. This is the largest generator change;
  the per-spwf variant stays as the fallback for non-batched callers
  (FAM paths, `fam_testing`).
- `apply_sphamil_block` calls one `apply_sphamil_batch(...)` for the
  whole block; wavefunction arrays `(mv,4,mb)` pass directly.
- Care: the sigma-matrix structure couples spinor components *within*
  one spwf; the emitted index arithmetic must respect the `(k,w)`
  layout (`k` is the second, `w` the third array index).
- Expected: attacks both the kernel share (wider GEMMs — the DGEMM
  s.p. dimension goes from fixed mesh widths to mesh-width x 4 x mb) and
  the sweep share (F-reuse); combined with A2/A3 this is the main prize.
- Risk: high — largest blast radius, all configs affected. Land behind
  the per-spwf fallback switch initially (e.g. block size 1 calls the
  old path) so bisecting a golden failure stays possible.

**B3. Batched `Derive_lap`** (per-column signs, mirroring the `Derive_X`
  batched interface), emitted by the generator for the N2LO-style blocks.
  Only as part of B2; standalone gain is small.

### Non-goals (reaffirmed)
- No cross-iteration potential caching / dirty flags (correctness risk
  rejected in discussion).
- No `F_N_N(a,b)` symmetry exploitation in the generator.
- No MPI work; changes stay MPI-plausible and are tested serially.


## Results (step 0, baseline timing)

Setup: gfortran 13.3 `-O3`, serial, EV8-symmetry configs. Small run: O16,
12^3 mesh (dx=1.0), as in `minimal.sh`. Representative run: Ca40, 20^3
mesh (dx=0.8), nwn=nwp=25, maxiter=100 — run from the throwaway directory
`timing_work/` at the repo root (executables and parameterizations copied
in; `SLy4.param` for NLO, `SN2LO1.param` for N2LO; untracked).

### Representing run: NLO (total 11.07 s)

```
19     Sp. Hamiltonian (5848 calls)          2.7276 s   24.6% of total
20       SpHamil: pot copy+combine            0.6310 s   23% of T_sphamil
21       SpHamil: on-the-fly derive           0.1130 s    4% of T_sphamil
22       SpHamil: potential action            1.9828 s   73% of T_sphamil
           (internal: sweeps/accumulation)    1.2722 s   64% of the action
23         Derive kernel time (35088)        0.7105 s   36% of the action
```

### Representing run: N2LO (total 24.07 s)

```
19     Sp. Hamiltonian (5752 calls)          6.4038 s   26.6% of total
20       SpHamil: pot copy+combine            1.3237 s   21% of T_sphamil
21       SpHamil: on-the-fly derive           0.1518 s    2% of T_sphamil
22       SpHamil: potential action            4.9274 s   77% of T_sphamil
           (internal: sweeps/accumulation)    3.1383 s   64% of the action
23         Derive kernel time (97784)        1.7891 s   36% of the action
```

### Reading

The 12^3 O16 smoke runs show the same splits (copy+combine 26% of
T_sphamil for NLO, kernel share of the action 32-35%), so the ratios are
robust across mesh sizes and functional orders.

### A1 implemented (and measured)

`apply_sphamil` split into the public wrapper (unchanged interface: copy
`Fin` + combine + dispatch) and `apply_sphamil_prec` (takes combined
potentials, starts all `T_sphamil` subtimers). All per-spwf callers
(`evolution.f90`: evolve_graddesc, evolve subspace paths,
`apply_sphamil_block`, `apply_sphamil_block_no_derivative_storage`,
spwf dispersions, max-spe estimate; `fam_testing.f90`:
`test_HFmatrices_me`) combine once per loop/block into a local
`Fcombined`/`Fwork`.

Ca40 20^3 (same call counts, energies bit-identical):

| | NLO before | NLO after | N2LO before | N2LO after |
|---|---|---|---|---|
| T_sphamil (s) | 2.7276 | 2.0540 | 6.4038 | 4.7242 |
| sphere: total run (s) | 11.074 | 10.409 | 24.066 | 22.118 |

i.e. T_sphamil -25% (NLO) / -26% (N2LO), whole run -6% / -8%; the
residual combination cost is once per block instead of once per spwf.
Validation: `minimal.sh NLO` passes (O16 energy/deformation vs reference);
N2LO Ca40/O16 energies unchanged; `fam_unit_tests.sh NLO` passes (covers
the modified FAM path); NLO-T builds (its `minimal.sh` input combination
is not supported — dies at input parsing, pre-existing). Pre-existing,
unrelated failures observed: `make CONFIG=BXL mf` references `F%F_Nm_Nm`
absent from the BXL PotentialVector (reproduces with pristine sources);
`store_derivatives.sh` harness needs `-p` and the `LC_ALL=C awk` fix.

### A2 implemented (and measured)

Change in `GenerateAction` (generator): the whole-array `temp = 0.0`
wipe before each contribution group is replaced by a slot
initialization (`=` instead of `+=` with the line sign folded as a
leading `-` if needed) on the first four-line group of each block.
The generated code otherwise keeps the old structure: one mesh sweep
per contribution group, temp accumulated across groups of a block,
derivatives, one final add to hpsi.

Measured (Ca40 20^3, two runs each, energies bit-identical to A1):

| | A1 | A2 | gain |
|---|---|---|---|
| NLO potential action (s) | 1.940 | 1.879 / 1.891 | ~3% |
| N2LO potential action (s) | 4.577 | 4.413 / 4.476 | ~3% |

Negative results, recorded so we do not retry them blind:

- **Sweep fusion across contribution groups of a block is a loss** (NLO
  action 1.94 -> 2.27 s, +17%): one large `do i=1,mv` body per block
  wrecks the compiler's per-group vectorization/register blocking.
  The per-group sweeps stay.
- **Folding derivative-free blocks directly into hpsi is neutral at
  best** (1.85 s without vs 1.87-1.89 s with, when restricted to
  single-sweep blocks; unrestricted it multiplies hpsi read-modify-
  write traffic and costs 13%). Removed entirely: the temp round trip
  is not the bottleneck at group granularity.

The remaining ~3% is the honest ceiling of the wipe/reshape family; the
next real gain is Stage B (F reuse and wider GEMM s.p. width across
spwfs in one sweep).

Validation for A2: `minimal.sh NLO` (O16 -128.5131254465, unchanged),
Ca40 NLO/N2LO energies bit-identical, `fam_unit_tests.sh NLO` passes,
NLO/NLO-T/N2LO mf and NLO fam build clean.



- T_sphamil is 25-27% of a mean-field-dominated run — the single largest
  sub-block of `Spwf evolution` (41%/34%), itself the largest timer.
- **candidate 1 (hoist pot copy+combine) targets 21-23% of T_sphamil**
  and is nearly free to eliminate (per block instead of per spwf).
- **candidates 3+4 (sweep fusion, accumulation passes) target 64% of the
  action**, i.e. `sweeps ≈ 2x kernels`. Consistent across configs.
- **candidate 2 (batching)** attacks both halves: kernels are only 36% of
  the action because the GEMM s.p. width is 4; the sweeps reload the same
  `F` rows per spwf and amortize completely under batching.
- On-the-fly `Derive_tot` is marginal in this workflow (2-4%): most calls
  come with precalculated derivatives.

Per-call overhead sanity check: 26 µs/call (copy+combine) at mv=8000 is of
the right order for a deep copy of ~50 potential fields; the timers'
own overhead (system_clock, ~0.1 µs per start/stop, ~10^5 calls) is
negligible.

Caveats: single machine, one run each, iteration cap at 100 (Ca40 not
fully converged — irrelevant for fractions), and these runs use the
`Rz,T,P,STy` reduced-symmetry builds only; broken-T runs use different
(complex) paths.

## Log

- Sprint document written; instrumentation planned (see above).- Instrumentation added: `T_sph_copy`, `T_sph_dtot`, `T_sph_action` in
  `src_orig/functional.f90` (`apply_sphamil`), `T_der_kernel` around the
  batched spwf kernels and the outer lap wrappers in
  `src_orig/derivatives.f90`; timer IDs registered in `src_orig/timing.f90`.
- `use timing` added to the `derivatives` module (needed for the kernel
  timers; previously the module did not use `timing`).
- Baseline timings recorded above (NLO + N2LO, gfortran, serial).
