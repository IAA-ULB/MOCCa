# AGENTS.sprint.md — Growing the spwf space in iterate_spectrum

Sprint: extend `scfscheme = 2` (frozen-potential spectrum iteration, see
`iterate_spectrum` in `src_orig/mocca.f90`) so that single-particle wavefunctions
are added *during* the iteration until the converged single-particle spectrum
covers a requested excitation energy above the Fermi level. Hartree-Fock only
(`pairingtype = 0`); no pairing-deformed objects this sprint.

Agreed decisions (2026 sprint kickoff):

- **Trigger**: convergence-gated. A new batch of spwfs is added once the
  previous batch is locked (dispersion below threshold). Debug knob first:
  fixed batch size 10, no cadence limit. Fixed-cadence fallback is NOT wanted.
- **Dev/test config**: `NLO` (time-reversal conserved), built with
  `INCLUDE=gnu-debug` for the inner loop; symmetry-marked code must also
  compile+run for NLO-T later in the sprint.
- **Reuse**: `TransformInput` in `src_orig/transform.f90` already performs the
  spwf-array extension (random-init new spwfs, per-block bookkeeping,
  loadbalance, dispersions, `rho_can`); the pairing branches are inert at
  `pairingtype = 0`. The new code in `mocca.f90` covers what it does not:
  `locked`, `sx/sy/sz` via `set_spwf_symmetries`, derivative arrays
  (`deallocate_derivatives` + `allocate_memory_derivatives`), and the per-spwf
  observable arrays (`HF_*`, `CAN_*`, `qp_*`).
- **Batch recipe**: the N extra spwfs are distributed equally over the
  symmetry blocks allowed by the build's `$TR/$NTR/$PBROKEN` guards, which
  `TransformInput` already polices.

Steps:

1. Extract the commented block inside `iterate_spectrum` into a routine
   (`add_spwfs`), driven by `TransformInput`. Milestone: one successful
   "add 10 spwfs at iteration k" run that keeps converging afterwards (NLO,
   gnu-debug).
2. Convergence-gated trigger in the `iterate_spectrum` loop.
3. Termination: stop adding (and let the iteration finish) once the highest
   locked spwf energy exceeds e_F + excitation cutoff; final convergence
   requires all below-cutoff spwfs locked.
4. User input: excitation-energy cutoff keyword(s) under the `scfscheme = 2`
   input; batch size as an optional knob.
5. Printing: `print_summary_spectrum` reports the added batch and the current
   spectrum ceiling.
6. Test: compare against a plain `scfscheme = 2` run started with all spwfs
   pre-loaded — spectrum below the cutoff must agree; then a golden reference
   (testing/integration, `.reference` convention).

Findings so far (session 1):

- Implementation lives in `src_orig/mocca.f90`: `iterate_spectrum` now houses
  the convergence-gated trigger (`nbatches`, `spwfs_per_block`, the temporary
  `max_growth_batches` cap), `spectrum_growth_blocks(nex)` builds the
  marker-gated batch recipe, `add_spwfs` performs the growth via
  `TransformInput`. Verified by Hephaestos dry-runs that the recipe is valid
  in all four symmetry scenarios (NLO / NLO-T / NLO-P / NLO-TP).
- `TransformInput` also resizes `hftransfo` and `sphamil` (its steps 6/7); do
  not redo that in `add_spwfs`.
- Critical: `TransformInput` zeroes `rho_can`; the (dt,mu) estimation
  (`IterativeEstimation`) and the weighted dispersions need HF occupations.
  Without a `SolvePairing` refill after growth, `relE` blows up to 1e8, dt
  collapses to ~1e-4 and the iteration stalls. `add_spwfs` therefore calls
  `SolvePairing` (HF only, pairingtype 0) after the growth.
- Observed on O16/SLy4/NLO 12^3 (gnu-debug): growth fires correctly at
  iteration ~110-170 (threshold-dependent), bookkeeping is clean (no crashes,
  bounds-check on), dispersions of the new random batch collapse by ~3 orders
  of magnitude in ~10 iterations. OPEN: the tail does not fully reconverge —
  with locking_treshold 1e-2 the max dispersion asymptotes around 2-3e-2 and
  creeps UP slowly; with 0.5 it bounces around 2.8e4. One state per batch
  (presumably the worst-separated eigenvalue of the frozen h, or a
  locked-spwf dispersion bookkeeping issue: locked spwfs keep stale
  dispersion values and `print_summary_spectrum` reports
  maxval over those) never locks, so the second batch never fires. Needs
  discussion before step 3: locking semantics vs stale dispersions, and
  whether the frozen-field diagonaliser needs a better init for new spwfs
  than random (e.g. iterate a few cheap refinement steps before merging
  them, or orthogonalise against the locked subspace only).

## Instability diagnosis (session 2)

- Not a locking-threshold artifact: identical exponential tail growth at
  thresholds 0.5, 1e-2 and 1e-3 (e-fold roughly every 250 iterations from a
  clean minimum ~0.37; at 0.5 the iteration lands in an exact 2-cycle with
  frozen gradn). Not model-space exhaustion either (50 of ~6900 lattice
  states, ceiling ~5 MeV).
- Supported mechanism: when most spwfs lock, the residual active set is
  1-2 states per block; heavy-ball (dt,mu) estimated for the full 50-state
  window makes the residual dynamics unstable, and locked chains are
  one-way, so blown-up residuals never re-lock.
- Fix implemented: `n_sacrificial` (evolution namelist) — the highest
  n_sacrificial states per block are never locked; tested at 1 and 3.
- Frozen stale dispersions of locked spwfs make the "Maximum dispersion"
  print misleading (see QOL_BACKLOG #4).

## dt/mu estimation for the active set (proposal, sprint step)

Current behaviour (`IterativeEstimation`, every iteration when EstimateParams):
- maxE: power iteration on a persistent random spwf (cheap after the first
  call since the field is frozen, but still repeated every iteration);
- relE (HF): the smallest gap between ANY occupied and ANY unoccupied
  spenergy, clamped at 0.2 MeV — a property of the CONVERGED set;
- dt = 4/(maxE+relE+2 sqrt(maxE relE)) hbar safety, mu from kappa = relE/maxE.

Two problems for scheme 2: (1) the window [global min spenergy, mesh maxE] is
the conditioning of the FULL problem, while only the UNLOCKED states evolve;
(2) relE is set by states that are already converged, not by the frontier gap
that actually limits progress of the active set.

Proposed estimator (per block, no extra h applications):
- maxE: compute once per frozen field (cache; refresh inside add_spwfs).
- Frontier relE: from the per-block active-set eigenvalues that
  apply_subspace_rotation's DSYEV already produces at zero cost:
  relE_B = smallest adjacent gap of {highest locked spenergy of block B}
  joined with the sorted active eigenvalues of block B (clamped at ~0.1 MeV).
  Either take relE = min over blocks per isospin, or go per-block with
  dt_B/mu_B from the same heavy-ball formulas in the block's own window.
- Momentum hygiene: reset mu = 0 (dt to the k=1 step) right after each
  growth and whenever a block's active count drops to ~1; optionally flush
  momentum periodically (every ~20 iterations) to kill parasitic modes.

Rationale: in the "explore new directions" mode (diagsphamil=.false.), the
gradient is orthogonalized against the spanned space, so its spectral
ceiling is still maxE, but the rate-limiting low end is the frontier
separation — not the occupied/unoccupied gap. This estimator measures
exactly that, for the subset that actually evolves, using eigenvalues the
loop already computes.

## Sacrificial states: result (session 2, final)

- `n_sacrificial` (evolution namelist): the highest n_sacrificial spwfs per
  symmetry block are never locked; the module exposes the logical array
  `sacrificial(:)` set during the locking pass, and `iterate_spectrum`
  tests convergence as ALL(locked .OR. sacrificial).
- Small n (1, 3) does NOT cure the tail instability: the runaway settles in
  the highest sacrificial state itself (observed d2h up to 8e3 at
  `n_sacrificial=3`).
- `n_sacrificial=10` (with `locking_treshold=1e-2`) WORKS end-to-end on the
  O16/SLy4/NLO 12^3 debug test: convergence-gated growth fires twice
  (30 -> 50 -> 70 spwfs), run converges in 81 iterations, final maxdisp
  stays bounded (sacrificial states hover at d2h ~0.1-0.2). Final energy
  -128.5375 MeV vs -128.5131 MeV initial: the 24 keV drift is the density
  residual of locking at 1e-2, consistent with the threshold.
- At `locking_treshold=1e-3` the gates still open but the batch-2 frontier
  cannot be resolved by the current global (dt,mu): tail instability
  returns. If 1e-3-quality spectra are needed, implement the frontier-window
  (dt,mu) estimator below before raising default thresholds.
- Workflow note: background runs started via `setsid ... &` from bash are NOT
  tracked by the process-notification system; use `tools.process.start` if
  terminal notifications are wanted.

## Batch-relative sacrificial semantics (session 2, final)

- `n_sacrificial` now dual-semantics (evolution namelist):
  - positive: the highest N states of every symmetry block never lock;
  - negative (-k): all spwfs added by the last k growth batches never lock.
- A per-spwf batch ledger `spwf_batch(:)` (wavefunctions module, 0 = present
  from start) is extended by `add_spwfs` after each TransformInput: within a
  block the new spwfs sit behind the old ones. The ledger is not written to
  wf files; a continuation run starts a fresh numbering.
- Test (O16/SLy4/NLO 12^3, gnu-debug, treshold 1e-2):
  `n_sacrificial=-2` fires both growths (30 -> 50 -> 70) and converges in
  254 iterations (slower than the absolute n=10 run's 81, since the
  pre-growth set has no sacrificial reserve and batch-1 alone is thinner);
  deep occupied spectrum reproduced exactly (-36.153, -20.576... as the
  reference run). Final-total-energy drift between runs with different
  sacrificial settings (~25-50 keV) comes from the final-printout density of
  threshold-locked and hovering states; must be quantified by the step-6
  cross-validation against a fully preloaded spectrum run.

## Excitation-cutoff keywords (session 3)

- New input (namelist &scfiteration, scfiteration.f90):
  - `spectrum_ecut` (MeV, real, default -1 = off): grow the spwf space,
    batch by batch, until the LOCKED spectrum reaches
    FermiEnergy + spectrum_ecut for BOTH isospins; scfscheme=2 only (stp
    otherwise). Ceiling per isospin = maxval(spenergies, mask=locked),
    locked states keep their eigenvalue at lock time.
  - `spectrum_batch` (integer, default 5): new spwfs per symmetry block
    per growth batch. Both are broadcast; sane errors; printed in the
    run header. Replaces the hardcoded spwfs_per_block/max_growth_batches.
- Bug fixed on the way: batch-relative n_sacrificial had an off-by-one —
  after the first growth, batch 0 (the initial states) was flagged
  sacrificial (cutoff = current+k instead of max(0, current+k)); the gate
  then saw empty locked sets (maxval -> -HUGE prints as ***).
- Toolchain note: 'make CONFIG=NLO INCLUDE=gnu-serial' only rebuilds files
  newer than their objects — make.inc is NOT in the dependency graph, so a
  debug-to-serial switch leaves most objects at -O0 until 'make clean'
  (global). After full -O3 rebuild: ~20 ms/iteration at nwt=50, ~2x that at
  nwt=150 (freeze WARNING: the whole batch machinery inherits this).
- Verified: gates fire only at genuine convergences (n_sacrificial=-2
  reserve intact), gate prints ceilings, growths chain correctly
  (30->50->70->90 / ->110 / ->150 observed).
- OPEN BLOCKER for the sprint goal: the locked ceiling does NOT advance
  across batches (stuck at the initial 3.553/4.667 n/p MeV in the O16 test,
  both for spectrum_batch=5 and =10). Batches appear to converge into the
  already-spanned band instead of capturing the next band of eigenvalues.
  Suspects to check next, in order:
  1. dump per-batch spenergies/lock status (extend the gate print or use
     the final table) — do batch-k states duplicate batch-(k-1) energies?
     (a near-degenerate pair -36.302/-36.153 was spotted in a final table
     that should have a single 0s neutron state; could also be a table-parse
     artifact mixing isospin tables);
  2. interplay of the sacrificial (never-locked) states with the active
     subspace rotation: do they "steal" the new directions the lockable
     batch states should capture?
  3. whether apply_subspace_rotation's active-submatrix diagonalization
     preserves the complement structure at all when most of the block is
     locked.

## Sacrificial machinery retired (session 3, per W.R.)

- The n_sacrificial/batch-ledger machinery is REMOVED (evolution, wavefunctions,
  mocca): locking is again the plain per-block chain, and the design is the
  agreed overshoot one: batches add spwfs; states above Fermi+spectrum_ecut
  are superfluous by construction (the last batch overshoots the cutoff).
- Three convergence semantics were measured (O16/SLy4/NLO 12^3, ecut=15,
  batch=10, treshold 1e-2):
  1. ALL(locked): honest; exactly one growth fires; the run then blocks
     forever because at least one fresh state hovers just above the locking
     treshold at the frontier (known tail problem).
  2. Spenergy-based exemption (spwf exempt if > Fermi+ecut): TREADMILL —
     fresh states' transient Rayleigh quotients (random start ~ +50..1000)
     look 'above target' before their first descent, convergence is
     declared while the new batch is garbage, and the gate adds batch after
     batch with the locked ceiling frozen at the pre-growth frontier
     (5 growths, ceiling constant 3.553/4.667). A spenergy of an unconverged
     state is NOT a predictor of its eventual eigenvalue; exemption cannot
     be based on it.
  3. (Earlier) k-batch lag + ALL(locked .or. sacrificial): converges and the
     ladder advances honestly, but the ceiling moved only 3.553 -> 3.664
     per batch of 10/block in this dense small box — viable but needs the
     ledger machinery W.R. wants removed.
- CONCLUSION/next step: the clean design (1) is the keeper; its blocker is
  the hovering frontier state. Cure = the frontier-window (dt,mu) estimator
  of the session-2 proposal (per-block heavy-ball parameters from the active
  set's eigenvalue window + momentum reset after each growth). That is the
  next work item.
- Note the QOL entry: 'Maximum dispersion' prints include stale locked values.

## iteration_status design (session 3, final state)

- `locked(:)` is replaced by `iteration_status(:)` (wavefunctions module):
  ITS_ACTIVE (evolve, must converge), ITS_LOCKED (frozen),
  ITS_SACRIFICIAL (evolves forever, exempt from convergence). Parameters
  ITS_ACTIVE/ITS_LOCKED/ITS_SACRIFICIAL live next to it. Sites touched:
  wavefunctions (Gram-Schmidt skips locked), IO (initial allocation),
  printing (LOCK column now Y/N/S), evolution (locking chain, evolve loop,
  calc_sphamil diagonal-fill, apply_subspace_rotation active set),
  mocca (convergence test = ALL(status /= ITS_ACTIVE), gate ceilings over
  ITS_LOCKED mask, add_spwfs rebuilds statuses).
- Sacrificial rule, positional only: sacrifice_from/prev_extra_blocks (2x8
  integers) fix the newest-TWO-batches boundary per block; batches append at
  block ends so no batch ledger is needed. No energy criterion anywhere.
- Printer now shows S for sacrificial spwfs.
- Acceptance trace (O16/SLy4 12^3, ecut=15, batch=10, treshold 1e-2,
  -O3): three growths 30->70->110->150, ceiling ADVANCING honestly
  (3.553/4.667 -> 3.664/5.094 measured over LOCKED states only), status
  counts sainne in the final table (20 S, Y below, N only in the fresh
  batch). Run reached maxiter at 150 spwfs before completing ecut=15 —
  a cost question (each batch reconvergence grows with nwt and the dense
  spectrum of this tiny box), not a correctness one. Next steps remain:
  frontier (dt,mu) estimator, then bigger-box validation + golden test.

## BLOCKER: state corruption from the third growth on (session 3, final)

- Symptom (3+ growths, nwt>=110): collapsed spwfs (spenergy EXACTLY 0,
  d2h ~1e-27, i.e. near-zero norm) among the sacrificial set, and LOCKED
  states with impossible eigenvalues (-119 .. -427 MeV) yet tiny
  dispersions and correct radii/J. Total energy off by 100s of MeV
  (-36.8 / -5.6 / -63.8 instead of -128.5).
- Bisection: PRE-DATES the iteration_status refactor — same corruption in
  the k-batch-lag era runs (accept: -5.58, 12 zero states; ecut25b5v2:
  -63.76). All runs with <= 2 growths were sane (brel -128.56). Nothing
  about the ladder/locking logic is implicated: the corruption lives in
  the growth bookkeeping, first bad quantity after the 3rd add_spwfs.
- ALSO explains two earlier anomalies: the input-path control run whose
  extraspwfs never evolved (frozen 1e9 dispersions), and the perpetual
  in-band 'hovering' (likely index-misaligned states, not physics).
- DO NOT TRUST any 3+-growth run result so far; the ceiling-advance
  measurements remain valid as *mechanism* tests only.
- Debug plan: instrument add_spwfs + TransformInput to dump per-spwf
  norms, spenergies and spwf_map/HFBlocks consistency after each growth;
  run 3-growth O16 case under gnu-debug bounds-check; suspects, in order:
  1) TransformInput's copying loops (the do b=1,8,2 odd-block-only
     spenergies/dispersions/hfpsi redistribution at repeated calls),
  2) the derivatives re-allocation/indexing after the 3rd growth,
  3) HFBlocks_global vs HFBlocks/spwf_map drift at repeated loadbalance.

## RESOLVED: growth-corruption root causes found and fixed (session 4)

- Root cause (the corruption): Evolve_momentum's update loop added
  momentum_updates(:,:,wave) to EVERY spwf of the block, with no status
  check; the gradient loop above it skips LOCKED states, so a locked spwf
  kept receiving its STALE momentum vector once per iteration — a linear
  walk away from norm 1 that GramSchmidt (skip locked) never corrected.
  Dormant with locking_treshold=-1 (upstream default); our sprint enabled
  locking. Fix: skip ITS_LOCKED in the update loop (evolution.f90).
- Root cause (real, secondary): Momentum_Updates keeps its first-call
  nwt_local size across add_spwfs growths while apply_subspace_rotation
  indexes it with post-growth indices — OOB writes. TransformInput
  resizes everything else. Fix: deallocate in add_spwfs; Evolve_momentum
  reallocates at the new size (old heavy-ball history lost, irrelevant).
- In passing: deallocate_derivatives said "$N3 allocate(candddPsi)"
  instead of deallocate — genuine latent bug for N3 functionals +
  pairingtype 2. Fixed.
- Verified (O16/SLy4 12^3, ecut=15, batch=10, treshold 1e-2, -O3):
  every NORMCHK deviation count 0, all norms exactly 1 across the whole
  run incl. both growths 30->70->110; total energy -128.513117 MeV.
  Stripped-binary rerun reproduces the instrumented run to the digit.
- REMAINING (not corruption): after the 2nd growth, 11 non-sacrificial
  states hover at dispersion ~4.3e-2 vs. locking_treshold 1e-2, so the
  run grinds to maxiter without a further growth — a locking-threshold
  calibration question for the added high-lying states. Next up.

## Hovering states: mechanism found, stagnation locks added (session 4b)

- Mechanism: once apply_subspace_rotation makes every active spwf an
  exact eigenvector of the projected h, a state whose residual h|psi> -
  e|psi> points only along LOCKED spwfs is caught in a fixed point: the
  heavy-ball update adds locked-direction components, GramSchmidt
  removes them again, the spwf NEVER moves (verified: dispersion frozen
  to 5 digits over 600 iterations, spwf i=35 E=3.53). all_lower_converged
  then blocks every state above it -> run stalls at maxiter. This is the
  old 'hovering' anomaly; textbook heavy-ball, nothing to do with bookkeeping.
- Fix: stagnation fallback, evolution namelist keyword stagnation_iter
  (default 0 = off). An active spwf whose dispersion has not set a new
  best (<0.9x) for stagnation_iter consecutive iterations becomes
  lockable and ladder-eligible; one audit line is printed per stagnation
  lock. Counters reset on growth (resize), so expect a burst of
  re-prints after each add_spwfs - cosmetic only (most fire on states
  already far below locking_treshold).
- Also: print_summary_spectrum's maximum dispersion is now masked to
  ITS_ACTIVE (was max over ALL incl. stale locked dispersions, printed
  a constant misleading value).
- Verified (O16/SLy4 12^3, ecut=15, batch=10, treshold 1e-2,
  stagnation_iter=150, -O3): 4 growths, ceilings staircase
  3.55->3.67->6.33->8.95 (n), run printing 'Convergence criteria
  satisfied' and finishing cleanly; total energy -128.513117 MeV,
  identical to before. CEILING VALUES ARE NOW TRUSTWORTHY (lock set
  clean) -> the next sprint item is the real acceptance scan over ecut.

## Evolve_spectrum extraction + scfscheme string keyword (session 4c)

- Evolve_spectrum(F, iteration) added to evolution.f90: the heavy-ball
  path used by scfscheme=2, extracted verbatim from Evolve_momentum as a
  first, behaviour-preserving step toward dedicated-spectrum machinery
  (deflated locking is next). run_mean_field now assigns
  Evolve_subspace => Evolve_spectrum for spectrum runs - this overrides
  the Strategy keyword of the evolution namelist for those runs (noted
  in the code and in the input docs below).
- BREAKING INPUT CHANGE: scfscheme is now a string,
  'MEANFIELD' (default) or 'SPECTRUM'. The legacy integers 0/1/2 no
  longer parse; the broken scheme 1 (linear mixing, never operational)
  is gone. No committed test input or reference used scfscheme, so
  only ad-hoc user inputs break. Internal comparisons now read
  trim(scfscheme) == 'SPECTRUM' (mocca, IO) / 'MEANFIELD' (precond
  branch). MPI_BCAST broadcasts the character variable.
- Regressions: O16/SLy4 12^3 spectrum acceptance case reproduces the
  pre-refactor output token-for-token (only the scheme header line
  differs); minimal.sh NLO converges -128.513125 MeV with the default
  MEANFIELD path.
- Known absent facilities on this branch, AGENTS.md's testing section
  notwithstanding: golden_diff.sh, moments/block-matrix unit suites are
  not in this checkout - do not assume they exist.
- NEXT (in the new routine): deflated-residual locking criterion
  ||P(h psi - eps psi)||^2 (P orth. to locked span), then joint
  Rayleigh-Ritz spectrum reporting at growth/convergence boundaries,
  then per-state dt/mu.
