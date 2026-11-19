# FOUND_BUGS.md — issues and inconveniences found, out of sprint scope

Issues (bugs, warts, inconveniences) encountered during the matrix_refactor
sprint that were deliberately NOT fixed or adjudicated within the sprint,
because they are pre-existing, out of scope, or runtime-irrelevant. Each
entry names a reproduction and a workaround. An entry is closed by fixing
the issue in a dedicated commit and flipping the status - do not delete
entries. In-scope sprint problems live in the design doc's amendment sections
and the commit messages, not here.

| # | Status | Where found |
|---|---|---|
| F-1 | OPEN | coulomb.f90 / debug NLO-P builds |
| F-2 | FIXED (2421931f) | wavefunctions.f90 dead kernels |
| F-3 | OPEN (adjudication pending) | linear_momentum arm |
| F-4 | OPEN (workaround) | Makefile / NLO-TP parallel build |
| F-5 | OPEN (harness design) | testing/integration/golden_diff.sh |
| F-6 | OPEN (trap) | testing/unit/fam_unit_tests.sh |
| F-7 | OPEN | fam_spurious_mode_subtraction scenario |
| F-8 | OPEN (cosmetic) | FAM output, parity-broken builds |
| F-9 | KNOWN LIMITATION | cross-build trace comparisons |

### F-1 Coulomb folding out-of-bounds in debug-profile NLO-P builds — OPEN

Pre-existing, found 2026-11-19 (angular-momentum arm attempt), NOT a sprint
regression. The debug-profile (-fbounds-check) NLO-P binary crashes in
coulomb.f90's obtain_folded_potentials (fold_form_factor_reverse, index 26 of
dimension 3 vs extent 24) during set_up_mean_field. An out-of-bounds in the
parity-broken Coulomb folding path that the -O3 CI artifacts never check.
Workaround: do not use debug-profile NLO-P builds; NLO-P gates run on the
frozen -O3 build (green).

### F-2 angmom_x_imag / spin_x_imag spin-pairing transcription bug — FIXED (2421931f)

Term-by-term audit against the validated kernels (design §14.2b) found the
spin down-block of these two caller-less kernels carrying sigma_y pairings
dressed with x-style signs instead of the honest
Im[conj(d_down) c_up] = wf2(3)wf1(2) - wf2(4)wf1(1). Fixed and
banner-annotated as caller-less; no production numerics changed (the arm's
purity print was angmom_x_imag's only reference; spin_x_imag had none).
Kept as a record because the same pattern bit the class-drafting once more,
in reverse: a class tag drafted from a legacy transposed-position fill
(nabla_x.T, §15.5) inherited a similar sign wart until the mirror oracle
flipped it - when wrapping legacy kernels, measure, never trust the tag.

### F-3 Derivative operators' opposite-half chamber residues (design §15.6) — adjudication pending

Found 2026-11-19 by the linear-momentum arm, deliberately UNADJUDICATED: the
honest nabla generators' opposite halves (eta-stripped cross-STy blocks) do
NOT vanish - O(0.2) against O(0.35) element scale, on maxiter=0 AND converged
states, on the crossed-component pairs - while the multipole full-box oracle
proves the basis STy-clean at 1e-16 on the same machinery. The angmom
analogue was adjudicated an artifact (DECISIONS.md record; design §14.2b);
by mechanism the artifact reading is favoured here too, but the derivative
finding stays unjudged until the derivative family gets its own full-box
reconstruction oracle - the angmom arm's check_fullbox is the ready
template. Production unaffected: every consumer contracts real densities and
the legacy agrees to 3e-15 exactly as filled.

### F-4 make -j4 races on missing .mod files (NLO-TP) — workaround: build serial

A parallel (-j4) mf build of NLO-TP terminated on missing .mod files
(2026-11-19, rotations gate). Workaround: build NLO-TP serially.

### F-5 golden_diff.sh teardown destroys its work/ directory — by design, with teeth

The harness owns testing/integration/work and destroys it at teardown; the
original row-0 scratch scenarios of the COM-correction migration were lost
that way (2026-11-19) and had to be recreated. Workaround: dev scratch lives
outside the harness tree (e.g. /tmp). Long-term option: give scratch
scenarios a committed home under testing/ if they must survive.

### F-6 fam_unit_tests.sh default --parameterisation t0t3 stops silently — a green-looking trap

Under an NLO fam build the default t0t3 parameterisation stops at the .func
mismatch with `stp` and exit status 0, without running any FAM test
(2026-11-19). A wrapper reading only the exit code sees green. Workaround:
always pass --parameterisation SLy4 for NLO fam builds. Consider making the
script fail loudly on input mismatch when touched next.

### F-7 fam_spurious_mode_subtraction scenario needs an input flip for NLO-P — OPEN

The committed scenario script's mf input needs nz=16 (unfolded z) before it
can run against NLO-P executables (2026-11-19); as committed it targets the
parity-conserving build, where parity-odd fields NaN by construction. The
PASS criterion itself was validated on NLO-P (spurious-mode subtraction
drives S_after to ~1e-17 for zcom and l=1 m=0). Unfixed.

### F-8 FAM output cosmetics in parity-broken builds — OPEN, cosmetic

In NLO-P strength tables the parity +/- rows are vestigial (all strength
books into parity +), and the box-size check's drho_antisym panel prints
NaN exactly as it does in healthy NLO runs (2026-11-19). Cosmetic; identical
behaviour in the parked pre-refactor anchor.

### F-9 Cross-build converged traces are environment-sensitive — known limitation

Converged cross-build trace comparisons showed ~8e-9 local but ~3e-6 in CI on
the parity-broken l=4 box-anisotropy rows (2026-11-19 CI failure): each build
stops at a slightly different point of a near-flat convergence landscape, so
the comparison is only meaningful on a shared state (maxiter=0) and through
basis-independent invariants. Handled by convention, see SPRINT.md; entry
kept as a warning against re-introducing converged cross-build comparisons.
