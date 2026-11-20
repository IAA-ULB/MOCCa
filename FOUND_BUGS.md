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
| F-1 | FIXED (0d883c19) | coulomb.f90 / debug NLO-P builds |
| F-2 | FIXED (2421931f) | wavefunctions.f90 dead kernels |
| F-3 | OPEN (adjudication pending) | linear_momentum arm |
| F-4 | OPEN (workaround) | Makefile / NLO-TP parallel build |
| F-5 | OPEN (harness design) | testing/integration/golden_diff.sh |
| F-6 | FIXED (b10916c8) | testing/unit/fam_unit_tests.sh |
| F-7 | OPEN | fam_spurious_mode_subtraction scenario |
| F-8 | OPEN (cosmetic) | FAM output, parity-broken builds |
| F-9 | KNOWN LIMITATION | cross-build trace comparisons |
| F-10 | FIXED (cbd45ac6) | namelist input, all modules |
| F-11 | FIXED (3c1af4b4) | CI pasta tests / fission_MOI defaults |

### F-1 Coulomb folding out-of-bounds in debug-profile NLO-P builds — FIXED (0d883c19)

Pre-existing, found 2026-11-19 (angular-momentum arm attempt), NOT a sprint
regression. obtain_folded_potentials sliced F%ExchangePotential - an array on
the original mesh - with the Coulomb-mesh offsets, so any build storing z
unfolded ($FULLZ: coul_offset_z = BC) overshoots the z slice to BC+nz,
crashing under -fbounds-check (index 26 vs extent 24). The direct
CoulombPotential slice is padded and was correct. Fixed by passing the
unpadded ExchangePotential. No production numerics were ever corrupted: the
branch requires coultreatment=1 with folding matrices allocated, a
combination rejected at startup (parameterization.f90); frozen -O3 NLO-P A/B
vs the pre-fix binary is token-identical. Debug-profile NLO-P builds work
again.

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
(2026-11-19, rotations gate). Not config-specific: reproduced on NLO-P
2026-11-19 - the first dispatched consumer whose provider is still
compiling dies with "Cannot open module file '...mod' for reading"
(there: IO.o vs io_wf.mod). Workaround: build serially. The same
missing-.mod-dependency wart bites SERIAL builds too after any module
interface change (observed 2026-11-19, the find_nml_error signature change):
consumers keep compiling against a stale .mod. A changed interface in a
module as widely used as geninfo demands a make clean of the CONFIG.

### F-5 golden_diff.sh teardown destroys its work/ directory — by design, with teeth

The harness owns testing/integration/work and destroys it at teardown; the
original row-0 scratch scenarios of the COM-correction migration were lost
that way (2026-11-19) and had to be recreated. Workaround: dev scratch lives
outside the harness tree (e.g. /tmp). Long-term option: give scratch
scenarios a committed home under testing/ if they must survive.

### F-6 fam_unit_tests.sh default --parameterisation t0t3 stops silently — FIXED (b10916c8)

Found 2026-11-19: under an NLO fam build the default t0t3 parameterisation
stopped at the .func mismatch with `stp` and exit status 0, without running
any FAM test, so a wrapper reading only the exit code saw green. Root cause:
serial `stp` ended in a bare `stop` (exit 0) while the MPI branch already
aborted with code 1 - `stp` now stops with code 1 on both branches, which
makes every abort loud to wrappers and CI. The one deliberate stp-as-halt
(test_RP_commutator's 'end of test', caller-less) keeps exit 0 through a
plain stop with identical stdout. Verified 2026-11-19: default t0t3 run of
fam_unit_tests.sh NLO now exits 1 with "ERROR: FAM calculation failed";
the SLy4 run keeps overall exit 0 with the full fam net green. The wrapper
workaround (--parameterisation SLy4 for NLO fam builds) stays valid.

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
### F-10 Unguarded namelist reads bypassed find_nml_error — FIXED (cbd45ac6)

Found 2026-11-19 during the F-6 work: only 3 of 17 input-namelist reads
(nucleus, skf, MomentConstraint) were guarded, so a typo in any other
namelist aborted inside the compiler runtime (exit 2, raw backtrace, no
indication of the offending input), and the one guarded call in geninfo
evaluated int(file_number) on an absent optional. All reads now carry
iostat+iomsg and route to find_nml_error, which prints the runtime
message (it names the unknown object) plus the offending line when the
input is a repositionable file; from STDIN the message stands alone.
Verified on NLO-P mf and NLO fam: unknown variables in &mesh (stdin and
file mode), in &skf and a truncated &pairing all report and exit 1, the
success path is unchanged, fam unit net (SLy4) green. The gennilsson
standalone tool's bare /nil read is out of the mocca.data path and left
as is.

### F-11 Pasta serial runs armed the fission inertia pipeline by default — FIXED

Found 2026-11-20 through the red boundary_conditions and spaghetti CI jobs.
Every serial run defaults N_inertia=4, whose hardcoded moment list
Q20/Q22/Q30/Q32 reaches l=3 - one past the multipoles a pasta build carries
(MaxMoment=2, moments.f90). The pre-sprint legacy Qlm_spme built those
out-of-scope inertia tables silently, so the pasta jobs were green; the
Qlm_spme_matrix guard (design 13, sprint step 3) aborts on them and, since
stp exits nonzero (F-6), pasta runs now die after the final printout but
before writewavefunction - the spaghetti restart then failed on the missing
tant_cyl.hdf5. Fix: pasta builds mirror the MPI default, N_inertia=0;
collective inertia is opt-in there (&IO/ with more than 4 moments and an
&inertia/ list within MaxMoment, validated loudly by read_inertia and by a
new ReadIOInput check rejecting 1..4). functional.f90's own
calc_collective_inertia call was already PASTA==0-gated - the end-of-run
block was the last unguarded consumer. Verified on BXL gfortran serial
builds: spag.sh (start/restart energies -489.574811/-489.574821),
boundary_conditions.sh (-130.858108 antiperiodic vs -130.858104 periodic)
and lattice_energy.sh all exit 0; a pasta input with N_inertia=4 stps at
input time. The same fix moved spag.sh's mocca_check1 capture off the
echo - it had masked the start-run abort as green, leaving only the
restart path to catch the failure.
