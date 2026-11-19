# AGENTS.md — MOCCa

Working notes for AI coding agents in this repository. Everything below was
verified against the code and tooling on 2026-11-16, on branch `matrix_refactor`
at commit 9960e8b2. If the file disagrees with the code, the code wins — update this file.

MOCCa is a 3D mesh solver for Skyrme mean-field and linear response physics.
There is no single "source code of the solver": the Fortran sources in
`src_orig/` are templates that Hephaestos specialises per build configuration.
User-facing documentation lives in `README.md` and `docs/`; this file is for
agents who edit and build the code.

## Maths and inconsistencies
 - This is a challenging refactor that blends a lot of coding and mathematics/physics. 
 - Do not hesitate to point out inconsistencies / errors in my prompts/conventions/...
   or to ask for further clarification of the user.
 - Do NOT try to rederive conventions or guess the precise meaning of operators. Ask the user.
## Editable vs generated

| Path | What it is | Edit? |
|---|---|---|
| `src_orig/*.f90` | The true sources, with `$`-markers for Hephaestos | YES — always here |
| `src_heph/*.py`, `src_heph/fortran_templates/*.py` | The Hephaestos generator itself | YES, but only for generator changes |
| `configs/*.py` | Build configurations (select symmetries, functional, ...). One file = one executable name | YES |
| `functionals/*.func` | EDF parameter sets referenced by configs | YES |
| `Makefile`, `make_include/make.inc.*` | Build machinery and machine templates | YES |
| `testing/**`, `scripts/**` | Integration/unit tests, hooks | YES |
| `build/<CONFIG>/{src_mf,src_fam}/` | Hephaestos output, per CONFIG | NO — generated |
| `build/<CONFIG>/{obj_*,mod_*}/` | Objects and module files | NO — generated |
| `exec/*.exe` | Compiled executables | NO — binary output |

`make clean` removes everything generated under `build/` (`.f90`, `.o`,
`.mod`) and it is regenerated on the next build. Nothing in `exec/` is
tracked by git (see `.gitignore`), nor is the local `make.inc`.

## Hephaestos and the `$`-markers

The Makefile generates each variant with

```
python3 Hephaestos.py <CONFIG> mf|fam <DENSUM> build/<CONFIG>/src_mf|src_fam <file>|dry-run
```

triggered automatically whenever any of these change: the matching
`src_orig/*.f90` file, `configs/<CONFIG>.py`, any `src_heph/*.py` /
`Hephaestos.py`, or the functional file referenced by the config.

`src_orig/` code carries `$`-markers, e.g. `$TR` (230 uses), `$NTR` (610),
`$PBROKEN` (155), `$OMP` (71), `$FAM` (43), `$REDUX/Y/Z`, `$N2*`/`$N3*`,
`$TAUSCALAR`/`$TAUTENSOR`. They gate code per symmetry scenario and
functional order: `$TR` regions only exist in time-reversal-conserving
builds, `$NTR` only in TR-broken ones, `$FAM` only in the FAM executable,
and so on.

Consequences when editing near markers:
- A change that compiles and passes tests for one CONFIG does not prove
  the others. The variants that matter most are kept/broken time-reversal
  (`NLO` vs `NLO-T`) and parity (`NLO-P`, `NLO-TP`); build more than one
  CONFIG when touching marked regions.
- Analyze what a marker activates with the dry-run: it processes all files
  and prints Hephaestos' report. Marker semantics live in
  `src_heph/heph_substitute.py` and `src_heph/heph_symmetries.py`.

## Building

- `make.inc` (root, gitignored, machine-local) holds the toolchain: `CXX`
  (gfortran here), `OPTFLAGS`, `USE_MPI`, `USE_HDF5`, `HDF5_LIB`,
  `LINEAR_ALGEBRA_LIB`, directory locations. Templates for other machines
  live in `make_include/make.inc.*`; `make INCLUDE=gnu-serial` builds
  without a root `make.inc`.
- Accelerate compilation when testing by using `make INCLUDE=gnu-debug`,
  the resulting executable will be slow.
- `make mf` and `make fam` build the mean-field and FAM executables;
  `make` defaults to `mf` and `CONFIG=default`.
- `make CONFIG=NLO mf` uses `configs/NLO.py` and produces
  `exec/MOCCa.NLO.exe`; the FAM variant produces `fam.NLO.exe`.
- Gotcha: `make.inc` is not part of the dependency graph. Changing
  optimisation flags or USE_MPI does NOT recompile anything. Run
  `make clean` (or touch the sources) after touching `make.inc`, or you
  will compare stale binaries against your edits.
- The version banner printed at runtime (commit, compiler, flags) is
  injected by the Makefile from `git show` and the compiler version.

## Style
  
- Formatting: fprettify with `.fprettify.rc` — indent 2, line length 132.
  Check: `fprettify -d -r -c .fprettify.rc src_orig/`; fix with `-r`.
- fprettify-version caveat (observed 2026-11-17): the check's output
  depends on the installed fprettify version. Recent versions impose
  normalization rules (end do / end if, double-quoted format strings)
  that the committed sources predate — their check flags essentially
  all of `src_orig/` (and fatally errors on `functional.f90` line 832),
  and running fix mode in-place reformats whole files beyond intent.
  The house style is the committed one: `enddo`, `endif`, single-quoted
  formats, indent 2, 132 columns. Match it; keep diffs minimal. CI's
  format step only produces a report artifact, it does not fail builds.
- Linting: fortitude with `.fortitude.toml` (currently only S201,
  superfluous implicit none, and C001).
  Check: `fortitude --config-file .fortitude.toml check src_orig/`.
- House comment style: `#`-fenced banners and `- - - -` separators, as
  found throughout the sources and `testing/`.
- Functions and subroutines should at least have a documentation section 
  explaining input and output.
- An example of the comment type expected:
  `subroutine set_spwf_symmetries(sx,sy,sz,blocks)
    !---------------------------------------------------------------------------
    ! Routine that assigns the correct "signs" under symmetry transformation 
    ! to all single-particle wavefunctions. 
    !
    ! Hephaestos fills in all practical signs here, using the following key
    !
    !  [dollar sign] SXab
    !
    !  where 
    !      X  = X;Y;Z depending on the direction
    !      a  = block index, i.e. 1-4 depending on the symmetries of the spwf
    !      b  = component index, i.e. 1-4 depending on the spinor component
    !           we are dealing with
    !---------------------------------------------------------------------------
    ! Input: blocks 
    !        The number of spwfs in every symmetry block. 
    !        This is an input of this routine, because this subroutine could be 
    !        called for a subset of spwfs when reading from file.
    !
    ! Output: sx, sy, sz
    !        The signs with respect to x/y/z reflection of the components of 
    !        the spwf spinors. 
    !---------------------------------------------------------------------------`

## Language 
- You should aim to be concise yet precise. This is particularly important 
  for important for in-code documentation.
- Where possible, use nuclear physics terminology - not terminology borrowed 
  from other fields. Example: "good quantum numbers" instead of "gradings". 
  
## Testing

All scripts are run from inside their containing directory:

```
cd testing/integration && bash minimal.sh NLO    # run one test
cd testing/integration && bash golden_diff.sh minimal NLO   # run + golden compare
```

- Integration tests (`testing/integration/*.sh`) follow the conventions of
  `testing/functions.sh`: an executable suffix argument; a transient
  `work/` directory destroyed at teardown; STDOUT logs written to
  `logs/<test>.<suffix>.out`. The `logs/` tree is gitignored except
  `*.reference.out`.
- Unit tests: `testing/unit/fam_unit_tests.sh` and
  `testing/unit/moments_unit_tests.sh` (the moments suite also runs in CI:
  the `moments_unit_tests` job of `.github/workflows/CI.yml`, on the
  gfortran NLO/NLO-T/NLO-P artifacts; the script takes the CI suffix scheme
  through an optional three-suffix interface,
  `moments_unit_tests.sh SUFFIX [SUFFIX_T SUFFIX_P]` - locally the -T/-P
  suffixes derive from the base).
  Also `testing/unit/block_matrix_unit_tests.sh`: the standalone unit
  suite of block_matrix_module itself (compiles the module on the fly,
  no MOCCa executable needed) including the 8 refusal scenarios of the
  rotation structure guards; its CI job `block_matrix_unit_tests` needs
  nothing but lint-and-format, so it gives the fastest feedback of the
  whole workflow.
- CI (`.github/workflows/CI.yml`) compiles gfortran and ifx builds in
  `ghcr.io/iaa-ulb/mocca-ci` containers and runs `minimal.sh` plus the
  integration batch scripts.

### Golden regression harness

`golden_diff.sh [TEST] [WORKSUFFIX] [REFSUFFIX]` reruns an integration test
with the working executable and compares its STDOUT, token by token and
number by number, against `logs/<TEST>.<REFSUFFIX>.out` — stripping the
volatile version banner and timers table, and allowing tolerances
`|a-b| <= atol + rtol*max(|a|,|b|)` (defaults 1e-8/1e-7) for refactor-induced
reassociation drift. `--no-run` re-diffs an existing log, `--update` promotes
the current log to the golden reference deliberately.

Golden references are ordinary test logs produced by a trusted executable
under the `<SUFFIX>.reference` convention, e.g.
`bash minimal.sh NLO.reference` — and they are committed (that is the point
of the `.gitignore` exception). Never silence a golden failure by widening
tolerances without understanding the diff it reports; the harness prints the
largest deviation and its line for exactly that.

Locale warning: on machines with a non-C `LC_NUMERIC` (e.g. `nl_BE.UTF-8`),
mawk truncates string-to-number conversion at the decimal point. Any awk
that parses numbers out of MOCCa output must be invoked as `LC_ALL=C awk`.
`golden_diff.sh` does this; keep doing it in new scripts.

## MPI

- Serial by default here (`USE_MPI=0` in `make.inc`); enabling MPI uses the
  same sources (`-DUSE_MPI` for the preprocessor) with an MPI-capable
  compiler.
- Data model: single-particle wavefunctions are rank-local; see
  `rank_map`, `spwf_map`, `spwf_inverse` in `src_orig/wavefunctions.f90`.
  Global quantities are merged over ranks.
- New code should stay MPI-plausible: guard communication through the
  existing patterns; do not introduce bare rank-dependent logic.
  Testing validity with MPI remains the developer's explicit choice.

## Module map

44 Fortran modules over 46 source files in `src_orig/`. Regenerate with:

From the repository root:

```
awk '/^[ \t]*module [A-Za-z0-9_]+/ && $2 != "procedure" {print $2" <- "FILENAME}' src_orig/*.f90 | sort
```

Grouped (module <- file where both are obvious):

- Drivers: MOCCa <- mocca.f90, SCFiteration, convergence, evolution,
  run_single.f90, fam (fam.f90), fam_run.f90, fam_testing, gmres <-
  fam_gmres.f90
- Numerics and mesh: constants, vectors, preconditioning <-
  precondition.f90, derivatives, folding, sphericalharmonics, nil8
- Wavefunctions and bases: wavefunctions, basis_transform, transform,
  hartreefock <- hartree-fock.f90, BCS, HFB, HFB_direct, HFB_gradient
- Symmetry-block matrices: block_matrix_module <- block_matrix_module.f90,
  the matrix_refactor sprint (self-contained: plain Fortran 2003+, no
  $-markers, no module dependencies)
- Densities, functional, pairing: densities, functional, parameterization,
  pairing, pairingcutoffs, pairing_strengths, Coulombmod <- coulomb.f90,
  cranking, moments <- moments.f90, momentsofinertia, fission_MOI
- I/O and infrastructure: IO, IO_aux, IO_wf, HDF5_auxiliary, Printing <-
  printing.f90, GenInfo <- geninfo.f90, compilation, timing, version

The ME / expectation-value territory of the code: `wavefunctions.f90`
(spwf storage, symmetry-block bookkeeping, the angmom_* / spin_* ME
routines), `momentsofinertia.f90` (J_mu matrices and inertia tensors),
`moments.f90` (multipole expectation values), `functional.f90`
(incl. CompCOMCorrection), `cranking.f90`, `basis_transform.f90`
(HF <-> canonical, HFBlocks usage), and `fam.f90` (drho, dkappa, FAM
perturbations).

## Current focus (updated 2026-11-19 — prune when stale)

- Refactor sprint on branch `matrix_refactor`: introducing a dense
  symmetry-block matrix type for efficient matrix-element calculation —
  at most 8 blocks per operator, selection rules resolved by one XOR at
  construction, BLAS-only runtime. The authoritative design is
  `mocca-symmetry-block-matrix-design-summary.md` at the repo root — not
  this file; read it first.
- Sprint state (updated 2026-11-19, end of sprint day). Steps 1-3 of the
  design's §10 - and all five rows of §13.3 - are COMPLETE on top of
  phases 1-2; every numerical gate is executed and green:
   * Phases 1-2: core sym_matrix types, init with goodqn mask (4 bits,
     $GOODQN), the classifiability interface (mandatory herm x sty
     classes, derived kind), fill machinery (block_filler callback,
     element_filler fallback, new_sym_matrix sugar). Commits 395bde9d,
     6a6d92c8, 55272151, 35bf016e. Unit suite:
     testing/unit/block_matrix_unit_tests (33506 checks at the time;
     53562 plus the 8-scenario refusal program after design 6a
     landed, see the sprint state below).
   * Step 3, the fm multipole generator: the verbatim Qlm_spme move
     (9faf8238), Qlm_spme_matrix/Qlm_spme_fm - BOTH halves: r^l Re Y_lm
     STY_EVEN/kind symmetric, r^l Im Y_lm STY_ODD/kind antisymmetric,
     chamber-weighted fill (810986ab, d7f42ff2), the dense() helper and
     the off-by-default spme_oracle hook (59ae4252, whose NLO-T/P/TP
     smoke runs and the parity keys predate 93e6e8de), the inertia
     switch (c0262789) and the AMENDED units flip (0b6de03f: the
     inertia chain is fm-based internally, every presented number stays
     in the historical b^(l/2) units - design §13.3 row 5 carries the
     amendment).
   * Gates all green: golden minimal NLO bit-identical at every row
     (4849 numbers); oracle agreement ~1e-15 relative (real halves,
     NLO/NLO-T/NLO-P); O16 inertia smoke byte-identical and inertia
     file bit-identical to the legacy binary; the heavy fission.sh
     (Pu240 HFB) validated on CI against the pre-move baseline
     (user-executed, 2026-11-19).
   * Moments test infrastructure: the moments_testing module +
     unit_test flag move to &mesh/ (cd8d9b8e), the three-configuration
     cross-build runscript (7dcc6118, NLO/NLO-T/NLO-P) and its CI job
     (57217a28). Conventions distilled below; usage in the Testing
     section.
  AMENDED STEP 4 LANDED (2026-11-19, 497b1d94): matvec and matmat are
  dropped from the type and the plan (consumer audit recorded in design
  §6/§10 and commit 1230f34a); the full-box reconstruction oracle for BOTH
  halves of Qlm_spme_matrix is in moments_testing (test_Qlm_fullbox:
  spinors extended across the fold planes with the sx/sy/sz component
  signs, operator regenerated at every chamber's coordinates, plain mesh
  quadrature - the chamber-weight table and the oc estimation bypassed;
  judgments against physics: identity incl. the dense() partner
  reconstruction, STy-characteristic vanishing of the opposite half,
  Hermitian conjugate). Gates green on NLO/NLO-T/NLO-P (frozen make.inc):
  identity reassociation-level (max 5.3e-13 / 1.8e-13 / 6.8e-14 absolute),
  opposite halves at <= 6.1e-16 everywhere, the odd-m signature-partner
  imaginary blocks of NLO-T covered at ~1e-14 - the wobbling-relevant
  content that had no reference path until now.
  STEP 5 LANDED (2026-11-19, 20fca2e5): get_qpme_op's two call sites
  feed from the sym_matrix generator - 'multipole' (l>0) from
  Qlm_spme_fm(l,m,.false.), 'zcom' from Qlm_spme_fm(1,0,.false.), in
  fm^l natural units. The b^(l/2)->fm rescales are gone; the sqrt(2)
  (m.ne.0) and zcom prefactors survive at the call sites. Gates green
  (frozen make.inc, serial): fam unit tests SLy4; A/B against the parked
  fam.NLO.reference.exe on a shared O18/SLy4/HFB state - Q20 and Q22
  strength famfiles bit-identical, GMRES residuals identical; fam
  builds of NLO / NLO-T / NLO-P compile. The Im half stays deferred
  (design §13.5).
  PARITY-ODD FAM FIELDS (understood 2026-11-19, user-confirmed): in
  parity-CONSERVING builds (NLO) the FAM machinery cannot represent a
  parity-breaking excitation - every parity-odd field NaNs (l=1
  multipole incl. the committed fam_spurious_mode_subtraction scenario
  and operator_type='zcom'), identically in the parked pre-refactor
  anchor: by construction, not a regression. Making such responses
  available is the long-term goal of this refactoring. For now the
  working scenario is the PARITY-BROKEN build; NLO-P has no parked fam
  reference, so the zcom call site's net is source-level (moments
  full-box oracle, Qlm_spme_fm(1,0), all builds) plus the NLO-P runs
  below. test_RP_commutator (fam_testing, zcom/zmomentum chain) exists
  but is commented out in fam_run.f90.
  NLO-P FAM VALIDATION (2026-11-19, post-step-5 fam.NLO-P.exe, frozen
  make.inc): O18/SLy4/HFB, nx=ny=8 dx=0.8, nz=16 (unfolded z), even
  nw=14, allowtransform. operator_type='zcom' and l=1 m=0 both: GMRES
  converges (17 iterations, residual 1.3e-08), ||b|| finite, strength
  tables internally consistent (S = -1/pi Im S_complex), <[R,P]> purely
  imaginary (-0.81 i; sign convention documented in the code, magnitude
  vs ideal 1i plausibly truncated-basis commutator inexactness), and
  the spurious-mode subtraction drives S_after to ~1e-17 for BOTH
  fields (isoscalar dipole = pure COM mode) - the fam_spurious_mode_
  subtraction PASS criterion holds on NLO-P. Known cosmetics: the
  strength table's parity + / - rows are vestigial in a parity-broken
  build (all strength books into parity +), and the box-size check's
  drho_antisym panel prints NaN exactly as in healthy NLO runs. The
  committed spurious script itself needs its mf input flipped (nz=16)
  before it can run against NLO-P executables.
  CLEANUP LANDED (2026-11-19, 4ca0f181): the last two legacy Qlm_spme
  call sites migrated to Qlm_spme_fm - functional.f90's
  verify_COM_motion and fam_testing.f90's test_RP_commutator. Audit
  corrections recorded in the design doc: the "functional
  (CompCOMCorrection, Q_10)" consumer label was wrong twice over -
  CompCOMCorrection never used Qlm_spme, and verify_COM_motion (the
  actual reference, a pushing-model diagnostic) is DEAD CODE: no call
  site anywhere in sources, markers or templates. The migration
  therefore changed no live numerics: golden minimal NLO 4849/4849
  numbers identical, fam unit tests and the moments unit net green,
  mf+fam compiles green on NLO/NLO-T/NLO-P (all frozen make.inc).
  verify_COM_motion additionally gained an explicit
  `use moments, only: Qlm_spme_fm` - it reached Qlm_spme before only
  via the transitive use fission_MOI -> use moments chain.
  STEP 6a LANDED (2026-11-19, design §6a - basis rotations, the
  intermediate infrastructure step before the §13.5 deferred items):
  the sym_rotation type in block_matrix_module holds one square real-
  orthonormal block per populated sector (the shape MOCCa's own
  producers make: HFB's Canonical diagonalises rho per 8-block, so
  cantransfo is exactly sector-diagonal; basis_transform's merged
  signature-pair loops mix nothing), built three ways: init (identity
  per populated sector), init_from_dense (extracts the sector-diagonal
  blocks from a full n x n rotation array in the basis_transform
  convention, REFUSING off-sector-nonzero or non-orthonormal input),
  set_sector (direct per-sector handover, orthonormality checked;
  tolerance default 1e-12 absolute on both paths). The in-place rotate
  of a sym_matrix transforms every stored block independently,
  A_b <- U_{srow}^T A_b U_{scol} (two matmuls per block, intrinsic
  matmul - the module stays free of explicit BLAS calls so the
  standalone unit suite links nothing), preserving every metadata
  field (herm/sty/kind/symmetry word/goodqn/block list/st): the same
  operator in the rotated basis, so the dropped-mat-mat decision stands
  - a similarity transform is not an operator product. Scope stands
  (user-confirmed): sector-diagonal real orthogonal ONLY - complex
  unitaries fall out of the STy-invariant model (STy antilinear) and
  signature-mixing rotations would need a signature-masked basis no
  config produces. No consumer migrated (first production candidates:
  the U^T O U factor of get_qpme_op's Bogoliubov image, canonical-
  basis Ksum/MOI chains, transform_mat call sites). Gates green
  (frozen make.inc): block_matrix unit net 53562 checks, 0 failures
  (was 33506), plus the new 8-scenario refusal program
  testing/unit/block_matrix_refusal_tests.f90 driven per exit status
  by block_matrix_unit_tests.sh (the structure guards stop, so they
  cannot be caught in-process); golden minimal NLO 4849/4849 numbers
  bit-identical; moments suite green on NLO/NLO-T/NLO-P incl. the
  cross-build comparison; mf builds of all four configs compile
  (NLO-TP had to go serial - the -j4 run races on missing .mod files).
  ANGULAR-MOMENTUM GENERATORS LANDED (2026-11-19, design §14.4 row 2;
  row 1 - the K_x fill fix - CANCELLED by measurement, §14.2): the
  marker-free module src_orig/angular_momentum.f90 hosts the five
  generators on sym_matrix - honest J_x (word 1, HERMITIAN x
  STY_EVEN), J_y (word 1, HERMITIAN x STY_ODD), J_z (word 0,
  HERMITIAN x STY_EVEN) and the TR-folded twins J_x.T (word 0,
  HERMITIAN x STY_EVEN) and J_y.T (word 0, ANTIHERMITIAN x STY_ODD,
  both halves kind +1 SYMMETRIC, measured) - wrapping
  wavefunctions.f90's angmom kernels verbatim as per-pair values. The
  folded generators exist in every build; callers gate ($TR-marked
  call sites only), the module is Hephaestos-marker-free. The arm
  src_orig/angular_momentum_testing.f90 (run_angular_momentum_tests,
  armed by the unit_test flag from Run_MOCCa, before
  run_moments_tests) is the kernel-transpose oracle: per object it
  verifies the block structure and stored values, measures
  K(wf_i,wf_j) vs K(wf_j,wf_i) per stored pair with the raw kernels and
  checks agree with the tagged kind, checks partner consistency
  A(s',s) = kind*A(s,s')^T, and PRINTS (does not judge) the purity of
  the opposite-half kernels. Arm verdicts on the O16/SLy4 unit state
  (maxiter=0): NLO - honest J_x/J_y structural zeros (0 blocks, the
  odd-m precedent), J_z/J_x.T/J_y.T symmetric at 2.5e-14/1.6e-14/
  1.1e-14; NLO-T - all five agree (J_y antisymmetric at 2.1e-15, both
  folds symmetric at 7.5e-15/2.1e-15). NLO-P via the frozen -O3 build
  only (see the coulomb open item). The §14.2 audit: the suspected
  K_x fill bug is REFUTED - the symmetric `$TR` fill is correct, an
  antisymmetric-fill A/B with symmetries.sh breaks exactly the
  x-consuming rows (Belyaev X 0.9295 -> 3.2213, the two <J^2>-<J>^2
  definitions split on X), source reverted, no production line
  changed.
  NEXT UP: the calcJ2andBelyaev_HFB consumer migration (design §14.4
  row 3), feeding the sym_rotation consumers of §6a. After that, the
  §13.5 deferred items; alternatively, first dispose of the single
  remaining Qlm_spme reference.
  OPEN ITEMS when resuming:
   * DEAD angmom KERNELS / THE jz-IS-REAL ERROR (found 2026-11-19 by
     the arm's purity print, recorded PRECISELY in design §14.2b):
     angmom_x_imag, angmom_y_real, angmom_z_imag are referenced
     nowhere in production (the unit arm is their first caller) and
     all three measure O(1) on the reduced mesh (same order as the
     validated halves), in NLO/NLO-T/NLO-P alike. This CONTRADICTS
     the code's own completeness premise - momentsofinertia.f90:774/
     :792: "jz_ij is always a real number", the §2a corollary that
     the stored real half is the complete operator. Two candidates,
     unadjudicated: (a) reduced-mesh artifact - the purity print
     integrates one chamber with dv multiplicity and the
     opposite-half cross-parity products would cancel across the fold
     planes of the full box (the multipole chamber-weight-0 case of
     design §13.1), so the true Im<jz> IS ~0 and only the naive
     single-chamber integral is O(1); (b) a genuine bug or packing
     mismatch in the dead kernel, which would put the premise in
     doubt. Deciding test (not built): full-box reconstruction of
     Im<i|J_z|j> per the §13.1 step-4 pattern. Nothing in production
     is affected either way; deciding (and then possibly deleting
     the dead kernels) is a user decision. (The symmetries.*.out
     logs have been rerun green on the frozen executables after the
     K_x A/B; the broken-experiment and healthy pre-change evidence
     of design §14.2 lives only in the doc's recorded numbers now.)
   * Qlm_spme deletion gate (design §13.2), now purely a DECISION: the
     routine's only remaining reference is the off-by-default spme
     oracle hook in fission_MOI.f90. Either Qlm_spme and the oracle go
     together, or both stay permanently as the production-scale cross-
     check of Qlm_spme_matrix. The moments unit net covers the
     generator either way.
   * verify_COM_motion and test_RP_commutator remain call-less dead
     code; enabling either (e.g. RP as an NLO-P zcom net) is a user
     decision, out of the migration's scope.
   * LATENT BUG (pre-existing, found 2026-11-19, NOT a sprint regression):
     the debug-profile (-fbounds-check) NLO-P binary crashes in
     coulomb.f90's obtain_folded_potentials (fold_form_factor_reverse,
     index 26 of dimension 3 vs extent 24) during set_up_mean_field -
     reproduced 2026-11-19 by the angular-momentum arm attempt (the
     debug moments_unit_tests NLO-P leg terminated on exactly this
     before the arm could run; the NLO-P arm gate is the frozen -O3
     build). An out-of-bounds in the parity-broken Coulomb folding path
     that the -O3 CI artifacts never check. Do not use debug-profile
     NLO-P builds until it is fixed; the frozen -O3 suite runs green.
   * exec/MOCCa.NLO-TP.exe was rebuilt with the frozen make.inc after
     the rotation commit (2026-11-19) - no longer the stale debug
     binary, but compiled only, no runtime-validated scenario; the
     mf executables NLO/NLO-T/NLO-P were rebuilt with the frozen
     make.inc on 2026-11-19 (including the oracle) and again after
     the rotation commit. fam.NLO predates the oracle commit.
- Conventions learned in this sprint (2026-11-19, user-confirmed; full
  reasoning in the commit messages and the design doc):
   * Folded meshes: the conservative/TR builds carry the z-fold in dv,
     so nz=12 describes a 24 fm z box; the parity-broken build stores
     the z-axis UNFOLDED and needs nz=24 for the same box (nz=12
     converges to a different truncated-box state, ~88 keV off for
     O16/SLy4 - a broken input, not new physics). Broken-symmetry
     inputs also demand even nwn/nwp and allowtransform=.true. in &IO/.
   * One-body Q^2 through sp matrices presumes the truncated basis
     RESOLVES Q|i> - false for any MOCCa basis (failure ~1e4 relative,
     worst for high l): no Qlm^2-vs-mesh-Squared test exists.
   * Odd-m multipole operators are structurally zero in TR-conserving
     builds (all stored sectors have q(Rz)=+i, their XOR partners stay
     empty); their wobbling-relevant content exists only once the
     TR-broken blocks are populated.
   * fam_unit_tests.sh with its default --parameterisation t0t3 stops
     silently at the .func mismatch under an NLO fam build (stp, exit
     0) and never runs the FAM tests - pass --parameterisation SLy4.
   * Cross-build comparisons are only meaningful through
     basis-independent invariants (traces) of operators common to the
     builds (even m), and only when the runs act on the SAME state: the
     moments unit suite (2026-11-19) compares the three builds on their
     common INITIAlisation with zero SCF iterations (maxiter=0). The
     CI failure that motivated this: each build stops at a slightly
     different point of the near-flat convergence landscape, and the
     parity-broken build's l=4 box-anisotropy rows showed converged
     traces at ~8e-9 local but ~3e-6 in CI - environment-sensitive by
     construction, not a machinery error (within each run, matrix trace
     and density-route moment agreed to 1e-13 in every environment).
   * The unit_test flag lives in &mesh/ (geninfo) and arms every
     executable's own suite (moments: run_moments_tests from
     Run_MOCCa; fam: run_FAM_tests from fam_run).
   * MOCCa's time reversal is FERMIONIC: T**2 = -1 on the spwfs
     (user-confirmed 2026-11-19; design §14). Consequences, MEASURED
     2026-11-19 by the kernel-transpose oracle (design §14.1/§14.2;
     the earlier asserted classes are superseded): the antilinear
     composite J_mu.T is HERMITIAN as an operator - the user-accepted
     chain (J.T)^dagger = -T**2 (J.T) = +(J.T) - and BOTH folded twins
     measure kind +1 SYMMETRIC storage; K_x = J_x.T is word 0,
     HERMITIAN x STY_EVEN; K_y = J_y.T is word 0, STY_ODD with herm =
     ANTIHERMITIAN as the kind-COMPLETING label (herm = kind * sty) -
     not a physical class: the §2a (herm,sty)->kind machinery does not
     transfer to antilinear composites. All kind signs are test
     assertions settled by the arm's oracle, never re-derived.

- Oracle-run recipes for the hook (local, gitignored): copy the
  parameterization (e.g. SLy4.param) next to the executable; TR-broken
  and parity-broken builds demand even nwn/nwp and
  `allowtransform=.true.` in &IO/; arm with `N_inertia>=5` in &IO/,
  an &inertia/ namelist (first four moments hardcoded Q20/Q22/Q30/Q32,
  worth adding Q40) and `spme_oracle=.true.`. Production-like NLO-T
  scenario: Mg24 HFB pairingscheme=1, omegaZ=0.3292333881 (cranking.sh
  style; pure HF inputs force diagsphamil on, which reshuffles the
  stored basis within blocks mid-run - judge oracle numbers on
  production-like inputs only). Caveat: with spme_oracle armed the
  hook runs per calc_collective_inertia call, which fires every
  iteration while unconverged - keep maxiter small for dev scratch
  runs.
- Local test artifacts for step 3 (gitignored, machine-local):
  `exec/fam.NLO.reference.exe` is the parked FAM anchor;
  exec/MOCCa.NLO.exe, MOCCa.NLO-T.exe, MOCCa.NLO-P.exe and fam.NLO.exe
  are freshly rebuilt and validated (frozen make.inc; fam.NLO
  includes the step-5 migration); fam.NLO-P.exe is validated on the
  zcom/l=1 NLO-P scenarios above; fam.NLO-T.exe is compiled (2026-11-19)
  but has no runtime-validated scenario yet;
  the mf and fam executables of NLO/NLO-T/NLO-P and the mf executable
  of NLO-TP were rebuilt with the frozen make.inc after the
  angular-momentum commit (2026-11-19, including the arm): NLO is
  golden-validated (minimal 4849/4849) and moments-suite-validated
  with the angular-momentum arm green on all of NLO/NLO-T/NLO-P;
  NLO-TP is compiled only (the debug-profile latent bug in the
  parity-broken Coulomb folding path is untouched, see the open
  items); fam.NLO is fam-unit-validated (SLy4) on the rebuilt binary;
  fam.NLO-T/fam.NLO-P are compiled only, without a
  runtime-validated scenario.
  fam_unit_tests NLO/SLy4 is the runtime green net for fam.NLO; the
  multipole A/B scratch scenarios live in /tmp only (O18/SLy4/HFB,
  Q20 and Q22 @ 25 MeV vs the parked anchor, bit-identical famfiles).
  `testing/integration/logs/fission.BXL.pre-move.out` was the inertia
  baseline produced by the pre-move BXL binary; the fission.sh comparison
  against it has been executed through CI (2026-11-19) - the batching
  recipe is retired.
  fission.sh (Pu240, HFB) is too slow to run per commit.
- Sprint workflow (2026-11-17): quick compile iterations use the debug
  profile - `make CONFIG=X clean` then `make CONFIG=X mf INCLUDE=gnu-debug`
  (-O0, bounds checks, no HDF5; none of this step's test nets needs
  HDF5). The root make.inc stays frozen (-O3 + HDF5) and must back every
  numerical gate or golden comparison; remember make.inc is outside the
  dependency graph, so an INCLUDE= switch without clean yields mixed
  binaries silently.
- Design-to-repo couplings: the module lives in `src_orig/` as plain
  Fortran 2003+ without `$`-markers, identical for every CONFIG; its
  constructor is the only place converting MOCCa's 1-based `HFBlocks`
  ordering to the 0-based sector word, and it must take global counts
  (`HFBlocks` is rank-local; `HFBlocks_global` is the global twin).
  Agreed encoding (2026-11-16): bit 0 = z-signature sign (0 = q(Rz) = +i),
  bit 1 = parity (0 = even), bit 2 = isospin (0 = neutron), so the sector
  word is simply the MOCCa block index minus one; see design summary §2a.
  Operator words handed to `init` are physical; the constructor additionally
  takes the good-quantum-number mask (design §2b/§2c, a 4-bit word: bits 0-2
  as the sector words, bit 3 = STy conserved), injected by Hephaestos into
  geninfo as $GOODQN, and the mandatory physical operator classes `herm`
  (HERMITIAN/ANTIHERMITIAN) and `sty` (STY_EVEN/STY_ODD) of design §2a,
  from which the storage kind is derived as kind = herm*sty — no
  `$`-marker selection logic in the routines that produce or consume
  sym_matrix objects.
- Hephaestos wiring: a new `src_orig/*.f90` file must be added both to
  `MF_SRC`/`FAM_SRC` in the Makefile and to `no_process_files` in
  `src_heph/preprocess.py`, or the per-file generation silently skips it
  (empty output file, compile failure).
- Storage is single-node by design (summary §1 and §12), not just a
  sprint convenience: keep the type free of rank-dependent logic and
  leave any gathering to whole-basis vectors at call sites.
- No matvec/matmat in sym_matrix (2026-11-19, user-confirmed after a
  consumer audit of the tree): nothing applies a one-body operator
  matrix to a basis vector - sp ME matrices are only traced against
  densities, contracted elementwise with other operator or amplitude
  matrices, or fed through within-block U/V basis transforms; fields act
  on wavefunctions on the mesh through the density route. Design §6/§10
  carry the amendment; the Qlm_spme deletion gate is redefined in §13.2.
- BLAS kernels (`dgemv`/`dgemm`) come from the existing
  `LINEAR_ALGEBRA_LIB` (`-llapack -lblas`); the summary's OpenBLAS/MKL
  suggestion (its §7) is deferred while `make.inc` stays frozen.
- The trusted reference executable is `exec/MOCCa.NLO.reference.exe`
  (serial, gfortran 13.3, `-O3 -fbacktrace -g`, HDF5 on); the golden
  `testing/integration/logs/minimal.NLO.reference.out` was generated from
  the un-refactored code at commit b719453e. Keep `make.inc` frozen while
  it is the trust
  anchor; changing flags invalidates golden comparability.
- The sprint is explicitly serial-only: MPI behaviour does not have to be
  preserved or validated during it.
- FAM anchor: `fam.NLO.reference.exe` is parked (built from the un-refactored
  tree, 2026-11-17) — the FAM-side twin of MOCCa.NLO.reference.exe and step
  0 of the design's §13. Design step 3 (§10/§13) migrates fam.f90's
  Qlm_spme consumer in its final commit and design step 9 targets FAM;
  FAM-side regressions are judged against the parked executable (no
  committed fam golden log exists, exec/ is untracked).
