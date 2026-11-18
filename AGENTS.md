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
- Unit tests: `testing/unit/fam_unit_tests.sh`.
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

## Current focus (dated 2026-11-17 — prune when stale)

- Refactor sprint on branch `matrix_refactor`: introducing a dense
  symmetry-block matrix type for efficient matrix-element calculation —
  at most 8 blocks per operator, selection rules resolved by one XOR at
  construction, BLAS-only runtime. The authoritative design is
  `mocca-symmetry-block-matrix-design-summary.md` at the repo root — not
  this file; read it first.
- Sprint status (implementation order in the design's §11): steps 1-2
  landed — core types and `init` (phase 1, commit 395bde9d), and the
  fill machinery (phase 2): `fill_blocks` with the `block_filler`
  abstract interface, the `element_filler` per-element fallback
  (`fill_elements`), and the `new_sym_matrix` constructor sugar. No
  `add_element` — the element-wise path is the callback only (see
  design summary §5). Unit tests: `testing/unit/block_matrix_unit_tests`.
  Next: step 4 (mat-vec) after the step-3 chain completes - see the
  re-plan below and design §§11/§13.3; step 3's remaining rows are the
  inertia switch (row 4), the units flip (row 5) and the consumer
  migration (row 6).
  Landed since (2026-11-17): the goodqn mask (commit 6a6d92c8: optional
  `goodqn` in `init`, `$GOODQN` injected into geninfo, masked unit tests),
  the pure, byte-for-byte move of `Qlm_spme` into moments.f90 with the
  ProcessMoments parity keys (commit 9faf8238; golden minimal NLO
  bit-identical), and the generator of design §13.3 row 3: the
  partner-reconstructing dense expansion `sym_matrix%dense()` in
  block_matrix_module (print_dense refactored onto it, unit-tested over
  patterns x words x kinds x masks), `Qlm_spme_matrix` (real symmetric
  half only, fm^l, marker-free, filled through the block_filler callback
  `Qlm_spme_fill_block`), `Qlm_spme_fm` (dense legacy shape) in
  moments.f90, and the off-by-default oracle hook in fission_MOI
  (`spme_oracle` in the &inertia/ namelist; prints legacy-vs-generator
  deviations plus symmetry self-checks). The `me_im` half that briefly
  shipped in row 3 was reverted as a chamber artifact - see the
  convention correction below; it returns at design step 6/7. Oracle
  runs NLO / NLO-T / NLO-P: real halves agree to ~1e-15 relative; all
  four mf configs + fam.NLO compile; golden minimal NLO still
  bit-identical after the commit. Next up: row 4 - the inertia switch
  (`Qsp` from `Qlm_spme_fm * 100^(-l/2)` in calc_collective_inertia).
- Landed since (2026-11-17): the y-time-simplex characteristic of design
  §2a in `block_matrix_module`, with `goodqn` extended to 4 bits (bit 3 =
  STy conserved, injected by Hephaestos via symdic['STy'] alongside Rz
  and P, default 15; NLO/NLO-T 15, NLO-P/NLO-TP 13).
- Landed since (2026-11-18): classifiability interface of design §2a -
  every operator is (anti-)Hermitian and STy-even/odd, so `init` and
  `new_sym_matrix` take the MANDATORY physical classes `herm`
  (HERMITIAN/ANTIHERMITIAN) and `sty` (STY_EVEN/STY_ODD) and DERIVE the
  storage kind as kind = herm*sty; `kind` is no longer an input and
  KIND_GENERAL/STY_NONE are gone. STy carries no block structure: the
  classes select the stored HALF of the complex matrix elements (real
  part even / imaginary part odd), never a block. goodqn bit 3 is not a
  gate - the mask is stored verbatim as the consumers' completeness
  criterion (matrix is the entire operator iff bit 3 set; only relevant
  for hypothetical SYMSTRINGs without STy). `Qlm_spme_matrix` now calls
  `init(HFBlocks_global, iop, HERMITIAN, STY_EVEN, goodqn=goodqn)`.
  Gates: unit suite reworked over the herm x sty cross products (33506
  checks, 0 failures), NLO + NLO-T mf and fam.NLO built with the frozen
  make.inc, golden minimal NLO bit-identical (4849 numbers identical).
- Landed since (2026-11-19): corrected classes and imaginary half of
  Qlm_spme_matrix (user instruction; supersedes the row-3 revert): MOCCa
  works with the real and imaginary multipole operators, both Hermitian,
  of opposite STy class - r^l Re Y_lm is STY_EVEN (stored numbers =
  Re <i|Q|j>, kind KIND_SYMMETRIC), r^l Im Y_lm is STY_ODD (stored
  numbers = Im <i|Q|j>, kind KIND_ANTISYMMETRIC), for every m. The fill
  callback is now a chamber-weighted sum over the folded mesh: per
  spinor-component pair the exact 0-or-1 weight
  prod_folded((1 + s(k,i)s(k',j) oc)/2) built from the sx/sy/sz extension
  signs and the operator's reflection characters oc (estimated numerically
  from mirrored-mesh harmonic generations in Qlm_spme_operator_characters).
  For every XOR-allowed block pair all weights are 1 - the real half stays
  bit-identical to the legacy fill and the oracle agreement (~1e-15
  relative) survives (NLO smoke run armed with spme_oracle); a weight 0 is
  the per-amplitude selection rule the sector word cannot express. NOTE:
  the 2026-11-17 revert message declared the within-block imaginary half
  "identically zero in every STy-conserving build" - that diagnosis does
  not survive: with MOCCa's linear antihermitian Rz, STy is the only
  antiunitary constraint on the stored basis and it fixes the MEs to be
  purely imaginary without annulling them. The imaginary half ships
  unvalidated until the step-6/7 odd-m oracle (deferral agreed
  2026-11-19). Gates: block_matrix_unit_tests 33506/0; NLO / NLO-T /
  NLO-P / NLO-TP mf + fam.NLO compile (debug profile, final source);
  frozen-make.inc NLO build; golden minimal NLO bit-identical.
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
  `exec/fam.NLO.reference.exe` is the parked FAM anchor; exec/MOCCa.NLO.exe
  fam.NLO.exe BXL/T/P/TP executables are freshly rebuilt and validated.
  `testing/integration/logs/fission.BXL.pre-move.out` is the inertia
  baseline produced by the pre-move BXL binary - batch the fission.sh
  comparison: rerun `bash fission.sh BXL` once the step's numerical gates
  are due, then `bash golden_diff.sh --no-run fission BXL BXL.pre-move`.
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
  0 of the design's §13. Design step 3 (§11/§13) migrates fam.f90's
  Qlm_spme consumer in its final commit and design step 9 targets FAM;
  FAM-side regressions are judged against the parked executable (no
  committed fam golden log exists, exec/ is untracked).
