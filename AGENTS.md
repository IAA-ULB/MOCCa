# AGENTS.md — MOCCa

Working notes for AI coding agents in this repository. Everything below was
verified against the code and tooling on 2026-11-16, on branch `matrix_refactor`
at commit 9960e8b2; the documentation was restructured on 2026-11-19 (see the
"Sprint documentation map"). If the file disagrees with the code, the code wins —
update this file.

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
`.mod`) and it is regenerated on the next build. This includes
`make CONFIG=<X> clean`: it is GLOBAL - it wipes every config's
`build/` tree, not just X (executables in `exec/` survive) - so never
run it in a loop over configs. The per-config incremental
`make CONFIG=<X> mf|fam` is safe after changes to module internals:
only the touched files and their dependents recompile. Nothing in
`exec/` is tracked by git (see `.gitignore`), nor is the local
`make.inc`.

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
- Adding a new `src_orig/*.f90` file: it must be listed in `MF_SRC`/`FAM_SRC`
  in the Makefile and in `no_process_files` in `src_heph/preprocess.py`
  (plus `FORTRANFILES` in `Hephaestos.py` for the dry-run report), or the
  per-file generation silently skips it (empty output file, compile
  failure).

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
- In-code documentation is bounded by information content, not routine
  size. The house banner format (`!----` fences, Input/Output sections)
  is reserved for routines whose semantics are genuinely non-obvious -
  symmetry sign conventions, units, quantum-number encodings, `$`-marker
  keys, operator classes. Everything else gets one to three `!` lines of
  purpose, plus Input/Output lines only where the argument names do not
  carry the meaning.
- Comment the why and the conventions, never the what. No revision
  dates, commit hashes, or design-summary restatements inside code
  comments - the commit message and the design summary are those
  records. Do not imitate a verbose legacy banner when editing near it,
  and do not rewrite existing banners without other cause (minimal diff
  wins); the conciseness rules bind NEW and MOVED code.
- Self-check before committing: if a comment line can be deleted without
  a routine's user losing information they need, delete that line.
- An example of the target economy (deliberately the same routine that
  used to appear here as a 25-line banner):
  `subroutine set_spwf_symmetries(sx,sy,sz,blocks)
    !---------------------------------------------------------------------------
    ! Assigns the x/y/z reflection signs of all spwf spinor components;
    ! Hephaestos fills them via the $SXab marker key:
    !   X = direction, a = symmetry-block index, b = spinor component.
    !---------------------------------------------------------------------------
    ! Input:  blocks - spwf counts per symmetry block (may be a subset
    !                  when reading from file)
    ! Output: sx, sy, sz - per-component reflection signs
    !---------------------------------------------------------------------------`
- Commits: subject line ≤ 60 chars, verb first, imperative mood.
  Body ≤ 5 lines, only the *why* — never restate the diff or the
  test output. No environment narration ("verified locally", etc.).

## Language 
- Be concise yet precise, in prose and especially in code comments; shorter
  than the Style example above means something was lost, longer means
  something was restated.
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
  no MOCCa executable needed) including the refusal scenarios of the
  rotation structure guards (8) plus the full_trace folded-composite
  refusal; its CI job `block_matrix_unit_tests` needs
  nothing but lint-and-format, so it gives the fastest feedback of the
  whole workflow.
- Sprint unit arms (moments suite, the angular-momentum and
  linear-momentum arms of the matrix_refactor refactor): armed by the
  `unit_test` flag in `&mesh/` and run per executable from `Run_MOCCa` /
  `mocca.f90`. Sprint status, scenario recipes and caveats live in
  `SPRINT.md`, known traps in `FOUND_BUGS.md`.
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

## Agent tooling gotchas (session-proven, 2026-11-20)

- Editing tracked docs with string-replacement tools: old_strings
  containing the section sign or Fortran backslash sequences
  (`\nabla`) fail to match on a byte-encoding mismatch, not a bad
  anchor. Anchor the edit on a nearby line without those characters,
  or use `sed -i`; keep the same characters out of the sed pattern
  too by matching surrounding words. After any line-number surgery on
  the sprint documents, re-read the region: one splice has mangled
  neighbouring lines before.
- Long builds and test runs: start them as background processes,
  single-quoted `bash -lc '...'`, no shell variables (an unquoted `$C`
  config loop once silently no-op'd), logging to `/tmp` files; the
  resulting exes and the gate logs outlive the session via
  `NOTES.local.md`.

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

49 Fortran modules over 51 source files in `src_orig/`. Regenerate with:

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
- Symmetry-block matrices (the matrix_refactor sprint: marker-free,
  plain Fortran 2003+, no module dependencies): block_matrix_module <-
  block_matrix_module.f90, angular_momentum + angular_momentum_testing,
  linear_momentum + linear_momentum_testing
- Densities, functional, pairing: densities, functional, parameterization,
  pairing, pairingcutoffs, pairing_strengths, Coulombmod <- coulomb.f90,
  cranking, moments <- moments.f90 (incl. moments_testing), momentsofinertia,
  fission_MOI
- I/O and infrastructure: IO, IO_aux, IO_wf, HDF5_auxiliary, Printing <-
  printing.f90, GenInfo <- geninfo.f90, compilation, timing, version

The ME / expectation-value territory of the code: `wavefunctions.f90`
(spwf storage, symmetry-block bookkeeping, the angmom_* / spin_* ME
routines), `momentsofinertia.f90` (J_mu matrices and inertia tensors),
`moments.f90` (multipole expectation values), `functional.f90`
(incl. CompCOMCorrection), `cranking.f90`, `basis_transform.f90`
(HF <-> canonical, HFBlocks usage), and `fam.f90` (drho, dkappa, FAM
perturbations).

## Sprint documentation map

The refactor sprint on this branch keeps its state in four companion
documents at the repo root. AGENTS.md carries evergreen repo facts only —
no step-by-step, no per-commit state:

| Document | Lifetime | Contents |
|---|---|---|
| `mocca-symmetry-block-matrix-design-summary.md` | permanent | The authoritative design: intent, conventions, amendments, measured records. Read first when touching `block_matrix_module` or any generator/consumer of `sym_matrix` |
| `SPRINT.md` | rewritten per step | Live status: one row per landed step (commits + gates), the NEXT UP pointer, sprint conventions, operational workflow |
| `DECISIONS.md` | until decided | Wouter's decision queue: open user decisions, pending adjudications, and the resolved record |
| `FOUND_BUGS.md` | until fixed | Issues and inconveniences found but out of sprint scope, each with reproduction and workaround |
| `NOTES.local.md` (gitignored) | per machine | Parked executables, validation status of local binaries, scratch locations. Never quote as repo fact |

## Documentation protocol (binding for agents)

1. Before any sprint step, read `SPRINT.md` (status, NEXT UP, conventions),
   `DECISIONS.md` (the queue — do not re-litigate closed entries) and the
   design sections your step touches. This file is not a source of sprint
   state.
2. Git is the record of what happened; the docs record why and what is
   decided. A landed step gets a commit plus ONE `SPRINT.md` table row. 
   Never restate the narrative in a second document.
3. `SPRINT.md` is rewritten, never appended to as a diary: flip the status
   column when a step lands, move the NEXT UP pointer, and leave the rest
   alone.
4. Measured conventions and amendments go into the design summary at the
   sections they amend, dated and commit-referenced. Every measured number
   lives in exactly one place; the SPRINT row points at it.
5. New out-of-scope findings get a `FOUND_BUGS.md` entry when found; new
   user decisions get a `DECISIONS.md` entry. Close entries by flipping
   their status — never delete them.
6. Machine-local facts (executables, which binaries are validated) belong in
   `NOTES.local.md`, never in tracked documents.
7. If a document disagrees with the code, the code wins — fix the document
   in the same commit that changed the code.
