# SPRINT.md — matrix_refactor sprint: status snapshot

Status of the symmetry-block-matrix refactor sprint. One row per landed step;
detailed narrative lives in the git commit messages and in the design document's
amendment sections - never restated here.

- Authoritative design: `mocca-symmetry-block-matrix-design-summary.md`
  (repo root). Read it before touching block_matrix_module or any generator.
- Sprint constraints and operational workflow: sections below.
- Open user decisions: `DECISIONS.md`. Out-of-scope findings: `FOUND_BUGS.md`.
- Machine-local state (executables, anchors, scratch dirs):
  `NOTES.local.md` (gitignored).

Maintenance rule ( AGENTS.md "Documentation protocol" is binding ):
this file is REWRITTEN when a step lands - flip the status column, move the
NEXT UP pointer - never appended to as a diary.

## Step table

| # | Step (design ref) | Status | Commit(s) | Gate(s) executed |
|---|---|---|---|---|
| 0 | Park FAM anchor `fam.NLO.reference.exe` (§13.3 row 0) | LANDED | tree at b719453e | none needed |
| 1 | Core sym_matrix types, classifiability init, fill machinery (§3-§5, §10 steps 1-2) | LANDED | 395bde9d, 55272151 | block_matrix unit net |
| 2 | goodqn mask + injected `$GOODQN` (§2c, §10 step 3 / §13.3 row 1) | LANDED | 6a6d92c8 | masked-init unit tests; 4 configs compile |
| 3 | Multipole generator Qlm_spme_matrix/fm, both halves, chamber-weighted fill, spme_oracle hook, inertia switch, units flip (§13.3 rows 2-5) | LANDED | 9faf8238, 810986ab, d7f42ff2, 59ae4252, c0262789, 0b6de03f | golden minimal NLO 4849/4849 bit-identical; oracle agreement ~1e-15 (NLO/NLO-T/NLO-P); O16 inertia bit-identical vs legacy; fission.sh via CI vs pre-move baseline |
| 4 | Full-box reconstruction oracle, both halves (§10 step 4, amended) | LANDED | 497b1d94 | test_Qlm_fullbox green on NLO/NLO-T/NLO-P (identity ~1e-13 reassociation level, opposite halves <= 6.1e-16) |
| 5 | FAM get_qpme_op fed from Qlm_spme_fm ('multipole' and 'zcom') (§10 step 5) | LANDED | 20fca2e5 | fam unit tests SLy4; A/B vs parked fam anchor: Q20/Q22 famfiles bit-identical, GMRES residuals identical |
| 6 | Legacy Qlm_spme call-site cleanup (verify_COM_motion, test_RP_commutator) (§13.2) | LANDED | 4ca0f181 | golden minimal identical; fam + moments nets green; compiles |
| 7 | Basis rotations: sym_rotation + rotate (§6a, §10 step 7) | LANDED | aa95fb26 | unit net 53562 checks + 8 refusal scenarios; golden minimal identical; moments suite green incl. cross-build; all 4 configs compile |
| 8 | Angular-momentum generators + kernel-transpose oracle arm (§14.4 row 2; row 1 K_x fix CANCELLED by measurement §14.2) | LANDED | 51fa4e00 | arm green NLO (debug), NLO-T (debug), NLO-P (frozen -O3); symmetries.sh A/B for the K_x audit |
| 9 | 14.2b adjudication: dead-kernel pairing fix + check_fullbox full-box oracle (§14.2b, §14.3) | LANDED | 2421931f, 0792d678 | opposite halves <= 2.1e-17 all builds; "jz is real" premise judged green; kernels fixed stay caller-less |
| 10 | Gradient generators (honest nabla trio + folded twins, get_nabla_block, tau contract) + arm (§15.5 row 1) | LANDED | de496fd0 | arm green NLO/NLO-T (debug, legacy agreement 2-5e-15), NLO-P (frozen -O3); NLO-TP compile |
| 11 | CompCOMCorrection fed from linear_momentum generators (§15.5 row 2) | LANDED | 0a36d261 | row-0 scenarios vs parked pre-migration binaries: energies bit-identical ($TR), log body bit-identical ($NTR cranked Mg24); golden minimal identical |
| 12 | verify_COM_motion migrated to Nabla_half_fm (§15.5 row 3) | LANDED | b68ba67b | source-level; NLO/NLO-T mf compiles |

Test-infrastructure commits alongside: cd8d9b8e (moments_testing + unit_test
flag to &mesh/), 7dcc6118 + 57217a28 (cross-build runscript + CI job),
b64a9707 (block_matrix CI job), ed8cac5f (fam golden reference for the
spurious-mode scenario).

## NEXT UP

1. calcJ2andBelyaev_HFB consumer migration (design §14.4 row 3) - feeds the
   sym_rotation consumers of §6a.
2. Then the §13.5 deferred items (odd-m / Im-operator FAM fields end to end,
   Im-operator Ksum/calc_Q20 generalisation), or first dispose of the
   remaining Qlm_spme reference (DECISIONS.md D-2).

Deferred / parked (not lost - do not re-litigate):
- COMMENT CLEANUP of the sprint's marker-free modules (raised 2026-11-19
  after the AGENTS.md style rules landed): block_matrix_module.f90,
  angular_momentum.f90, linear_momentum.f90 and their *_testing.f90
  arms measure 28-54% comment lines - they inherited the verbose banner
  scale the old AGENTS.md example taught. Slim them to the new economy
  (banners only for signs, units, marker keys, quantum-number
  encodings; no dates, commit hashes, design restatements). Scope:
  these ~7 files ONLY - no mass rewrite of legacy sources (minimal
  diff wins there). Pure-comment change, no numerics touched, so the
  gate is a frozen-build compile per touched config plus the block
  matrix unit net; do it in one dedicated commit per module, never
  mixed with a physics change.
- §15.6 derivative opposite-half adjudication: the linear-momentum family
  still needs its own full-box reconstruction oracle; the angmom arm's
  check_fullbox (step 9) is the ready template. Artifact reading favoured by
  mechanism; production unaffected.
- §15.4 DEN-vs-cantransfo audit; sym_rotation production consumers
  (candidates in design §6a).
- FAM rows beyond step 5 (CompNablaMelements deletion, test_RP_commutator)
  - deferred by user direction; see DECISIONS.md D-3/D-4.
- Spin operators beyond J (the spin_* kernels) - no design section yet.

## Conventions (user-confirmed; full reasoning in commit messages and design)

- Folded meshes: TR/conservative builds carry the z-fold in dv, so nz=12 is a
  24 fm z box; parity-broken builds store z UNFOLDED and need nz=24 for the
  same box (nz=12 converges to a different truncated-box state, ~88 keV off
  for O16/SLy4 - broken input, not new physics). Broken-symmetry inputs also
  demand even nwn/nwp and allowtransform=.true. in &IO/.
- MOCCa's time reversal is FERMIONIC: T**2 = -1 (design §14). The antilinear
  composite J_mu.T is HERMITIAN as operator; both folded twins MEASURE kind +1
  symmetric storage; K_y = J_y.T carries herm = ANTIHERMITIAN purely as the
  kind-completing label (herm = kind * sty). All kind signs are test
  assertions settled by oracle, never re-derived (design §14.1/§14.2).
- Odd-m multipole operators are structurally zero in TR-conserving builds
  (stored sectors all q(Rz)=+i, XOR partners empty). Their wobbling-relevant
  content exists only when TR-broken blocks are populated.
- One-body Q^2 through sp matrices presumes the truncated basis RESOLVES Q|i>
  - false for any MOCCa basis (failure ~1e4 relative, worst at high l). No
  Qlm^2-vs-mesh Squared test exists or can exist.
- Cross-build comparisons are meaningful only through basis-independent
  invariants (traces) of common operators (even m) on the SAME state: the
  moments suite compares the three builds on their common initialisation with
  maxiter=0 (motivated by a CI failure where each build stopped at a different
  point of a near-flat convergence landscape - environment-sensitive by
  construction, not a machinery error).
- Parity-odd FAM FIELDS in parity-conserving builds NaN by construction (the
  FAM machinery cannot represent a parity-breaking excitation; identical in
  the parked pre-refactor anchor). The working scenario is the parity-broken
  build. NLO-P fam validation scenarios: see FOUND_BUGS.md F-7 and NOTES.local.
- The unit_test flag lives in &mesh/ (geninfo) and arms every executable's
  own suite: moments from Run_MOCCa, fam from fam_run, plus the
  angular-momentum and linear-momentum arms (wired in mocca.f90, before
  run_moments_tests).
- Storage single-node by design (design §1/§12); keep block_matrix_module
  free of rank-dependent logic and of `$`-markers; the sprint is serial-only.

## Operational workflow

- Quick compile iterations: debug profile -
  `make CONFIG=X clean` then `make CONFIG=X mf INCLUDE=gnu-debug` (-O0,
  bounds checks, no HDF5). Caveat: debug-profile NLO-P builds crash on a
  pre-existing Coulomb bug (FOUND_BUGS.md F-1) - use the frozen -O3 build
  for NLO-P gates.
- The root `make.inc` stays frozen (-O3 + HDF5) and must back every numerical
  gate or golden comparison (it is the trust anchor for the golden logs and
  parked executables). make.inc is outside the dependency graph: an
  INCLUDE= switch without clean yields mixed binaries silently.
- Ancestors/anchors live in exec/ (gitignored, machine-local): see
  NOTES.local.md for what is validated and what is stale.
- Oracle-run recipe (spme_oracle hook in fission_MOI): copy the
  parameterization (e.g. SLy4.param) next to the executable; even nwn/nwp +
  allowtransform=.true. for broken-symmetry builds; arm with N_inertia>=5 in
  &IO/, an &inertia/ namelist (first four moments hardcoded Q20/Q22/Q30/Q32)
  and spme_oracle=.true. Production-like NLO-T scenario: Mg24 HFB
  pairingscheme=1, omegaZ=0.3292333881 (cranking.sh style; pure HF inputs
  force diagsphamil on, which reshuffles the stored basis within blocks
  mid-run - judge oracle numbers on production-like inputs only). The hook
  fires per calc_collective_inertia call - keep maxiter small in dev runs.
- fam_unit_tests.sh needs --parameterisation SLy4 under NLO fam builds - the
  t0t3 default stops silently without running anything (FOUND_BUGS.md F-6).
