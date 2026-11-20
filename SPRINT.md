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
| 1 | Core block_matrix types, classifiability init, fill machinery (§3-§5, §10 steps 1-2) | LANDED | 395bde9d, 55272151 | block_matrix unit net |
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
| 13 | D-6 dead-surface simplification: transposed flag / new_block_matrix / fill_elements dropped, Nabla_basis_pointers inlined (design §3/§5/§6 amended) | LANDED | 2ccb5f4c, 7771f3be | block unit net 50254 checks green; 4-config frozen -O3 mf compile; golden minimal 4849/4849 identical |
| 14 | COMMENT CLEANUP of the sprint's production modules (parked deferral discharged; the testing arms keep their verification-contract banners) | LANDED | f3fbe9f1, 38cb86f6, df48b5c4 | frozen -O3 NLO mf compile per module; block matrix unit net green |
| 15 | block_matrix developer reference: docs/user-guide/block-matrix.md + nav | LANDED | b9483982 | none (documentation; mkdocs build blocked locally by missing pymdownx) |
| 16 | D-7 row 1 (design §15.9): type layer - `$GOODQN` bit 4 in Hephaestos, tau/folded components + TAU_EVEN/TAU_ODD, mandatory tau at every init call site, the bit-4-gated honest-trace contract full_trace (+ the folded-composite refusal) | LANDED | ccf485d4 | block unit net 50260 checks 0 failures, refusal scenarios 1-9 all refused; 4-config frozen -O3 mf compile 0 errors, per-config masks as recorded in design §2c; golden minimal NLO 4849/4849 identical |
| 17 | D-7 row 2 (design §15.9): interface collapse - the get_J_block accessor in the get_nabla_block shape (F-12 fix: effective-word gate in both accessors), Jt_spme_matrix deleted, NablaT_spme_matrix private behind Nabla_half_fm, arms build their own folded twins and assert the accessor packaging + the production lane | LANDED | 67a0d31d | arms green NLO/NLO-T/NLO-P with moments net + cross-build 0 (fold agreement 0.0 on NLO and NLO-P, tau both-ways 4.4e-16/1.1e-16, halffm lane 0.0); NLO-TP mf compile; NLO mf+fam compile 0 errors; golden minimal NLO 4849/4849 identical; block unit net 50260/0 |
| 18 | D-7 row 3 (design 15.9): debug windows - get_J_complete/get_nabla_complete, the eager doubled-basis complete-matrix materialisation per family (bit-4-cleared init, honest kernels on materialised T-images); arms assert window/accessor agreement per block, the full_trace contract and the bit-4-clear degeneration; the gradient fullbox_nabla_me + check_fullbox discharges the 15.6 opposite-half adjudication | LANDED | 7179a85c | arms green NLO/NLO-T/NLO-P with moments net + cross-build 0 (window vs accessor 0.0-4.4e-16 with image-image blocks at FULL tau coverage, degeneration exact, trace contract holds; 15.6 verdict: residues, no kernel bug, numbers in design 15.6); NLO-TP mf compile; NLO mf+fam compile; golden minimal NLO 4849/4849 identical; block unit net 50260/0 |
| 19 | NAME UNIFICATION (user-directed cleanup, not a design step): sym_matrix type + *_sym_matrix procedures -> block_matrix everywhere; docs/user-guide/sym-matrix.md -> block-matrix.md; zero remaining sym_matrix occurrences | LANDED | 963e6d3a | block unit net + refusal scenarios green; moments net + cross-build green NLO/NLO-T/NLO-P; NLO-TP + NLO fam compile; golden minimal NLO 4849/4849 PASS (0.17% tolerance) |
| 20 | Consumer migration (design 14.4 row 3, D-8): calcJ2andBelyaev_HFB fills jx/jy/jz from the scenario-blind Jmu_spme - honest trio in T-broken builds, honest J_z plus the folded pair in T-conserving ones, the choice on goodqn bit 4 inside the producer; j-matrix fill source marker-free, kernel arithmetic and contraction gates unchanged | LANDED | e3501416, 7aab4b63 (rename) | symmetries.sh SLy4 (four builds, cross-build consistency 0) with pre-migration worktree A/B: all four logs bit-identical in physical content, timers and 1e-15 commutator residues excepted; golden minimal NLO 4849/4849 identical; block unit net 50260/0; NLO fam + 4-config frozen mf compile |

Test-infrastructure commits alongside: cd8d9b8e (moments_testing + unit_test
flag to &mesh/), 7dcc6118 + 57217a28 (cross-build runscript + CI job),
b64a9707 (block_matrix CI job), ed8cac5f (fam golden reference for the
spurious-mode scenario).

## NEXT UP

1. The sym_rotation consumers of §6a: the cantransfo legs of
   calcJ2andBelyaev_HFB through sym_rotation (design 6a), with
   fission.sh as the heavy gate at the end of the chain; also D-8
   follow-ups (the folded flag of Nabla_half_fm and the remaining
   `$TR`/`$NTR` contraction gates in the migrated routine are D-5-adjacent).
2. Then the §13.5 deferred items (odd-m / Im-operator FAM fields end
   to end,
   Im-operator Ksum/calc_Q20 generalisation), or first dispose of the
   remaining Qlm_spme reference (DECISIONS.md D-2).

Deferred / parked (not lost - do not re-litigate):
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
  bounds checks, no HDF5). The debug profile is again usable for every
  CONFIG now that the Coulomb crash of FOUND_BUGS.md F-1 is fixed.
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
