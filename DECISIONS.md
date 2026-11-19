# DECISIONS.md — decision queue for Wouter

One numbered entry per decision. Entries are closed by replacing status with
DECIDED, the date and one line of resolution - the entry then stays as the
record. Do not delete entries; do not open duplicates of closed ones.
Background lives in the design summary sections referenced - read those
before deciding.

## Open - user to decide

### D-1 Folded-composite class labelling (JxT/JyT) — OPEN (raised 2026-11-19)

The measured folded-kinds resolution (design §14.1/§14.2/§2d) is operational,
not principled: both folded kernels measure symmetric kind +1 storage, and the
class bookkeeping that expresses this inside the type is a convention pinned
by measurement - J_x.T is tagged HERMITIAN x STY_EVEN, while J_y.T combines
symmetric kind with STY_ODD and therefore carries herm = ANTIHERMITIAN purely
as the kind-completing label (herm = kind * sty), not as a physical adjoint
class (the §2a kind machinery does not transfer to antilinear composites).

Question: is overloading `herm` as a storage label acceptable as the permanent
resolution, or should the type grow an explicit folded-composite constructor
path that names the antilinear composites honestly?

Constraint whichever way: the kernel-transpose oracle (and for derivatives
the mirror-kind oracle) stays the arbiter of every kind sign - never re-derive,
re-measure.

### D-2 Qlm_spme deletion gate — OPEN (design §13.2)

Purely a decision now. The routine's only remaining reference is the
off-by-default spme_oracle hook in fission_MOI.f90. Either Qlm_spme and the
oracle go together, or both stay permanently as the production-scale
cross-check of Qlm_spme_matrix. The moments unit net covers the generator
either way.

### D-3 CompNablaMelements deletion (design §15.5 row 4) — OPEN, deferred to the FAM phase

The remaining references are exactly: fam.f90 zmomentum, fam_testing.f90
test_RP_commutator (both user-deferred FAM phase), the definition in
densities.f90, and linear_momentum_testing's A/B oracle. Take the deletion
decision together with the FAM work.

### D-4 Dead code: verify_COM_motion and test_RP_commutator — OPEN

Both remain call-less (verify_COM_motion is a pushing-model diagnostic with
no call site anywhere; test_RP_commutator is commented out in fam_run.f90).
Enabling either - e.g. RP as an NLO-P zcom net - is a user decision, out of
the migration's scope. Both have been migrated to the generators so they are
build-clean if ever enabled.

## Open - technical adjudications (deciding test named; no user input needed)

### D-5 fac/4 and 2x TR factors in CompCOMCorrection (design §15.6) — UNADJUDICATED

The factors are empirically certified as preserved (A/B bit-identical,
§15.5 row 2) but whether they are mechanistically the §15.8 mechanism-1
T-pair mirroring is not settled. Settle if the mechanism is ever needed for a
future consumer, not before.

## Decided / resolved (record)

- 2026-11-19: matvec and matmat are DROPPED from sym_matrix and the plan -
  the consumer audit found no caller (design §6/§10, commit 1230f34a).
- 2026-11-19: basis rotations are sector-diagonal real orthogonal ONLY
  (complex unitaries fall out of the STy-invariant model; signature-mixing
  rotations need a basis no config produces) - design §6a.
- 2026-11-19: the suspected K_x fill bug is REFUTED by A/B measurement -
  the symmetric `$TR` fill is correct; design §14.4 row 1 cancelled,
  §14.2 carries the audit.
- 2026-11-19: the angmom opposite-half O(1) purity prints were
  reduced-mesh chamber residues, not kernel errors (design §14.2b RESOLVED);
  the "jz_ij is always a real number" premise is judged green. One genuine
  transcription bug existed in a dead kernel and was fixed (caller-less,
  FOUND_BUGS.md F-2).
- 2026-11-19: parity-odd FAM fields in parity-conserving builds NaN by
  construction, identically in the parked pre-refactor anchor - not a
  regression; making them available is a long-term goal of this refactor.
- 2026-11-19: parity-broken merged blocks: only blocks 1 and 5 are populated
  at runtime, every stored state q(Rz) = +i - the s = block - 1 encoding
  stands in every build (design §2b/§2d).
