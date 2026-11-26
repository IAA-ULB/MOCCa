# QOL backlog — small issues for a future quality-of-life sprint

Small papercuts noticed while working on other things. None block current
sprints; collect and fix in one batch later. Delete entries as they are fixed.

1. `writewavefunction` (src_orig/IO.f90, ~line 749) slices
   `OutputFileName(len_trim-3:)`: any `OutputFilename` shorter than 4
   characters crashes with a runtime substring error in bounds-checking
   builds (observed with `OutputFilename='mf'`, NLO gnu-debug). All
   integration-style runs with short names dodge it by luck.

2. The extended-input printout in `TransformInput` (src_orig/transform.f90)
   has a format typo: it prints the header `nwn, nwn` where the second field
   should be `nwp`.

3. `iterate_spectrum` (src_orig/mocca.f90) carries a batch of unused-import
   warnings (`nwt`, `hfblocks`, `hfblocks_global`, `hfpsi`, ... in the
   wavefunctions use-lists) — leftovers from code that moved into
   `add_spwfs`. Harmless, but noisy in debug builds.

4. Locked spwfs are skipped in the evolution and keep their last computed
   dispersion value; every debug print that reports
   `maxval(dispersions)` (e.g. `print_summary_spectrum`) therefore reports a
   stale value once anything locks. If this is ever the basis for a
   convergence decision (instead of the per-block locking chain), clear the
   dispersions of locked spwfs when they lock — or annotate the print.

5. `TransformInput` zeroes `rho_can` wholesale when extending the spwf set
   (input path). Anyone adding spwfs at input must know to call SolvePairing
   afterwards; a comment or a refill inside the routine would help.
