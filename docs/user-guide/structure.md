# Input Namelists

MOCCa reads its input configuration from STDIN using Fortran namelist format. The namelists must appear **in a specific order**, and most are required for every calculation. Below is the complete list of namelists that the code accepts.

For more in-depth documentation of all options available in each namelist, see the corresponding files in [namelists](namelists/).

## Overview

The input file (typically named `mocca.in`) consists of a sequence of namelists, each starting with `&name` and ending with `/`. The **order of namelists is critical** — they must appear in the sequence shown below.

Most namelists are mandatory, but three are **conditionally optional**:

- `&indices` — only required when `BlockNumber > 0` in the `&pairing` namelist
- `&Inertia` — only required when `N_inertia > 0` in the `&IO` namelist  
- `&MomentConstraint` — only required when `MoreConstraints = .true.` in the preceding `&MomentParam` or `&MomentConstraint` namelist

Additionally, `&MomentConstraint` is the **only namelist that can be repeated**. It servers to specify multiple multipole constraints.

---

## Namelist Reference

| Namelist | Required | Description | Details |
|----------|----------|-------------|---------|
| [`&nucleus`](#nucleus) | Yes | Nucleus properties and convergence criteria | [&nucleus](namelists/nucleus.md) |
| [`&mesh`](#mesh) | Yes | Spatial mesh parameters (grid dimensions, spacing) | [&mesh](namelists/mesh.md) |
| [`&func`](#func) | Yes | Energy density functional options | [&func](namelists/func.md) |
| [`&pairing`](#pairing) | Yes | Pairing treatment options (HF, BCS, HFB) and blocking settings | [&pairing](namelists/pairing.md) |
| [`&indices`](#indices) | Optional[1] | Indices for blocking quasiparticles | [&indices](namelists/indices.md) |
| [`&evolution`](#evolution) | Yes | Parameters for the evolution of single-particle wavefunctions | [&evolution](namelists/evolution.md) |
| [`&scfiteration`](#scfiteration) | Yes | Options for self-consistent field iteration (mixing, preconditioning) | [&scfiteration](namelists/scfiteration.md) |
| [`&wfs`](#wfs) | Yes | Number of s.p. wavefunctions and transformation options | [&wfs](namelists/wfs.md) |
| [`&IO`](#io) | Yes | Input/output filenames and symmetry transformation options | [&IO](namelists/io.md) |
| [`&Inertia`](#inertia) | Optional[2] | Collective inertia calculation parameters | [&Inertia](namelists/inertia.md) |
| [`&MomentParam`](#momentparam) | Yes | General settings for multipole constraints | [&MomentParam](namelists/momentparam.md) |
| [`&MomentConstraint`](#momentconstraint) | Optional[3] | Specification of a single multipole constraint | [&MomentConstraint](namelists/momentconstraint.md) |
| [`&Cranking`](#cranking) | Yes | Cranking frequencies and angular momentum options | [&Cranking](namelists/cranking.md) |


1. Required only if `BlockNumber > 0` in `&pairing`
2. Required only if `N_inertia > 0` in `&IO`
3. Required only if `MoreConstraints = .true.` in `&MomentParam` or preceding `&MomentConstraint`; **can be repeated!**

---

## Basic descriptions

 We provide a list here of basic functionality in each namelist; this documentation file certainly does not exhaust all options. For full details, see the dedicated file for each namelist in [namelists](namelists/).


### &nucleus

See [&nucleus](namelists/nucleus.md) for a complete list of options.

**Purpose:** Specifies the nucleus to be studied and the convergence criteria.

**Basic parameters:**

- `neutrons`, `protons` — number of neutrons and protons
- `energy_prec`, `moment_prec`, `disp_prec`, `fermi_prec` — convergence thresholds for energy, multipole moments, dispersion, and Fermi energy

---

### &mesh

See [&mesh](namelists/mesh.md) for a complete list of options.

**Purpose:** Defines the spatial discretization grid for the calculation.

**Basic parameters:**

- `nx`, `ny`, `nz` — number of mesh points in each Cartesian direction
- `dx` — mesh spacing in femtometers (fm)

---

### &func

See [&func](namelists/func.md) for a complete list of options.

**Purpose:** Specifies which energy density functional parameterization to use.

**Basic parameter:**

- `name_param` — name of the parameterization, which must match a corresponding `.param` file in the current directory (e.g., `SLy4.param` or an entry in `forces.param`)

---

### &pairing

See [&pairing](namelists/pairing.md) for a complete list of options.

**Purpose:** Controls the treatment of pairing correlations and blocking in the calculation.

**Key parameters:**

- `Type` — pairing approximation: `HF` (Hartree-Fock, no pairing), `BCS`, or `HFB` (Hartree-Fock-Bogoliubov)
- `BlockType` — blocking scheme (0-4: no blocking, true blocking by indices/energy, EFA by indices/energy)
- `BlockNumber` — number of quasiparticles to block; if > 0, triggers reading of `&indices` namelist
- `pairingscheme` — solver type (0: direct, 1: gradient)
- `Bogofromfile` — whether to continue Bogoliubov transformation from input file

---

### &indices

See [&indices](namelists/indices.md) for a complete list of options.

**Purpose:** Specifies which single-particle states to block in pairing calculations.

**Required only if:** `BlockNumber > 0` in `&pairing`

**Basic parameters:**

- `BlockIndices` — list of state indices (Hartree-Fock basis) to block
- `BlockLowest` — alternative specification using string codes like `'n+'`, `'p-'`, etc. to block lowest states by quantum numbers

---

### &evolution

See [&evolution](namelists/evolution.md) for a complete list of options.

**Purpose:** Controls the evolution of single-particle wavefunctions during the self-consistent iteration.

**Basic parameters:**
- `maxiter` — maximum number of mean-field iterations
- `printiter` — interval for verbose output (collective corrections are only recalculated at these intervals)
- `Strategy` — evolution strategy: `HEAVYBALL` (default) or `IMTIME` (imaginary time)

---

### &scfiteration

See [&scfiteration](namelists/scfiteration.md) for a complete list of options.

**Purpose:** Controls the self-consistent field iteration procedure, particularly density and potential mixing/preconditioning.

**Basic parameters:**

- `scfscheme` — type of SCF evolution (0: potential preconditioning, 1: density mixing)
- `denmix` — density mixing parameter (default: 0.75)
- `preconfactor` — potential preconditioning strength

---

### &wfs

See [&wfs](namelists/wfs.md) for a complete list of options.

**Purpose:** Specifies the number of single-particle wavefunctions and options for their initialization.

**Basic parameters:**

- `nwn`, `nwp` — number of neutron and proton wavefunctions
- `osc_freq` — harmonic oscillator frequencies for Nilsson orbital generation (when starting from scratch)

**Important:** These numbers must match the input `.wf` file if restarting from a previous calculation (modulo symmetry transformation rules).

---

### &IO

See [&IO](namelists/io.md) for a complete list of options.

**Purpose:** Input and output file management, plus symmetry transformation options.

**Basic parameters:**

- `Inputfilename` — name of the `.wf` file to read (use `'INIT'` to start from scratch)
- `Outputfilename` — name for the output `.wf` file
- `BXLFIT` — prefix for output files used by Brussels fitting codes (e.g., `${BXLFIT}zXXXnYYY.out`)
- `denfile`, `potfile`, `sphffile`, `spcanfile`, `inertfile` — filenames for optional detailed output
- `N_inertia` — number of collective degrees of freedom for inertia calculation; if > 0, triggers reading of `&Inertia` namelist
- `Extraspwfs` — 8 integers for adding random wavefunctions in each isospin-parity-signature block
- `AllowTransform` — enables mesh, wavefunction, or symmetry transformations

---

### &Inertia

See [&Inertia](namelists/inertia.md) for a complete list of options.

**Purpose:** Specifies which collective inertias to calculate.

**Required only if:** `N_inertia > 0` in `&IO`

**Basic parameters:**

- `Inertia_l`, `Inertia_m` — lists of (l, m) values for which to calculate collective inertias

---

### &MomentParam

See [&MomentParam](namelists/momentparam.md) for a complete list of options.

**Purpose:** Global settings for multipole moment constraints and calculations.

**Basic parameters:**

- `MoreConstraints` — if `.true.`, signals that `&MomentConstraint` namelists follow
- `ContinueAll` — if `.true.`, use multipole constraint data from input `.wf` file

---

### &MomentConstraint

See [&MomentConstraint](namelists/momentconstraint.md) for a complete list of options.

**Purpose:** Specifies a single multipole moment constraint.

**Required only if:** `MoreConstraints = .true.` in `&MomentParam`

**Can be repeated:** Yes — use multiple instances to constrain multiple multipole moments.

**Basic parameters:**

- `l`, `m` — multipole indices
- `Constraint` — target value for the constrained multipole moment
- `Iteration` — number of iterations to apply constraint (-1 for entire run)
- `MoreConstraints` — `.true.` if more constraints follow, `.false.` for last one
- `Multfromfile` — if `.true.`, use Lagrange multiplier from input file

**Note:** For quadrupole constraints, an alternative parameterization using `iq1`, `iq2` is available.

---

### &Cranking

See [&Cranking](namelists/cranking.md) for a complete list of options.

**Purpose:** Controls cranking calculations for studying rotational states.

**Basic parameters:**

- `omegaX`, `omegaY`, `omegaZ` — cranking frequencies in MeV (hbar x omega)
- `crankX`, `crankY`, `crankZ` — target expectation values for angular momentum components
- `CrankTypeX/Y/Z` — type of cranking constraint (0: constant omega, 1: projection on feasible space)
- `continuecrank` — if `.true.`, use omega values from input `.wf` file

**Important:** Non-zero cranking frequencies require time-reversal breaking. Only `omegaZ` can be non-zero in standard code versions.
