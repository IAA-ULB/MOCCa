# MOCCa Test Suite

This document provides an overview of all operational tests available in the MOCCa testing framework, located in the `testing/` directory. The tests are organized into three main categories: **unit tests**, **integration tests**, and **regression tests**.

All tests use the common boilerplate functions defined in `testing/functions.sh` for setup, teardown, and result extraction.

---

## Test Infrastructure

### Common Functions (`testing/functions.sh`)

The testing framework provides the following helper functions for test scripts:

| Function | Purpose |
|----------|---------|
| `setup_test_env` | Sets up working directory, copies executables and parameter files for mean-field calculations |
| `setup_test_env_fam` | Sets up for FAM (Finite Amplitude Method) calculations |
| `setup_test_env_dep` | Sets up for tests with external potential files |
| `teardown_test_env` | Cleans up working directory after test completion |
| `get_total_energy_stdout` | Extracts final total energy from MOCCa output |
| `get_coulomb_energy_stdout` | Extracts Coulomb energy from MOCCa output |
| `get_Z_stdout` | Extracts proton number from MOCCa output |
| `get_rms_stdout` | Extracts RMS radius from MOCCa output |
| `get_B20_stdout`, `get_B22_stdout` | Extracts quadrupole deformation parameters |
| `get_muz_stdout` | Extracts z-component of magnetic moment |
| `get_strength` | Extracts strength from FAM output at given energy |
| `get_Belyaev_stdout` | Extracts Belyaev moment of inertia |
| `get_DJ2_stdout` | Extracts dispersion of J^2 |
| `get_inertia_components` | Extracts inertia tensor components (2020, 2030, 3030) |
| `get_neck_stdout` | Extracts neck value from MOCCa output |
| `extract_strength_components` | Extracts strength components from FAM file |
| `check_convergence` | Checks if MOCCa calculation converged |
| `compare_floats` | Compares two floating-point numbers within tolerance |

---

## Compilation Tests

The CI/CD workflow (`.github/workflows/workflow.yml`) systematically tests compilation across multiple configurations, compilers, and build types to ensure MOCCa can be built reliably in various environments.

### Compilation Matrix

#### Mean-Field Executables (NUCLEI)

| Category | Options | Notes |
|----------|---------|-------|
| **Configuration** | LO, LO-T, NLO, NLO-T, BXL, BXL-T, BXL-P, BXL-TP, N2LO, N2LO-T, N2LO-P, N2LO-TP, BXL-N2LO | 13 different symmetry configurations |
| **Compiler** | gfortran, ifx | GNU Fortran and Intel Fortran Compiler |
| **MPI** | SERIAL, MPI | Serial and parallel builds |
| **HDF5** | Enabled for gfortran, disabled for ifx | HDF5 support for data storage |


#### Mean-Field Executables (PASTA)

| Category | Options | Notes |
|----------|---------|-------|
| **Configuration** | BXL | PASTA-specific configuration |
| **Compiler** | gfortran, ifx | GNU Fortran and Intel Fortran Compiler |
| **MPI** | SERIAL, MPI | Serial and parallel builds |
| **HDF5** | Enabled for gfortran | HDF5 support for data storage |


#### FAM Executables (NUCLEI)

| Category | Options | Notes |
|----------|---------|-------|
| **Configuration** | LO, LO-T, LO-P, NLO, NLO-T, NLO-P, BXL | FAM-specific configurations |
| **Compiler** | gfortran, ifx | GNU Fortran and Intel Fortran Compiler |
| **MPI** | SERIAL | Currently only serial FAM builds |
| **HDF5** | Enabled for gfortran | HDF5 support for data storage |


### Special Compilation Cases

In addition to the matrix builds, the workflow includes specific targeted builds:

- **ifx with MPI**: NLO and BXL configurations are explicitly tested with Intel compiler and MPI enabled
- **Compiler independence**: Key configurations (LO, NLO, BXL) are built with both gfortran and ifx to verify compiler compatibility
- **HDF5 verification**: gfortran builds include HDF5 support by default

### Build Process

Each compilation job:

1. Checks out the code
2. Copies source files from `src_orig/`
3. Selects the appropriate `make.inc` file based on compiler family and MPI settings
4. Compiles with the specified CONFIG and CALCTYPE
5. Uploads the resulting executable as an artifact for subsequent test jobs

The executables are named following the pattern:
```
MOCCa.{CONFIG}.{COMPILER}.{CALCTYPE}.{SERIAL|MPI}.exe
fam.{CONFIG}.{COMPILER}.{CALCTYPE}.{SERIAL|MPI}.exe
```

For example:

- `MOCCa.NLO.gfortran.NUCLEI.SERIAL.exe`
- `MOCCa.BXL-T.ifx.NUCLEI.MPI.exe`
- `fam.LO.gfortran.NUCLEI.SERIAL.exe`

---

## Unit Tests

Unit tests verify the correctness of specific code components. 

### FAM Unit Tests

| Test Script | Description | Default Parameters | Verification |
|-------------|-------------|-------------------|--------------|
| [`fam_unit_tests.sh`](../testing/unit/fam_unit_tests.sh) | Runs FAM solver unit tests for a 18O nucleus in a minimal box | Type: HFB, Parameterization: t0t3, nw: 15 | Tests defined in Fortran code; checks exit status |

**Usage:**
```bash
bash fam_unit_tests.sh [EXESUFFIX] [--pairing HF|HFB] [--parameterisation PARAM] [--nw NW]
```

---

## Integration Tests

Integration tests verify that different components of MOCCa work together correctly and produce physically meaningful results.

### Mean-Field Calculations

| Test Script | Description | Nucleus | Verification |
|-------------|-------------|--------|--------------|
| [`minimal.sh`](../testing/integration/minimal.sh) | Spherical calculation of 16O with SLy4 in minimal box | 16O | Total energy: -128.513125 MeV (±1 keV), β20 = 0.0 (±0.0001), β22 = 0.0 (±0.0001) |
| [`fission.sh`](../testing/integration/fission.sh) | Tests fission barrier calculations | 240Pu | Convergence and energy criteria |
| [`lattice_energy.sh`](../testing/integration/lattice_energy.sh) | Tests lattice energy calculations | - | Energy comparison |
| [`spag.sh`](../testing/integration/spag.sh) | Tests single-particle spectrum and gap properties | - | Energy and gap verification |

### Symmetry and Invariance Tests

| Test Script | Description | Verification |
|-------------|-------------|--------------|
| [`symmetries.sh`](../testing/integration/symmetries.sh) | Tests equivalence of different symmetry modes (T, P, TP) for 24Mg constrained to triaxial shape | Total energy, Belyaev moment of inertia, dispersion of J2 agree across symmetry modes (±1 keV, ±1e-3 ħ2/MeV) |
| [`fam_rotational_invariance.sh`](../testing/integration/fam_rotational_invariance.sh) | Tests rotational invariance in FAM calculations | Strength values match expected symmetry properties |
| [`fam_sc_symmetry.sh`](../testing/integration/fam_sc_symmetry.sh) | Tests sign-conjugation symmetry in FAM | Strength values maintain SC symmetry |
| [`fam_strength_symmetry.sh`](../testing/integration/fam_strength_symmetry.sh) | Tests strength function symmetry properties | Strength values respect symmetry constraints |
| [`fam_pairing_zero_mode.sh`](../testing/integration/fam_pairing_zero_mode.sh) | Tests pairing zero mode in FAM calculations | Zero-energy mode properties |
| [`fam_linearity.sh`](../testing/integration/fam_linearity.sh) | Tests linear response in FAM | Linear response properties |
| [`fam_t0t3.sh`](../testing/integration/fam_t0t3.sh) | Tests FAM with t0-t3 parameterization | FAM results validation |

### Continuation and Restart Tests

| Test Script | Description | Verification |
|-------------|-------------|--------------|
| [`continuation_hdf5.sh`](../testing/integration/continuation_hdf5.sh) | Tests continuation using HDF5 file format | Successful restart and convergence |
| [`continuation_wf.sh`](../testing/integration/continuation_wf.sh) | Tests continuation using wavefunction files | Successful restart and convergence |

### Advanced Physics Tests

| Test Script | Description | Nucleus | Verification |
|-------------|-------------|--------|--------------|
| [`blocking.sh`](../testing/integration/blocking.sh) | Tests 4-step blocking workflow: false vacuum → EFA → EFA without T → full blocking | 55Cr | Energy and deformation comparison at each step |
| [`magmoment_Sc41.sh`](../testing/integration/magmoment_Sc41.sh) | Tests magnetic dipole moment calculation for Sc41 | 41Sc | μ_z = 5.80 (±0.01) |
| [`cranking.sh`](../testing/integration/cranking.sh) | Tests cranking calculations for rotational states | - | Angular momentum and energy verification |
| [`hfb_direct_gradient.sh`](../testing/integration/hfb_direct_gradient.sh) | Tests HFB with direct and gradient pairing schemes | - | Successful completion with both schemes |
| [`boundary_conditions.sh`](../testing/integration/boundary_conditions.sh) | Tests different boundary condition treatments | - | Convergence and energy stability |
| [`store_derivatives.sh`](../testing/integration/store_derivatives.sh) | Tests storage and retrieval of derivative quantities | - | Data consistency checks |

### Memory and Performance Tests

| Test Script | Description | Verification |
|-------------|-------------|--------------|
| [`memory_leak.sh`](../testing/integration/memory_leak.sh) | Tests for memory leaks in extended calculations | Memory usage remains stable |
| [`parallel.sh`](../testing/integration/parallel.sh) | Tests parallel execution capabilities | Successful execution, scalability |

---

## Regression Tests

Regression tests ensure that new code changes do not break existing functionality.

| Test Script | Description | Verification |
|-------------|-------------|--------------|

---

## Running Tests

### Prerequisites

- Compiled MOCCa executables in `exec/` directory
- Parameter files in `parameterizations/` directory
- Required parameterizations (e.g., SLy4, SLy5s1, t0t3, BSkG1, BSkG3) must be available

### Basic Usage

Most test scripts follow this pattern:
```bash
bash test_script.sh [EXESUFFIX] [PARAMETERS...]
```

Where `EXESUFFIX` is the suffix of the executable (e.g., "NLO" for `MOCCa.NLO.exe`).

### Running Specific Tests

**Unit test:**
```bash
cd testing/unit
bash fam_unit_tests.sh NLO --pairing HFB --parameterisation t0t3 --nw 15
```

**Integration test:**
```bash
cd testing/integration
bash minimal.sh NLO
bash symmetries.sh NLO NLO-P NLO-T NLO-TP
```


### Test Output

All tests write logs to their respective `logs/` subdirectories and use a temporary `work/` directory that is cleaned up automatically (unless the test fails).

---
