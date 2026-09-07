# Breaking symmetries: practical aspects

We cover the basic workflow of running MOCCa with broken symmetries, i.e. in less-than-maximally-symmetric mode.

[ADD LINK TO CONFIG.md]

# Understanding MOCCa's maximally symmetric mode 

First some preliminaries on the maximally symmetric mode - i.e. the most used, essentially `default` mode. 






## Available symmetry modes

The code provides four distinct types of executables, each corresponding to a different combination of imposed symmetries:

| Executable | Source folder |ons | Typical use case | Spatial representation |
|------------|---------------|----------------|----------------|------------------|-------------------------|
| `MOCCa.XXX.exe` | `src/` | (Nx/2, Ny/2, Nz/2) | (nwn/2, nwp/2) | Starting point for most calculations | Full symmetry (PRS) |
| `MOCCa.XXX.T.exe` | `src_T/` | (Nx/2, Ny/2, Nz/2) | (nwn, nwp) | Nuclei with odd-N/Z | Full symmetry (PRS) |
| `MOCCa.XXX.P.exe` | `src_P/` | (Nx/2, Ny/2, Nz) | (nwn/2, nwp/2) | Octupole deformation and fission barriers | Reduced symmetry (RS) |
| `MOCCa.XXX.TP.exe` | `src_TP/` | (Nx/2, Ny/2, Nz) | (nwn, nwp) | Octupole for odd-Z/N | Reduced symmetry (RS) |

Here, `XXX` refers to your configuration name (e.g., `NLO` for a specific functional configuration).

*Note: The original manual included diagrams showing the spatial representation for different symmetry combinations (e.g., Boxes_PRS.eps for fully symmetric, Boxes_RS.eps for parity-broken modes). These visual aids may be helpful to include in future versions of this documentation.*

The symmetries referred to are:
- **Rz**: Rotational symmetry around the z-axis
- **T**: Time-reversal symmetry
- **P**: Parity (spatial reflection) symmetry
- **STy**: Simulated time-reversal symmetry (for odd systems)

## Compilation

All symmetry variants can be compiled using the main Makefile. To compile all variants:

```bash
make MC=mf CONFIG=YourConfig CALCTYPE=NUCLEI
```

To compile a specific symmetry variant, use the corresponding make target:

*Note: The exact make targets may vary between different MOCCa versions. The examples below follow the current MOCCa convention. In the original Tantalus code, the targets were slightly different (e.g., `make single`, `make single_T`, etc.).*

```bash
# Default (all symmetries conserved)
make single MC=mf CONFIG=YourConfig CALCTYPE=NUCLEI

# Time-reversal breaking
make single_T MC=mf CONFIG=YourConfig CALCTYPE=NUCLEI

# Parity breaking  
make single_P MC=mf CONFIG=YourConfig CALCTYPE=NUCLEI

# Both parity and time-reversal breaking
make single_TP MC=mf CONFIG=YourConfig CALCTYPE=NUCLEI
```

All compiled executables will be placed in the `exec/` directory.

## Symmetry transformation hierarchy

Practical constraints limit the code to:
1. Breaking only **one symmetry in any given run**
2. Being capable of breaking only **one single symmetry at a time**

This results in the following transformation hierarchy:

```
                        Maximally symmetric mode
                                     |
                   |----------<<<<<----------|---------->>>>>>----------|
                   |           Time-reversal broken                Parity broken          |
                   |                   |                                    |
                   |                   X                                    |
                   |                                                    |
                   |                          Parity & time-reversal broken       |
                   |__________________________________________________________|
```

The allowed transformation paths are:
- From **maximally symmetric** → **time-reversal broken**
- From **maximally symmetric** → **parity broken**
- From **parity broken** → **parity & time-reversal broken**

**Note**: The transformation from time-reversal broken to PT-broken is **NOT SUPPORTED** in the current repository.

## Warmstarting from more symmetric calculations

Wavefunction files from calculations with **more imposed symmetries** can be used to warmstart calculations with **less imposed symmetries**. This is particularly useful for:

- Starting from a converged symmetric calculation and then breaking a symmetry
- Reducing computation time by building upon previous results
- Systematically exploring symmetry breaking effects

### How to perform symmetry transformations

To transform from a more symmetric to a less symmetric calculation:

1. **Prepare your input file**: Start with a converged wavefunction file (`.wf`) from the more symmetric calculation
2. **Set the transformation flag**: Add `AllowTransform = .true.` to the `&IO` namelist
3. **Use the appropriate executable**: Run the calculation with the symmetry-breaking executable
4. **Ensure consistency**: Verify that the new degrees of freedom match the symmetry-broken version

### Transformation requirements

When breaking specific symmetries, the input parameters must satisfy certain conditions:

#### Breaking time-reversal symmetry

The total number of wavefunctions must be doubled:
```
nwn = 2 × nwn_file
nwp = 2 × nwp_file
```

#### Breaking parity symmetry

The total number of mesh points in the z-direction must be doubled:
```
nz = 2 × nz_file
```

## Important caveats

1. **Single operation constraint**: Only **one** of the following actions is allowed in any given run:
   - Modifying the number of mesh points
   - Adding single-particle wavefunctions
   - Breaking a single self-consistent symmetry

2. **Safety flag**: The `AllowTransform = .true.` flag removes several sanity checks on the input. It is **very easy** to make mistakes when this flag is activated.

3. **Wavefunction addition constraint**: When adding extra single-particle wavefunctions using the `extraspwfs` flag in the `&IO` namelist, ensure that `nwn = nwn_file + number of new neutron spwfs` and similarly for protons. The new wavefunctions are initialized randomly and set at very high energy.

4. **User responsibility**: When removing points from the mesh, **you are responsible** for verifying that the box remains sufficiently large. The code will not verify this automatically.

5. **Alternative symmetry breaking method**: Self-consistent symmetries can also be broken using the `&MomentConstraint` namelist with a specific `Iteration` parameter. If `Iteration` is set to a positive value (instead of -1), the constraint will only be applied for that number of iterations, which can be useful for breaking self-consistent symmetries.

## Practical workflow

### Example: Starting from symmetric and breaking parity

1. **Run symmetric calculation**:
   ```bash
   ./exec/MOCCa.NLO.exe < symmetric_input.in > symmetric_output.log
   ```

2. **Prepare parity-breaking input**: Create a new input file with:
   ```fortran
   &IO
     AllowTransform = .true.
     InWaveFunction = 'symmetric_result.wf'
     OutWaveFunction = 'parity_broken_result.wf'
     /&
   &mesh
     nz = 2 * nz_file  ! Double the z-direction points
     ! Other mesh parameters remain the same
     /&
   ```

3. **Run parity-breaking calculation**:
   ```bash
   ./exec/MOCCa.NLO.P.exe < parity_input.in > parity_output.log
   ```

### Example: From parity-broken to fully broken

1. **Starting from parity-broken wavefunction**: `parity_broken_result.wf`

2. **Prepare input for full symmetry breaking**:
   ```fortran
   &IO
     AllowTransform = .true.
     InWaveFunction = 'parity_broken_result.wf'
     OutWaveFunction = 'fully_broken_result.wf'
     /&
   ```

3. **Run fully broken calculation**:
   ```bash
   ./exec/MOCCa.NLO.TP.exe < fully_broken_input.in > fully_broken_output.log
   ```

## Verification

The code will always output information about the transformation performed in the header. Look for output like:

```
___________________________________________________________
| Transformation of the input                              |
| On file:                                                 |
|    nx, ny, nz    =    12   12   24                       |
|    dx            =  0.8000000 (fm)                       |
|    nwn, nwp      =      20     20                        |
|    (n+,n-)       = (   20    0    0    0)                |
|    (p+,p-)       = (   20    0    0    0)                |
| This calculation:                                        |
|    nx, ny, nz    =    12   12   28                       |
|    dx            =  0.8000000 (fm)                       |
|    nwn, nwp      =      20     20                        |
|    (n+,n-)       = (   20    0    0    0)                |
|    (p+,p-)       = (   20    0    0    0)                |
|__________________________________________________________|
```

This output helps you verify that the transformation was performed as expected.
