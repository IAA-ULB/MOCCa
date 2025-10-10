# Development of a Python based MOCCa

At the 4th MOCCa user meeting wouter presented an example script of how a MOCCaPy script could look like. That seems like a good starting point for a top-down approach.

## f90 code overview

| .f90 file              | module             | use                                                                                 |
|------------------------|--------------------|-------------------------------------------------------------------------------------|
| basis_transform.f90    | basis_transform    | use geninfo                                                                         |
|                        |                    | use wavefunctions, only : HFblocks, nwt, spwf_map, nwt_local, hfblocks_global       |
|                        |                    | use timing                                                                          |
| BCS.f90                | BCS                | use vectors, only: PotentialVector                                                  |
|                        |                    | use wavefunctions                                                                   |
|                        |                    | use pairingcutoffs                                                                  |
| compilation.f90        | compilation        |                                                                                     |
| constants.f90          | constants          | use compilation                                                                     |
| convergence.f90        | convergence        | use Geninfo                                                                         |
| coulomb.f90            | Coulombmod         | use geninfo,          only: nx, ny, nz, dx, dp, stp                                 |
|                        |                    | use densities,        only: densityvector, potentialvector                          |
|                        |                    | use parameterization, only: dp, e2, pi, dv, coulorder, coultreatment                |
|                        |                    | use timing                                                                          |
| cranking.f90           | cranking           | use compilation                                                                     |
|                        |                    | use geninfo                                                                         |
|                        |                    | use derivatives                                                                     |
|                        |                    | use wavefunctions                                                                   |
|                        |                    | use nil8                                                                            |
|                        |                    | use pairing                                                                         |
|                        |                    | use densities                                                                       |
| densities.f90          | densities          | use compilation                                                                     |
|                        |                    | use geninfo                                                                         |
|                        |                    | use vectors                                                                         |
|                        |                    | use wavefunctions                                                                   |
|                        |                    | use pairing                                                                         |
|                        |                    | use derivatives                                                                     |
|                        |                    | use preconditioning                                                                 |
|                        |                    | use basis_transform                                                                 |
|                        |                    | use timing                                                                          |
| derivatives.f90        | derivatives        | use geninfo                                                                         |
| evolution.f90          | evolution          | use wavefunctions                                                                   |
|                        |                    | use functional                                                                      |
|                        |                    | use preconditioning                                                                 |
|                        |                    | use timing                                                                          |
| fam.f90                | fam                | use densities                                                                       |
|                        |                    | use moments                                                                         |
|                        |                    | use fission_MOI                                                                     |
|                        |                    | use evolution                                                                       |
| fam_gmres.f90          | gmres              | use geninfo, only : dp                                                              |
| fam_testing.f90        | fam_testing        | use densities                                                                       |
|                        |                    | use moments                                                                         |
|                        |                    | use Coulombmod, only : solve_coulomb                                                |
|                        |                    | use pairing,    only :  rho_can, pairingtype, rho_pairing, kappa_pairing            |
|                        |                    | use fission_MOI                                                                     |
|                        |                    | use functional                                                                      |
|                        |                    | use evolution                                                                       |
|                        |                    | use fam                                                                             |
|                        |                    | use gmres                                                                           |
| fission_MOI.f90        | fission_MOI        | use geninfo                                                                         |
|                        |                    | use densities                                                                       |
|                        |                    | use parameterization                                                                |
|                        |                    | use moments                                                                         |
|                        |                    | use timing                                                                          |
| folding.f90            | folding            | use geninfo                                                                         |
| functional.f90         | functional         | use compilation                                                                     |
|                        |                    | use geninfo                                                                         |
|                        |                    | use densities                                                                       |
|                        |                    | use parameterization                                                                |
|                        |                    | use pairing                                                                         |
|                        |                    | use transform                                                                       |
|                        |                    | use vectors                                                                         |
|                        |                    | use Cranking                                                                        |
|                        |                    | use pairing_strengths                                                               |
|                        |                    | use timing                                                                          |
|                        |                    | (use HDF5)                                                                          |
| geninfo.f90            | Geninfo            | use compilation , only : dp                                                         |
|                        |                    | (use MPI)                                                                           |
| hartree-fock.f90       | hartreefock        | use vectors, only: DensityVector, PotentialVector                                   |
|                        |                    | use wavefunctions                                                                   |
| hdf5_auxiliary.f90     | HDF5_auxiliary     | use HDF5                                                                            |
|                        |                    | use geninfo, dp, stp                                                                |
| HFB.f90                | HFB                | use geninfo                                                                         |
|                        |                    | use vectors, only: PotentialVector                                                  |
|                        |                    | use wavefunctions                                                                   |
|                        |                    | use pairingcutoffs                                                                  |
|                        |                    | use HFB_direct                                                                      |
|                        |                    | use HFB_gradient                                                                    |
| HFB_direct.f90         | HFB_direct         | use wavefunctions                                                                   |
|                        |                    | use parameterization                                                                |
| HFB_gradient.f90       | HFB_gradient       | use wavefunctions                                                                   |
| IO.f90                 | IO                 | use geninfo                                                                         |
|                        |                    | use wavefunctions                                                                   |
|                        |                    | use pairing                                                                         |
|                        |                    | use functional                                                                      |
|                        |                    | use momentsofinertia                                                                |
|                        |                    | use moments                                                                         |
|                        |                    | use Coulombmod                                                                      |
|                        |                    | use transform                                                                       |
|                        |                    | use fission_MOI                                                                     |
|                        |                    | use IO_wf, only: SYM_CODE, TRANS_CODE, allowtransform, extraspwfs                   |
|                        |                    | use IO_wf, only: version_number, file_version                                       |
|                        |                    | (use HDF5)                                                                          |
| IO_aux.f90             | IO_aux             |                                                                                     |
| IO_wf.f90              | IO_wf              | use compilation                                                                     |
|                        |                    | use GenInfo,       only: nx, ny, nz, dp, NPROCS, MPI_RANK, stp                      |
|                        |                    | use wavefunctions, only: HFBLOCKS                                                   |
|                        |                    | use functional,    only: ini_name_param, pairingtype, BCSGaps, HFBGaps, FermiEnergy |
|                        |                    | use transform,     only: sym_transfo_needed                                         |
| moments.f90            | moments            | use geninfo                                                                         |
|                        |                    | use sphericalharmonics                                                              |
|                        |                    | use Densities                                                                       |
| momentsofinertia.f90   | momentsofinertia   | use geninfo                                                                         |
|                        |                    | use densities                                                                       |
|                        |                    | use parameterization                                                                |
| nil8.f90               | nil8               | use compilation                                                                     |
| pairing.f90            | pairing            | use compilation                                                                     |
|                        |                    | use geninfo                                                                         |
|                        |                    | use wavefunctions                                                                   |
|                        |                    | use hartreefock                                                                     |
|                        |                    | use BCS                                                                             |
|                        |                    | use HFB                                                                             |
|                        |                    | use HFB_gradient                                                                    |
|                        |                    | use pairingcutoffs                                                                  |
|                        |                    | use parameterization                                                                |
|                        |                    | use pairing_strengths                                                               |
|                        |                    | use timing                                                                          |
|                        |                    | (use MPI)                                                                           |
| pairing_strengths.f90  | pairing_strengths  | use geninfo                                                                         |
|                        |                    | use parameterization                                                                |
| pairingcutoffs.f90     | pairingcutoffs     | use wavefunctions                                                                   |
| parameterization.f90   | parameterization   | use iso_fortran_env                                                                 |
|                        |                    | use pairingcutoffs                                                                  |
| precondition.f90       | preconditioning    | use derivatives                                                                     |
| printing.f90           | Printing           | use geninfo                                                                         |
|                        |                    | use pairing                                                                         |
|                        |                    | use wavefunctions                                                                   |
|                        |                    | use convergence                                                                     |
| scfiteration.f90       | SCFiteration       | use functional                                                                      |
|                        |                    | use densities                                                                       |
| sphericalharmonics.f90 | sphericalharmonics | use geninfo                                                                         |
| tantalus.f90           | Tantalus           | use geninfo                                                                         |
| timing.f90             | timing             | use compilation,     only : dp                                                      |
|                        |                    | use geninfo,         only : MPI_RANK, NPROCS, MPI_BLOCK_ASSIGNMENTS                 |
|                        |                    | use iso_fortran_env, only : int64, real64                                           |
|                        |                    | (use MPI           , only : MPI_COMM_WORLD, MPI_BARRIER)                            |
|                        |                    | use densities                                                                       |
| transform.f90          | transform          | use geninfo                                                                         |
|                        |                    | use wavefunctions                                                                   |
|                        |                    | use pairing                                                                         |
| vectors.f90            | vectors            | use geninfo                                                                         |
| version.f90            | version            | use GenInfo, only: SYMSTRING, MPI_RANK, NPROCS, reduX, reduY, reduZ                 |
|                        |                    | use IO,      only : SYM_CODE, TRANS_CODE                                            |
| wavefunctions.f90      | wavefunctions      | use derivatives                                                                     |
|                        |                    | use nil8                                                                            |
|                        |                    | (use MPI)                                                                           |

![module dependency graph](module_dependencies.png)