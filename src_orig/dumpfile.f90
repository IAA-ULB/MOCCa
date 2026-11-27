!===============================================================================
! Debug dumpfile: standalone HDF5 snapshot of the full single-particle state,
! strictly more than what a restart needs. For debugging and external
! post-processing only; not read back by MOCCa.
!
! Content (see also write_mocca_hdf5, whose group layout this mirrors):
!   /                       run metadata via write_hdf5_attributes
!   /wavefunctions/states   spwfs on the mesh, (mv,4,nwt)
!   /wavefunctions/hpsi     action of h on every spwf, (mv,4,nwt)
!   /wavefunctions/dpsi     first derivatives, (mv,3,4,nwt)
!   /wavefunctions/ddpsi    second derivatives, (mv,6,4,nwt)
!   /fields/densities       mean-field densities via write_hdf5_densities
!   /fields/potentials      mean-field potentials via write_hdf5_potentials
!
! Derivatives are recomputed on the fly with apply_sphamil(onthefly=.true.);
! store_derivatives=.true. is not required in the input. Serial HDF5 only.
!===============================================================================

#if(USE_HDF5>0 && $FAM == 0)
module dumpfile
  use HDF5
  use geninfo,       only : dp, mv, stp
  use wavefunctions, only : HFpsi, sx, sy, sz, nwt, nwn
  use functional,    only : potentials, apply_sphamil, PotentialVector
  use densities,     only : Density, write_hdf5_densities
  use functional,    only : write_hdf5_potentials
  use IO_wf,         only : write_hdf5_attributes

  implicit none
  private
  public :: write_dumpfile

contains

  subroutine write_dumpfile(ofn)
    !---------------------------------------------------------------------------
    ! Writes the debug dumpfile: spwfs, h|psi>, derivatives, densities and
    ! potentials. Created or truncated.
    !
    ! Input:
    !   ofn : character, filename to write to.
    !---------------------------------------------------------------------------
    character(len=*), intent(in)   :: ofn

    integer(HID_T)                 :: file_id, group_id, subgroup_id
    integer(HID_T)                 :: space_id, plist_id, dset_id
    integer(hsize_t)               :: dims(3), dims4(4)
    integer                        :: wave, iso, h5ferr
    real(KIND=dp)                  :: hpsi_dump(mv,4,nwt)
    real(KIND=dp)                  :: dpsi_dump(mv,3,4,nwt), ddpsi_dump(mv,6,4,nwt)
    real(KIND=dp)                  :: dpsi(mv,3,4), ddpsi(mv,6,4)
$N3 real(KIND=dp)                  :: dddpsi(mv,10,4)
    type(PotentialVector)          :: F

#if(USE_MPI > 0)
    call stp('write_dumpfile: the debug dumpfile is not implemented for MPI runs.')
#endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Apply h to every spwf, collecting hpsi and the recomputed derivatives
    F = potentials
    do wave=1,nwt
      iso = -1
      if(wave .gt. nwn) iso = +1

      hpsi_dump(:,:,wave) = apply_sphamil(HFpsi(:,:,wave), dpsi, ddpsi,       &
$N3     &                                 dddpsi,                            &
      &                                 sx(:,wave), sy(:,wave), sz(:,wave),&
      &                                 iso, .true., F)

      dpsi_dump(:,:,:,wave)   = dpsi
      ddpsi_dump(:,:,:,wave)  = ddpsi
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Open the file and lay out the groups
    call h5open_f(h5ferr)
    call h5fcreate_f(ofn, H5F_ACC_TRUNC_F, file_id, h5ferr)
    if(h5ferr .ne. 0) call stp('Error creating dumpfile.')

    call h5gcreate_f(file_id, 'wavefunctions', group_id, h5ferr)            ; call h5gclose_f(group_id, h5ferr)
    call h5gcreate_f(file_id, 'fields', group_id, h5ferr)
    call h5gcreate_f(file_id, 'fields/densities' , subgroup_id, h5ferr)     ; call h5gclose_f(subgroup_id, h5ferr)
    call h5gcreate_f(file_id, 'fields/potentials', subgroup_id, h5ferr)     ; call h5gclose_f(subgroup_id, h5ferr)
    call h5gclose_f(group_id, h5ferr)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Run metadata on the root group
    call write_hdf5_attributes(file_id)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Spwf datasets: states, hpsi, derivatives
    dims = (/ mv, 4, nwt /)
    call h5screate_simple_f(3, dims, space_id, h5ferr)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, h5ferr)
    call h5dcreate_f(file_id, '/wavefunctions/states', H5T_NATIVE_DOUBLE,   &
    &                space_id, dset_id, h5ferr, plist_id)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, HFpsi, dims, h5ferr)
    if(h5ferr .ne. 0) call stp('Error writing states to dumpfile.')
    call h5dclose_f(dset_id, h5ferr)
    call h5dcreate_f(file_id, '/wavefunctions/hpsi', H5T_NATIVE_DOUBLE,     &
    &                space_id, dset_id, h5ferr, plist_id)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, hpsi_dump, dims, h5ferr)
    if(h5ferr .ne. 0) call stp('Error writing hpsi to dumpfile.')
    call h5dclose_f(dset_id, h5ferr)
    call h5pclose_f(plist_id, h5ferr)
    call h5sclose_f(space_id, h5ferr)

    dims4 = (/ mv, 3, 4, nwt /)
    call h5screate_simple_f(4, dims4, space_id, h5ferr)
    call h5dcreate_f(file_id, '/wavefunctions/dpsi', H5T_NATIVE_DOUBLE,     &
    &                space_id, dset_id, h5ferr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dpsi_dump, dims4, h5ferr)
    if(h5ferr .ne. 0) call stp('Error writing dpsi to dumpfile.')
    call h5dclose_f(dset_id, h5ferr)
    call h5sclose_f(space_id, h5ferr)

    dims4 = (/ mv, 6, 4, nwt /)
    call h5screate_simple_f(4, dims4, space_id, h5ferr)
    call h5dcreate_f(file_id, '/wavefunctions/ddpsi', H5T_NATIVE_DOUBLE,    &
    &                space_id, dset_id, h5ferr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ddpsi_dump, dims4, h5ferr)
    if(h5ferr .ne. 0) call stp('Error writing ddpsi to dumpfile.')
    call h5dclose_f(dset_id, h5ferr)
    call h5sclose_f(space_id, h5ferr)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Mean-field densities and potentials
    call write_hdf5_densities(file_id, Density)
    call write_hdf5_potentials(file_id, F)
    ! Convenience link, as in write_hdf5_fields
    call H5Lcreate_hard_f(file_id, '/fields/densities/D_I_I', file_id,      &
    &                     '/fields/rho', h5ferr)

    call h5fclose_f(file_id, h5ferr)
    call h5close_f(h5ferr)

  end subroutine write_dumpfile

end module dumpfile
#endif
