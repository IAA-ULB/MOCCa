submodule (functional) functional_io
!==============================================================================
!
! Implementation submodule of module functional: functional_io
!  Interface block lives in src_orig/functional.f90.
!
!==============================================================================
implicit none
contains

  subroutine WritePotentials(chan, F)
    !---------------------------------------------------------------------------
    !  Subroutine writing the different potentials to file.
    !---------------------------------------------------------------------------
    type(PotentialVector), intent(in) :: F
    integer, intent(in) :: chan
    integer             :: io

    ! Signalling how many potentials have been stored.
    write(chan, iostat=io) $POTENTIALNUMBER

    ! Then, for every potential write the
    ! * Name
    ! * Value
    ! Note that the name is written as a length-30 string, padded with spaces.
    ! If not, the unformatted in/out cannot correctly determine the end of a
    ! string and comparisons can not be made.
$WRITEPOTENTIALS
  end subroutine WritePotentials
 
#if(USE_HDF5>0 && $FAM == 0)
  ! No need to write potentials to file when dealing with FAM

  subroutine write_hdf5_potentials(file_id, F)
    use HDF5
    !---------------------------------------------------------------------------
    !  Subroutine writing the different potentials to the hdf5 file.
    !
    !  Note that all potentials will get stored as a completely flat array 
    !  because they can all have different ranks.
    !
    ! Input:
    !   file_id : integer(hid_t), identifier of the hdf5 file
    !   F       : potentialvector to write
    !---------------------------------------------------------------------------
    use HDF5
    use HDF5_auxiliary, only : hdf5_write_dataset_1d

    type(PotentialVector), intent(in) :: F
    integer(hid_t), intent(in)        :: file_id
    real(KIND=dp)                        :: pot(mv,4)

    ! Hephaestos fills in a call to hdf5_write_dataset_1d for every Skyrme potential
$WRITEPOTENTIALS_HDF5

    ! a little bit of extra work since the Coulomb potential is potentially
    ! defined on a bigger mesh
    pot = transfer_coulomb_mesh(F, .false.) 
    call hdf5_write_dataset_1d(file_id, 'coulomb_potential', pot, &
    &                     size(pot), groupname='fields/potentials')

  end subroutine write_hdf5_potentials
#endif


  function ReadPotentials(chan, filenx, fileny, filenz, sym_transfo_needed) &
  & result(F)
    !---------------------------------------------------------------------------
    ! Subroutine that reads the different mean-field potentials from a 
    ! wavefunction file created by a previous run of the code.
    ! Notes:
    !  1. this function does not rely on MPI I/O and simply reads everything 
    !     with rank 0 and then does a bunch of MPI_BCASTS.
    !  2. this function should not be confused with the readpotentials_separate
    !     routine below
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   chan                  : integer, channel number for input
    !   filenx, fileny,filenz : integers, number of mesh points in every
    !                           direction for the quantities on file
    !   sym_transfo_needed    : logical, if a symmetry transformation is
    !                           needed (.true.) or not (.false.)
    !                  .false.: use the potentials as read from
    !                           file, transforming only the number of mesh
    !                           points if needed.
    !                  .true. : use the potentials from file for further
    !                           calculations. This means just reading them here
    !                           and trusting the rest of the program to do the
    !                           the rest.
    ! Output:
    !   F                      : a potential-vector, read from file
    !---------------------------------------------------------------------------
    integer, intent(in)   :: chan, filenx, fileny, filenz
    logical, intent(in)   :: sym_transfo_needed
    type(PotentialVector) :: F, F_temp
    integer               :: io, potnumber, potcount, it, filemv
    character(len=30)     :: potname

#if(USE_MPI > 0)
    integer             :: mpi_err
#endif

    filemv = filenx * fileny * filenz

    ! Checking how many potentials have been stored on file
    if(MPI_RANK .eq. 0) read(chan, iostat=io) potnumber
#if(USE_MPI > 0)
    call MPI_BCAST(potnumber, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
#endif

    do potcount = 1,potnumber
        ! Read the name of the potential stored
        if(MPI_RANK .eq. 0) read(chan, iostat=io) potname
#if(USE_MPI > 0)
        call MPI_BCAST(potname,30, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpi_err)
#endif
        ! Then select which potential we are going to be reading
#if( $FAM == 0)        
        select case(trim(potname))
$READPOTENTIALS
        CASE DEFAULT
          ! The potential is not in this program, forget about it
          if(MPI_RANK .eq. 0 ) read(chan, iostat=io)
        end select
#else        
        ! Do NOT read potentials from file for FAM; these are real on file and complex in FAM!
        if(MPI_RANK .eq. 0 ) read(chan, iostat=io)
#endif 
    enddo

  end function ReadPotentials

#if(USE_HDF5 > 0)

  function read_potentials_hdf5(file_id, filenx, fileny, filenz, sym_transfo_needed) &
   & result(F)
    use HDF5
    !---------------------------------------------------------------------------
    ! Subroutine that reads the different mean-field potentials from a 
    ! wavefunction file created by a previous run of the code.
    ! Notes:
    !  1. now only sequental reading
    !  2. this function should not be confused with the readpotentials_separate
    !     routine below
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   chan                  : id of the group potentials
    !   filenx, fileny,filenz : integers, number of mesh points in every 
    !                           direction for the quantities on file
    !   sym_transfo_needed    : logical, if a symmetry transformation is 
    !                           needed (.true.) or not (.false.)
    !                  .false.: use the potentials as read from 
    !                           file, transforming only the number of mesh 
    !                           points if needed. 
    !                  .true. : use the potentials from file for further 
    !                           calculations. This means just reading them here
    !                           and trusting the rest of the program to do the
    !                           the rest.
    ! Output:
    !   F                      : a potential-vector, read from file
    !---------------------------------------------------------------------------
    use HDF5
    use HDF5_auxiliary

    integer(hid_t), intent(in)   :: file_id
    integer, intent(in)          :: filenx, fileny, filenz
    logical, intent(in)          :: sym_transfo_needed
    character(len=21), parameter :: groupname='/fields/potentials'

    integer(hid_t)             :: group_id
    type(PotentialVector)      :: F, F_temp
    real(kind=dp), allocatable :: Ftmp(:) 
    integer                    :: filemv, it, h5ferr
#if(USE_MPI > 0)
    integer                    :: mpi_err
#endif

    filemv = filenx * fileny * filenz
    
    if(MPI_RANK .eq. 0) then
      call h5gopen_f(file_id,groupname,group_id,h5ferr)
      if(h5ferr.ne.0) then
        write(*,*) 'Error opening group ',trim(groupname)
        call stp('ERROR: reading potentials in hdf5 format')
      endif
    endif

$READPOTENTIALS_HDF5

    if(MPI_RANK.eq.0) then
      !close the potentials group
      call h5gclose_f(group_id, h5ferr)
    endif

  end function read_potentials_hdf5


  subroutine hdf5_readpot(id, name, dset, n)
    use HDF5
    !-----------------------------------------------------------------
    ! Subroutine that reads a double precision potential of length n
    ! with some name in the hdf5 file.
    ! 
    ! Input:
    !   id    : integer(hid_t), identifier of the hdf5 file or group
    !   name  : character(len=*), name of the dataset to read
    !   n     : integer, length of the dataset to read
    !
    ! Output:
    !   dset  : real(kind=dp), intent(inout), array of length n to
    !           store the read dataset
    !-----------------------------------------------------------------
    use HDF5

    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: id
    integer, intent(in)          :: n
    real(kind=dp), intent(inout) :: dset(n)
    integer(hid_t)               :: dset_id !identifiers
    integer                      :: error
    integer(hsize_t), dimension(1) :: dims, data_dims

    dims=(/n/)
    data_dims(1)=n
    ! open dataset, "dset_id" is returned
    call h5dopen_f(id, name, dset_id, error)
    ! read dataset 
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dset, data_dims, error)
    ! Close access to dataset 
    call h5dclose_f(dset_id, error)

    if (error.ne.0) then
      call stp('ERROR: reading potentials in hdf5 format')
    endif

  end subroutine hdf5_readpot
#endif


  function readpotentials_separate(chan, ifn) result(F)
    !---------------------------------------------------------------------------
    ! Read mean-field potentials from a separate user-provided file that is 
    ! NOT a wavefunction file. This routine should not be confused with the 
    ! readpotentials subroutine from the functional module.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input: 
    !  * chan : integer, channel number to read the file
    !  * ifn  : input filename (will be checked for existence)
    !
    ! Output:
    !      F  : a set of mean-field potentials
    !
    ! Caution: this routine is currently foreseen for a specific application, 
    !          limited to maximally symmetric calculations and .func files
    !          for which F_Nm_Nm and G_I_NS potentials are defined.
    ! 
    !---------------------------------------------------------------------------
    use Coulombmod ! module explicitly 'used' in order to be able to place the 
                   ! values of the direct and exchange Coulomb potentials 
                   ! correctly on the mesh
  
    type(PotentialVector)        :: F
    integer, intent(in)          :: chan
    character(len=*), intent(in) :: ifn

    logical :: exists
    integer :: i,j,k,io,  mu, nu, ox, oy, oz, headercount
    real(KIND=dp), allocatable :: Vc(:), Ec(:)
    real(KIND=dp)              :: x,y,z
    character(len=200)         :: temp

    inquire(file=ifn, exist=exists)
    if(.not.exists) then
      print *, 'Input file specified does not exist!'
      stop
    endif
    
    open (chan,file=ifn)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! We need to skip any header lines (indicated by #).
    ! For a MOCCa-created file, there are 14 of them by default but other
    ! people might write a different amount
    io = 0; headercount = -1
    do while(io.eq.0) 
      headercount = headercount + 1
      read(chan, iostat=io, fmt='(a200)') temp
      if(temp(1:1) .ne. '#') io = 1
    enddo  
    ! We've found an error; we have counted the number of header lines!
    rewind(chan)
    ! ... and now we skip this number of lines    
    do i=1,headercount
        read(chan, fmt=('()'))
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Allocate the relevant potentials    
    allocate(F%F_I_I   (nx*ny*nz,4))     ; F%F_I_I    = 0.0d0
    allocate(F%FP_I_I  (nx*ny*nz,2))     ; F%FP_I_I   = 0.0d0
$TAUSCALAR    allocate(F%F_Nm_Nm(nx*ny*nz,4))   ; F%F_Nm_Nm = 0.0d0
$TAUTENSOR    allocate(F%F_N_N(nx*ny*nz,3,3,4)) ; F%F_N_N   = 0.0d0
$N1    allocate(F%G_I_NS (nx*ny*nz,3,3,4)) ; F%G_I_NS  = 0.0d0
    allocate(Vc(nx*ny*nz))              ; VC        = 0.0d0
    allocate(Ec(nx*ny*nz))              ; EC        = 0.0d0

    ! We assume the points on the file are correctly ordered in 
    ! FORTRAN fashion, such that we do not have to worry about looping 
    ! separately over x/y/z and can just loop once over all mesh points.
    ! This also means the coordinate information is not used.
    do i=1,nx*ny*nz
      read(chan, fmt='(7es25.12)', iostat=io, advance='no')        & 
      &                            x,y,z,                          & !unused
      &                            F%F_I_I(i,1),F%F_I_I(i,2),      & ! U(r)
      &                            Vc(i),  Ec(i)                     ! Coulomb
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! It is easy to read the kinetic potential from ETFSI calculations if 
      ! D_Nm_Nm is defined as a contracted density
$TAUSCALAR      read(chan, fmt='(2es25.12)', iostat=io, advance='no')    & 
$TAUSCALAR      &                            F%F_Nm_Nm(i,1), F%F_Nm_Nm(i,2)    ! kinetic
!      ! If not, then we have to do some reorganisation
$TAUTENSOR      read(chan, fmt='(2es25.12)', iostat=io, advance='no')    & 
$TAUTENSOR      &                            F%F_N_N(i,1,1,1), F%F_N_N(i,1,1,2)! kinetic
$TAUTENSOR      F%F_N_N(i,2,2,:) = F%F_N_N(i,1,1,:)/3
$TAUTENSOR      F%F_N_N(i,3,3,:) = F%F_N_N(i,1,1,:)/3
$TAUTENSOR      F%F_N_N(i,1,1,:) = F%F_N_N(i,1,1,:)/3
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Pairing field
      read(chan, fmt='(7es25.12)', iostat=io, advance='no')                    &
      &                            F%FP_I_I(i,1),F%FP_I_I(i,2)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! All components of the spin-orbit field
      do mu=1,3
        do nu=1,3
$N1          read(chan, fmt='(2es25.12)', advance='no', iostat=io) &
$N1          &               F%G_I_NS(i,mu,nu,1), F%G_I_NS(i,mu,nu,2)
        enddo
      enddo
      
      read(chan, *) ! Advance to new line
      
      if(io.ne.0) then
        print *, 'Problem encountered reading potential file ', ifn
        print *, 'IOSTAT = ', io
        stop
      endif
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Make sure the Coulomb module is configured with the right array dimensions
    call setup_coulomb(F)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! The index juggling is ugly, but necessary, because the Coulomb 
    ! potential has a different size than the Lagrange mesh.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ox = coul_offset_x ; oy = coul_offset_y ; oz = coul_offset_z
      
    do k=1,nz
      do j=1,ny
        do i=1,nx
          F%CoulombPotential(i+ox,j+oy,k+oz)  = Vc(meshindex(i,j,k))
          F%ExchangePotential(i+ox,j+oy,k+oz) = Ec(meshindex(i,j,k))
        enddo
      enddo
    enddo
    ! Ensure that the folded potentials are allocated if needed
    call Obtain_folded_potentials(F, sx_rho, sy_rho, sz_rho)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Make sure isospin combinations are made correctly for all potentials
    F%F_I_I(:,3) = F%F_I_I(:,1) + F%F_I_I(:,2)
    F%F_I_I(:,4) = F%F_I_I(:,1) - F%F_I_I(:,2)

$TAUSCALAR    F%F_Nm_Nm(:,3)   = F%F_Nm_Nm(:,1)   + F%F_Nm_Nm(:,2)
$TAUSCALAR    F%F_Nm_Nm(:,4)   = F%F_Nm_Nm(:,1)   - F%F_Nm_Nm(:,2)
$TAUTENSOR    F%F_N_N(:,:,:,3) = F%F_N_N(:,:,:,1) + F%F_N_N(:,:,:,2)
$TAUTENSOR    F%F_N_N(:,:,:,4) = F%F_N_N(:,:,:,1) - F%F_N_N(:,:,:,2)

$N1    do mu=1,3
$N1      do nu=1,3
$N1        F%G_I_NS(:,mu,nu,3) = F%G_I_NS(:,mu,nu,1) + F%G_I_NS(:,mu,nu,2)
$N1        F%G_I_NS(:,mu,nu,4) = F%G_I_NS(:,mu,nu,1) - F%G_I_NS(:,mu,nu,2)
$N1      enddo
$N1    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Close channel after succesfull IO operations.
    close(chan)
  end function readpotentials_separate


end submodule functional_io
