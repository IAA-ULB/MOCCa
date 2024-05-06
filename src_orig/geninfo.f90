module GenInfo
  !=============================================================================
  !_________ _______  _       _________ _______  _                 _______ 
  !\__   __/(  ___  )( (    /|\__   __/(  ___  )( \      |\     /|(  ____ \
  !   ) (   | (   ) ||  \  ( |   ) (   | (   ) || (      | )   ( || (    \/
  !   | |   | (___) ||   \ | |   | |   | (___) || |      | |   | || (_____ 
  !   | |   |  ___  || (\ \) |   | |   |  ___  || |      | |   | |(_____  )
  !   | |   | (   ) || | \   |   | |   | (   ) || |      | |   | |      ) |
  !   | |   | )   ( || )  \  |   | |   | )   ( || (____/\| (___) |/\____) |
  !   )_(   |/     \||/    )_)   )_(   |/     \|(_______/(_______)\_______)
  !                                                                       
  !  Copyright W. Ryssens & M. Bender
  !
  !=============================================================================
  ! Hephaestos keywords:
  ! 
  !   NUMSYM   : $NUMSYM
  !     Number of spatial symmetries, needed to calculate the volume element dv.
  !
  !   SYMSTRING: $SYMSTRING
  !   SYMLEN   : $SYMLEN
  !     String detailing the generators of the single-particle group desired.
  !     Only used for printing and human inspection.
  !
  !   REDUX    : $REDUX
  !   REDUY    : $REDUY
  !   REDUZ    : $REDUZ
  !     Integers indicating the reduction of the Cartesian axes.
  !
  !   LINESIZEX : $LINESIZEX
  !   LINESIZEY : $LINESIZEY
  !   LINESIZEZ : $LINESIZEZ
  !     Relations between nx/ny/nz and the actual number of mesh points in the 
  !     simulated system. This information is relevant here when using 
  !     periodic boundary conditions.
  !=============================================================================

  use compilation
  implicit none

  save

  character(len=$SYMLEN), parameter      :: SYMSTRING  = "$SYMSTRING"
  integer, parameter                     :: reduX      = $REDUX
  integer, parameter                     :: reduY      = $REDUY
  integer, parameter                     :: reduZ      = $REDUZ
  !-----------------------------------------------------------------------------
  !Number of points in every direction and total number of points
  integer :: nx=30,ny=30,nz=30, mv
  !-----------------------------------------------------------------------------
  !Number of protons and neutrons in the nucleus
  real(KIND=dp) :: Neutrons=10, Protons=10
  ! .... or alternatively a fixed chemical potential/fermi energy
  real(KIND=dp) :: mun = -10d8, mup = -10d8
  ! .... which is signalled by this particular flag
  logical       :: fixfermi = .false.
  !-----------------------------------------------------------------------------
  ! Line element and volume element of the box. dx is in fm, dv in fm^3.
  real(KIND=dp)  :: dx=0.8_dp
  real(KIND=dp)  :: dv=(0.8_dp**3)*(2**$NUMSYM)
  !-----------------------------------------------------------------------------
  ! Pi is always practical (delicious) to have.
  real(KIND=dp), parameter  :: pi=4.0_dp*atan2(1.0_dp,1.0_dp)
  !-----------------------------------------------------------------------------
  ! k_sh -- shift of the wavefunctions to be periodic
  real(KIND=dp) :: k_shx, k_shy, k_shz 
  !-----------------------------------------------------------------------------
  ! Maximum number of iterations and number of iterations to skip printing of
  ! the code in the evolve subroutine
  integer :: MaxIter=100, PrintIter=100
  !---------------------------------------------------------------------------
  ! This is the number of iterations during which the selfconsistent 
  ! potentials are not changed. The TOTAL number of iterations remains
  ! MaxIter; (MaxIter - FreezeIter) iterations DO change the potentials.
  integer :: FreezeIter = 0
  !---------------------------------------------------------------------------
  ! Coordinates of the mesh points for the calculation as well as the 
  ! coulomb calculation
  real(KIND=dp), allocatable         :: meshx(:), meshy(:), meshz(:)
  real(KIND=dp), allocatable, target :: meshgrid(:,:)
  ! Coordinates of the mesh points in the inertial frame of the nucleus, i.e.
  ! with the origin at the center-of-mass.
  real(KIND=dp), allocatable         :: meshx_shifted(:), meshy_shifted(:)     
  real(KIND=dp), allocatable         :: meshz_shifted(:)
  real(KIND=dp), allocatable, target :: meshgrid_shifted(:,:)  
  !-----------------------------------------------------------------------------
  ! Inverse temperature Beta = (k_b T)^{-1}.
  ! Negative values are used to indicate an infinite value, i.e. T = 0.
  real(KIND=dp) :: inversetemp = -1
  !-----------------------------------------------------------------------------
  ! Convergence criteria. See the module convergence for additional remarks.
  !
  !   Keyword         Default    Quantity
  !  ------------    ---------  -------------
  !   energy_prec     1d-9     abs((E^(i) - E^(i-1))/E^(i))     < energy_prec 
  !                            
  !   moment_prec     1d-3     abs((Q2m^(i) - Q2m^(i))/Q2m^(i)) < moment_prec
  !                                         if abs(beta_2m^(i)) > 0.01 
  !                            
  !
  !   disp_prec       1d-5     abs(sum_i v^2_i <psi|h^2|psi> - epsilon^2)
  !                                     < disp_prec
  !
  !   gradient_prec   1d+0     |s.p. gradient|  <    gradient_prec  
  !
  !   fermi_prec      1d-3     abs(lambda^(i) - lambda^(i-1)) < fermi_prec
  !                                    for both nucleon species
  !
  !   angmom_prec     1d-3     abs(<J_mu>^(i) - <J_mu>^(i-1)) < angmom_prec
  !                                    for all cartesian directions
  !-----------------------------------------------------------------------------
  real(KIND=dp) :: energy_prec = 1d-11, moment_prec = 1d-5, disp_prec = 1d-5
  real(KIND=dp) :: fermi_prec = 1d-3, angmom_prec   = 1d-3, gradient_prec=1d+0
  !-----------------------------------------------------------------------------
  ! Pairing tolerance
  ! Tolerance passed into the pairing solver. What exactly this determines 
  ! depends on the solver used, but for the default (Brent) solver, this 
  ! determines the relative precision on the Fermi energy itself.
  !
  ! Be very careful if you change this, as reducing this precision can lead to
  ! nonconverging calculations, especially when doing blocked calculations.
  real(KIND=dp) :: pairing_prec = 1d-15
  !-----------------------------------------------------------------------------
  ! Minimum number of iterations to perform before the code can stop itself
  ! for convergence detection. Default value = -1, in which case the multipole
  ! moments module modifies this number.
  integer :: min_iter_conv = -1
  !=============================================================================
  ! MPI parallelization variables
  !  NPROCS   = the number of MPI processes we are working with
  !  MPI_RANK = the rank of the current core
  ! Note that MPI_ranks are indexed starting at zero. 
  !
  ! NPROCS=1, MPI_RANK= 0 corresponds to a sequential calculation.
  !-----------------------------------------------------------------------------
  integer :: NPROCS = 1, MPI_RANK    = 0 
  !-----------------------------------------------------------------------------
  ! Load balancing strategy for the MPI ranks
  ! (0) : naive 1d block distribution of spwfs among ranks
  ! (1) : give entire symmetry blocks to MPI ranks
  integer :: balancing_strategy = 0
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Additional MPI communicator for the assigned symmetry block
  integer              :: MPI_COMM_BLOCK   ! communicator of the local team
  integer              :: MPI_SYM_BLOCK    ! assigned symmetry block
  integer              :: MPI_BLOCK_SIZE   ! number of spwfs for this block
  integer              :: MPI_BLOCK_RANK   ! rank inside the local team
  integer              :: MPI_BLOCK_NPROCS ! size of the local team
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! BLACS information for the communication between 1D and 2D grids
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! default BLACS context regrouping all ranks
  integer :: blacs_cntxt
  ! BLACS context for 1D spwf calculations within the current symmetry block
  integer :: blacs_cntxt_1D
  ! BLACS context for 2D spwf calculations within the current symmetry block
  integer :: blacs_cntxt_2D
  ! Size of the BLACS 2D layout within the current symmetry block
  integer :: NROW_2D, NCOL_2D, NROW_1D, NCOL_1D
  ! The coordinates of this MPI rank within the 2D BLACS layout
  integer :: MYROW_2D, MYCOL_2D, MYROW_1D, MYCOL_1D
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! SCALAPACK information
  ! blocking factors for rows and columns
  integer :: block_factor_row = 2 
  integer :: block_factor_col = 2 
  ! descriptors generated by DESCINIT
  integer :: desc_psi_1d(10) ! 1D layout of the spwfs
  integer :: desc_psi_2d(10) ! 2D block-cyclic layout of the spwfs
  integer :: desc_mat_2d(10) ! 2D block-cyclic layout of matrices in spwf-space
  !=============================================================================

contains

  subroutine ReadGenInfo(file_number)
    !---------------------------------------------------------------------------
    ! Read some of the general information needed.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   file_number : optional integer. If present, read from (open) channel
    !                 with this number. If absent, read from STDIN.
    !---------------------------------------------------------------------------
    integer(dp), intent(in), optional   :: file_number   
    integer                             :: io
#if(USE_MPI>0)
     integer                            :: mpi_err
#endif

    Namelist /nucleus/ neutrons,protons, inversetemp, mun, mup, fixfermi,      &
    &                  energy_prec, moment_prec, disp_prec, pairing_prec,      &
    &                  fermi_prec, balancing_strategy
    Namelist /mesh/    nx,ny,nz, dx

    if(MPI_rank .eq. 0) then    
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Reading the information on the nucleus by the first MPI rank
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      if(present(file_number)) then
        read (unit=file_number, nml=nucleus, iostat=io)
      else
        read (unit=*, nml=nucleus, iostat=io)
      endif

      ! Reading information on the mesh
      if(present(file_number)) then
        read (unit=file_number, nml=mesh)
      else
        read (unit=*, nml=mesh)
      endif   
      
      if((balancing_strategy .ne. 0) .and. (balancing_strategy .ne. 1)) then
        call stp('Invalid value for balancing_strategy.')
      endif

      if(fixfermi .and. (mun.eq.-10d8 .or.mup.eq.-10d8) )then
        call stp( 'You should fix an appropriate Lambda_N and Lambda_P.')
        stop
      endif
      
      ! Sanity check on the number of mesh points
      if(mod(nx,2) .ne. 0) then
        call stp('NX must be even')
      endif 
      if(mod(ny,2) .ne. 0) then
        call stp('NY must be even')
      endif
      if(mod(nz,2) .ne. 0) then
        call stp('NZ must be even')
      endif
    endif

#if(USE_MPI > 0)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Broadcasting all information from MPI_rank 0 to the rest
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! content of /mesh/ namelist
    call MPI_BCAST(nx, 1, MPI_integer, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(ny, 1, MPI_integer, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(nz, 1, MPI_integer, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(dx, 1, MPI_REAL8,   0, MPI_COMM_WORLD, mpi_err)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! content of /nucleus/ namelist
    ! a) particle numbers
    call MPI_BCAST(protons , 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(neutrons, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    ! b) fermi level options
    call MPI_BCAST(mun     , 1, MPI_REAL8  , 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(mup     , 1, MPI_REAL8  , 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(fixfermi, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_err)
    ! c) convergence parameters
    ! Note: convergence checking is likely to be done by a single MPI_RANK
    !       but this duplication just makes future programming errors 
    !       less likely.
    call MPI_BCAST(energy_prec , 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(moment_prec , 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(disp_prec   , 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(pairing_prec, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_BCAST(fermi_prec  , 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)

    call MPI_BCAST(balancing_strategy,1,MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
#endif
  
    ! Some bookkeeping operations, to be executed by all MPIranks 
    mv = nx * ny * nz
    dv = (dx**3)*(2**$NUMSYM)    
    ! NS: Shift of the wavefunctions applied
    k_shx=pi/($LINESIZEX*dx)
    k_shy=pi/($LINESIZEY*dx)
    k_shz=pi/($LINESIZEZ*dx)        
    
    call inimesh(meshx, meshy, meshz, nx, ny,nz, meshgrid,0.0d0,0.0d0,0.0d0)
        
  end subroutine ReadGenInfo

  function vector_product( mu ) result(indices)
    !---------------------------------------------------------------------------
    ! Function that returns the indices of the vector product with index mu
    ! meaning that in the expression 
    ! 
    !          (v1 x v2)_mu  = sum_(nu kappa) eps_{mu nu kappa} v1_nu v2_kappa
    !   
    ! It returns the indices of v1 and v2 on the rhs of this equation with the 
    ! positive Levi-Civita symbol. 
    !
    ! Hence, we can write
    !       i = indices(1)
    !       j = indices(2)
    !       (v1 x v2)_mu =  v1_i v2_j - v1_j v2_i
    !---------------------------------------------------------------------------
    integer :: indices(2)
    integer, intent(in) :: mu
        
    select case(mu)
    case(1)
      ! (v1 x v2)_x = v1_y v2_z - v1_z v2_y
      indices = (/ 2, 3 /)
    case(2)
      ! (v1 x v2)_y = v1_z v2_x - v1_x v2_z
      indices = (/ 3, 1 /)
    case(3)
      ! (v1 x v2)_z = v1_x v2_y - v1_y v2_z
      indices = (/ 1, 2 /)
    end select

  end function vector_product
  
  subroutine inimesh(x,y,z, mx, my, mz, mesh, shiftx, shifty, shiftz)
    !---------------------------------------------------------------------------
    ! Generate the coordinates of the mesh points for the Lagrange mesh.
    !
    ! Input: 
    !  mx,my,mz   : number of points on the mesh in every direction that need
    !               to be represented
    !  shiftx/y/z : Coordinates of the nuclear c.o.m. with respect to the box
    !               Defaults to (0,0,0)
    ! Output:
    !  x,y,z    : 1D arrays containing the coordinate values of the mesh points
    !  mesh     : 3D array containing the coordinate values of the mesh points
    !
    ! All distances in units of [fm].
    !---------------------------------------------------------------------------
    integer                                         :: i,j,k
    integer, intent(in)                             :: mx, my, mz
    real(KIND=dp), intent(in)                       :: shiftx, shifty, shiftz
    real(KIND=dp), intent(out), allocatable         :: x(:), y(:), z(:)
    real(KIND=dp), intent(out), allocatable, target :: mesh(:,:)

    real(KIND=dp)          :: startx, starty, startz
    real(KIND=dp), pointer :: gridx(:,:,:), gridY(:,:,:)  , gridZ(:,:,:)    

    allocate( x(mx), y(my),z(mz))
    allocate(mesh(mx*my*mz,3))
    
    if(reduX .eq.1) then
      startX = 1/2.0_dp
    else
      startX = -(mx/2-1/2.0_dp) - shiftx/dx     
                                ! divided by dx, because we will multiply after
    endif
    
    if(reduY .eq.1) then
      startY = 1/2.0_dp
    else
      startY = -(my/2-1/2.0_dp) - shifty/dx
                                ! divided by dx, because we will multiply after
    endif
    
    if(reduZ .eq.1) then
      startZ = 1/2.0_dp
    else
      startZ = -(mz/2-1/2.0_dp) - shiftz/dx
                                ! divided by dx, because we will multiply after
    endif
    
    do i=1,mx
      x(i) = (startx +(i-1))*dx
    enddo    

    do i=1,my
      y(i) = (starty +(i-1))*dx
    enddo
    do i=1,mz
      z(i) = (startz +(i-1))*dx
    enddo

    mesh = 0
    gridx(1:mx,1:my,1:mz) => mesh(1:mx*my*mz,1)
    gridy(1:mx,1:my,1:mz) => mesh(1:mx*my*mz,2)
    gridz(1:mx,1:my,1:mz) => mesh(1:mx*my*mz,3)

    do k=1,mz
      do j=1,my
        do i=1,mx
          gridx(i,j,k) = x(i)
          gridy(i,j,k) = y(j)
          gridz(i,j,k) = z(k)
        enddo 
      enddo
    enddo
    
  end subroutine inimesh
  
  integer function meshindex(i,j,k)
      !-------------------------------------------------------------------------
      ! The code relies on two types of mesh storage
      !   1) (i,j,k): indices used for arrays stored on a three-dimensional 
      !               mesh such as for example the Coulomb potential. 
      !               we have 1 <= i <= nx
      !                       1 <= j <= ny
      !                       1 <= k <= nz
      !   2) (i)    : one-dimensional indices that are used for efficiency
      !               to index the whole mesh.
      !                       1 <= i <= nx*ny*nz
      ! 
      ! This routine translates a set of indices (i,j,k) into the corresponding
      ! index in a one-dimensional mapping. 
      ! 
      ! Input:
      !   i,j,k : x/y/z mesh-indices in a three-dimensional mapping
      ! Output:
      !   meshindex : the equivalent index in a 1D mapping in FORTRAN order
      !               i+(j-1)*nx+(k-1)*ny*nx 
      !-------------------------------------------------------------------------
      integer, intent(in) :: i,j,k
      
      meshindex = i+(j-1)*nx+(k-1)*ny*nx
  
  end function meshindex

  subroutine find_nml_error(nmlname, iunit)
    !---------------------------------------------------------------------------
    ! Complain about an error in a namelist input, using the backspace command
    ! to find the offending line in an opened file. Note: this means this cannot
    ! be used to find errors in STDINPUT.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input: 
    !     nmlname : namelist name, to tell the user.
    !     iunit   : unit of the opened file to backspace.
    !
    !---------------------------------------------------------------------------
    character(len=*)    :: nmlname
    integer, intent(in) :: iunit
    
    character(len=1000) :: line

    backspace(iunit)
    read(iunit,fmt='(A)') line
    
    print *, '--------------------------------------------------------------'
    print *, 'Input problem encountered for namelist ', nmlname
    print *, 'This is the offending line:' 
    print *, ' > ', trim(line)
    print *, 'It likely contains a variable the code does not know about.'
    print *, '--------------------------------------------------------------'
    
    ! Stop the program
    call stp('')
  end subroutine find_nml_error 

  subroutine stp(msg, routine)
    !---------------------------------------------------------------------------
    ! A routine for stopping the entire code elegantly. For a single-core run
    ! it is somewhat trivial to type "print *, 'some error' ; stop" but this 
    ! does not translate well to MPI runs.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !     msg    : character, message to display.
    !     routine: character, optional. Name of the routine from which stp 
    !              was called.
    !---------------------------------------------------------------------------
    character(len=*), intent(in)           :: msg
    character(len=*), intent(in), optional :: routine

#if(USE_MPI>0)    
    integer :: mpi_err
    
    print *, 'RANK ', MPI_RANK, ' reports the following error: ', msg
#endif
    if(present(routine)) print *, "Error occurred in routine ", routine
#if(USE_MPI > 0)
    call MPI_ABORT(MPI_COMM_WORLD,1,mpi_err) ! force all MPI ranks to stop
                                             ! with error code 1
    !call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
    !call MPI_FINALIZE(mpi_err)
#endif
    stop ! simple stop

  end subroutine stp

  pure integer function LeviCivita(i,j,k)
    !---------------------------------------------------------------------------
    ! This function is a quick & dirty implementation of the LeviCivita symbol
    ! epsilon_{ijk}
    !---------------------------------------------------------------------------
    integer, intent(in) :: i,j,k

    if((i.eq.j).or.(j.eq.k).or.(k.eq.i)) then
        LeviCivita=0

    elseif(((i.eq.1).and.(j.eq.2).and.(k.eq.3)) &
     & .or.((i.eq.3).and.(j.eq.1).and.(k.eq.2)) &
     & .or.((i.eq.2).and.(j.eq.3).and.(k.eq.1))) then
        LeviCivita=1
    else
        LeviCivita=-1
    endif

    return
  end function LeviCivita
  
  function to_upper (str) result (string)
    !---------------------------------------------------------------------------
    ! Subroutine that changes a string to uppercase.
    !---------------------------------------------------------------------------
    character(*)        :: str
    character(len(str)) :: string

    Integer :: ic, i
    !Ugly but effective and independent of platform and implementation.
    character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    string = str
    do i = 1, len_trim(str)
    ic = INDEX(low, str(i:i)) !Note that ic = 0 when substring is not found
    if (ic > 0) then
    string(i:i) = cap(ic:ic)
    else
    string(i:i) = str(i:i)
    endif
    end do

  end function to_upper

  function to_lower (str) result (string)
    !---------------------------------------------------------------------------
    ! Subroutine that changes a string to lowercase.
    !---------------------------------------------------------------------------
    character(*)        :: str
    character(len(str)) :: string

    Integer :: ic, i
    !Ugly but effective and independent of platform and implementation.
    character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    string = str
    do i = 1, len_trim(str)
    ic = INDEX(cap, str(i:i)) !Note that ic = 0 when substring is not found
    if (ic > 0) then
    string(i:i) = low(ic:ic)
    else
    string(i:i) = str(i:i)
    endif
    end do

  end function to_lower

  function rps(string,length) result(r)
    !--------------------------------------------------------------------------
    ! function rps (right-padded-string) to add blancs to a string such that 
    ! it is printed left adjusted. Inspired by 
    ! http://computer-programming-forum.com/49-fortran/45c9683fdbd85176.htm
    !--------------------------------------------------------------------------
    character(len=*) :: string
    integer          :: length
    character(len=length) :: r

    r = adjustl(string)
  end function rps 

  subroutine clean_geninfo()
    !---------------------------------------------------------------------------
    ! Deallocate all allocated arrays, to exit in a clean fashion.
    !
    !---------------------------------------------------------------------------
  
    if(allocated(meshx)) then
      deallocate(meshx, meshy, meshz)
      deallocate(meshgrid)
    endif

  end subroutine clean_geninfo
end module GenInfo
