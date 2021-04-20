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
  !=============================================================================

  use compilation
  implicit none

  save

  character(len=$SYMLEN), parameter      :: SYMSTRING  = "$SYMSTRING"
  integer, parameter                     :: reduX      = $REDUX
  integer, parameter                     :: reduY      = $REDUY
  integer, parameter                     :: reduZ      = $REDUZ
  !---------------------------------------------------------------------------
  !Number of points in every direction and total number of points
  integer :: nx=30,ny=30,nz=30, mv
  !---------------------------------------------------------------------------
  ! Total number of single-particle wave-functions
  integer :: nwt=12
  !---------------------------------------------------------------------------
  !Number of protons and neutrons in the nucleus
  real(KIND=dp) :: Neutrons=10, Protons=10
  ! .... or alternatively a fixed chemical potential/fermi energy
  real(KIND=dp) :: mun = -10d8, mup = -10d8
  ! .... which is signalled by this particular flag
  logical       :: fixfermi = .false.
  !---------------------------------------------------------------------------
  ! Line element and volume element of the box. dx is in fm, dv in fm^3.
  real(KIND=dp)  :: dx=0.8_dp
  real(KIND=dp)  :: dv=(0.8_dp**3)*(2**$NUMSYM)
  !---------------------------------------------------------------------------
  ! Pi is always practical to have.
  real(KIND=dp), parameter  :: pi=4.0_dp*atan2(1.0_dp,1.0_dp)
  !---------------------------------------------------------------------------
  ! Coordinates of the mesh points for the calculation as well as the 
  ! coulomb calculation
  real(KIND=dp), allocatable         :: meshx(:), meshy(:), meshz(:)
  real(KIND=dp), allocatable, target :: meshgrid(:,:)
  !---------------------------------------------------------------------------
  ! Inverse temperature Beta = (k_b T)^{-1}.
  ! Negative values are used to indicate an infinite value, i.e. T = 0.
  real(KIND=dp) :: inversetemp = -1
  !---------------------------------------------------------------------------
  ! Convergence criteria
  !      Name       Default       Implementation   
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  !   energy_prec     1d-1     abs((E^(i) - E^(i-1))/E^(i))     < energy_prec
  !   moment_prec     1d-5     abs((Qlm^(i) - Qlm^(i))/Qlm^(i)) < moment_prec
  !                                if Qlm^(i) is large enough
  !   disp_prec       1d-5     abs(sum_i v^2_i <psi|h^2|psi> - epsilon^2)
  !                                     < disp_prec
  !   fermi_prec      1d-3     abs(lambda^(i) - lambda^(i-1)) < fermi_prec
  !                                    for both nucleon species
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  real(KIND=dp) :: energy_prec = 1d-11, moment_prec = 1d-5, disp_prec = 1d-5
  real(KIND=dp) :: fermi_prec = 1d-3
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
  ! Counter variables for the MPI implementation
  integer :: Counter = 1, Run = 1

contains

  subroutine ReadGenInfo(file_number)
    !---------------------------------------------------------------------------
    ! Read some of the general information needed.
    !---------------------------------------------------------------------------
    integer(dp), intent(in), optional   :: file_number   

    Namelist /nucleus/ neutrons,protons, inversetemp, mun, mup, fixfermi,      &
    &                  energy_prec, moment_prec, disp_prec, pairing_prec,      &
    &                  fermi_prec
    Namelist /mesh/    nx,ny,nz, dx
    
    ! Reading the information on the nucleus
    if(present(file_number)) then
      read (unit=file_number, nml=nucleus)
    else
      read (unit=*, nml=nucleus)
    endif

    if(fixfermi .and. (mun.eq.-10d8 .or.mup.eq.-10d8) )then
        print *, 'You should fix an appropriate Lambda_N and Lambda_P'
        stop
    endif

    ! Reading information on the mesh
    if(present(file_number)) then
      read (unit=file_number, nml=mesh)
    else
      read (unit=*, nml=mesh)
    endif   
    mv = nx * ny * nz
    dv = (dx**3)*(2**$NUMSYM)
    
    call inimesh(meshx, meshy, meshz, nx, ny,nz, meshgrid)
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
    case DEFAULT
      print *, 'Unknown value for mu in vector_prod.'
      stop
    end select

  end function vector_product
  
  subroutine inimesh(x,y,z, mx, my, mz, mesh)
    !---------------------------------------------------------------------------
    ! Generate the coordinates of the mesh points for the Lagrange mesh.
    ! Severe modification for Hephaestos will be necessary.
    !---------------------------------------------------------------------------
    integer                                         :: i,j,k
    integer, intent(in)                             :: mx, my, mz
    real(KIND=dp), intent(out), allocatable         :: x(:), y(:), z(:)
    real(KIND=dp), intent(out), allocatable, target :: mesh(:,:)

    real(KIND=dp), pointer :: gridx(:,:,:), gridY(:,:,:)  , gridZ(:,:,:)    

    allocate( x(mx), y(my),z(mz))
    allocate(mesh(mx*my*mz,3))
    
    do i=1,mx
      x(i) = (1/2.0_dp +(i-1))*dx
    enddo    
    do i=1,my
      y(i) = (1/2.0_dp +(i-1))*dx
    enddo
    do i=1,mz
      z(i) = (1/2.0_dp +(i-1))*dx
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
