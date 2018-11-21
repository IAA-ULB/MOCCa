module GenInfo

    use compilation
    implicit none

    save
    !---------------------------------------------------------------------------
    !Number of points in every direction and total number of points
    integer :: nx=30,ny=30,nz=30, mv
    !---------------------------------------------------------------------------
    ! Total number of single-particle wave-functions
    integer :: nwt=12
    !---------------------------------------------------------------------------
    !Number of protons and neutrons in the nucleus
    real(KIND=dp)  :: Neutrons=10, Protons=10
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
    real(KIND=dp), allocatable :: meshx(:), meshy(:), meshz(:)
    real(KIND=dp), allocatable :: coulmeshx(:), coulmeshy(:), coulmeshz(:)

    !---------------------------------------------------------------------------
    ! Inverse temperature Beta = (k_b T)^{-1}.
    ! Negative values are used to indicate an infinite value, i.e. T = 0.
    real(KIND=dp) :: inversetemp = -1

    !---------------------------------------------------------------------------
    ! Convergence criteria
    !      Name       Default         
    !   energy_prec     1d-9     abs((E^(i) - E^(i-1))/E^(i))     < energy_prec
    !   moment_prec     1d-3     abs((Qlm^(i) - Qlm^(i))/Qlm^(i)) < moment_prec
    !                                if Qlm^(i) is large enough
    !   disp_prec       1d-5     abs(sum_i v^2_i <psi|h^2|psi> - epsilon^2)
    !                                     < disp_prec
    real(KIND=dp) :: energy_prec = 1d-11, moment_prec = 1d-5, disp_prec = 1d-7
    
contains

  subroutine ReadGenInfo
    !---------------------------------------------------------------------------
    ! Read some of the general information needed.
    !---------------------------------------------------------------------------
    Namelist /nucleus/ neutrons,protons, inversetemp
    Namelist /mesh/    nx,ny,nz, dx
    
    ! Reading the information on the nucleus
    read (unit=*, nml=nucleus)

    ! Reading information on the mesh
    read (unit=*, nml=mesh)
    
    mv = nx * ny * nz
    dv = (dx**3)*(2**$NUMSYM)
    
    call inimesh
  end subroutine ReadGenInfo
  
  subroutine inimesh
    !---------------------------------------------------------------------------
    ! Generate the coordinates of the mesh points for the demanded Lagrange mesh
    ! Severe modification for Hephaestos will be necessary.
    !---------------------------------------------------------------------------
    integer       :: i
    
    allocate(    meshx(nx  ),     meshy(ny  ),     meshz(nz  ))
    allocate(coulmeshx(nx+2), coulmeshy(ny+2), coulmeshz(nz+2))
    
    do i=1,nx
      meshx(i)     = (1/2.0_dp +(i-1))*dx
    enddo
    
    do i=1,nx+2
      coulmeshx(i) = (1/2.0_dp +(i-1))*dx
    enddo
    
    do i=1,ny
      meshy(i) = (1/2.0_dp +(i-1))*dx
    enddo
    
    do i=1,ny+2
      coulmeshy(i) = (1/2.0_dp +(i-1))*dx
    enddo
    
    do i=1,nz
      meshz(i) = (1/2.0_dp +(i-1))*dx
    enddo
    
    do i=1,nz+2
      coulmeshz(i) = (1/2.0_dp +(i-1))*dx
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

    if(len(str) .ne. len(string)) then
        stop('Strings of different length in to_upper')
    endif
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
end module GenInfo
