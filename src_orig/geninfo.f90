module GenInfo

    use CompilationInfo

    implicit none

    save
    !---------------------------------------------------------------------------
    !Number of points in every direction and total number of points
    ! Default = 16
    integer, public :: nx=16,ny=16,nz=16, nv=16**3
    !---------------------------------------------------------------------------
    !Number of protons and neutrons in the nucleus
    real(KIND=dp), public :: Neutrons, Protons
    !---------------------------------------------------------------------------
    ! Line element and volume element of the box. dx is in fm, dv in fm^3.
    real(KIND=dp), public  :: dx=0.8_dp
    real(KIND=dp), public  :: dv=(0.8_dp**3)
    !---------------------------------------------------------------------------
    !Time step in 10^-22 s
    real(KIND=dp), public  :: dt=0.012_dp, ReadjustTime=0.95_dp
    real(KIND=dp), parameter, public  :: pi=4.0_dp*atan2(1.0_dp,1.0_dp)
    !---------------------------------------------------------------------------
    !Maximum number of iterations and number of iterations to skip printing of
    ! the code in the evolve Subroutine
    integer, public :: MaxIter=0, PrintIter=0
    !---------------------------------------------------------------------------
    !Convergence parameters.
    real(KIND=dp),public  :: MomentPrec=1d-4, EnergyPrec=0.01E-08
    real(KIND=dp),public  :: PairingPrec=1d-4, CrankPrec=1d-4
    
contains

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
  
  subroutine to_upper (str, string)
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

  end subroutine to_upper
end module GenInfo
