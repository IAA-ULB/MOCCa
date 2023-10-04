module pairing_strengths
 !==============================================================================
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
 !==============================================================================
 !
 ! 
 !
 !==============================================================================
 ! Hephaestos keywords: [NONE at the moment]
 !
 !==============================================================================
 
 use geninfo
 use parameterization
 use pairingcutoffs

 implicit none

 !------------------------------------------------------------------------------
 ! This abstract interface should, I believe, not be necessary to give context
 ! to our interpolation routines since there are explicit interfaces around. 
 ! GFORTRAN compiles fine, but for unknown reasons IFORT fails catastrophically
 ! when this is removed.
 abstract interface 
    function inter_abstract(deltans, deltann, deltanp, eta, iso) result(Delta)
      real*8, intent(in) :: deltans(:), deltann(:)
      real*8, intent(in) :: deltanp(:), eta(:)
      integer, intent(in):: iso 
      real*8, allocatable:: delta(:)
    end function
 end interface
 !------------------------------------------------------------------------------

contains 

 subroutine print_micro_pairing_info(ptype, intertype)
  !-----------------------------------------------------------------------------
  ! Print some information on the type of microscopic pairing type detected.
  !
  ! Input:
  ! ------- 
  !   ptype     : selection of gap to get to
  !   intertype :  selection of INM interpolation routine
  !
  !   ptype   type of gap     
  !   -----   -----------     
  !     0     Cao             
  !
  !
  ! intertype   interpolation routine         Interpolation approach
  ! ---------   ----------------------       ---------------------------------
  !  0          standard_interpolation  N. Chamel et al., PRC 80, 065804 (2009).
  !  1          linear interpolation    D=(1-|eta|)D_sym + |eta|D_{q,pure}
  !-----------------------------------------------------------------------------
  integer, intent(in) :: ptype, intertype
  
  1 format (' Microscopic treatment of the pairing active')
  2 format ('      Vmic prescription :  ', a20)
  3 format ('      INM interpolation :  ', a20)
 
  print 1
  select case(ptype)
  case(0)
    print 2, 'Cao et al., PRC 74 064301 (2006)'
  case DEFAULT
    call stp('PTYPE not recognized in pring_micro_pairing_info.')
  end select 

  select case(intertype)
  case(0)
    print 3, ' "Standard interpolation" from N. Chamel et al., PRC 80, 065804 (2009).'
    print *, '     ATTENTION: this INM interpolation is NOT recommended. '
  case (1)
    print 3, ' Linear interpolation '
    print *, '     Delta_q = (1 - |eta|) Delta_sym + |eta| Delta_{q,pure}'
  case DEFAULT
    call stp('intertype not recognized in print_micro_pairing_info.')
  end select 

 end subroutine print_micro_pairing_info

 function vmicro(rho, F_Nm_Nm, iso, ptype, intertype) 
  !-----------------------------------------------------------------------------
  ! Select the right routine for calculation the (position-dependent)
  ! microscopic pairing strength from among possible options. 
  ! 
  ! Input:
  !   rho     : density   
  !   iso     : isospin (1 or 2 for neutrons or protons) 
  !   F_NM_NM : potential associated with D_Nm_Nm, for calculating
  !             the position-dependent effective mass.
  !   ptype   : select the prescription for microscopic pairing strength
  !   intertype: select the prescription for INM matter interpolation
  !   
  ! Output:
  !   vmicro: deduced pairing strength for both isospin species
  !-----------------------------------------------------------------------------
 
  integer, intent(in)       :: ptype, iso, intertype
  real(KIND=dp), intent(in) :: rho(mv,4), F_Nm_Nm(mv,4)
  real(KIND=dp)             :: vmicro(mv)
  procedure(inter_abstract), pointer :: interpolation
  
  call start_timer(T_microscopic_pairing)
 
  ! Select the right type of interpolation
  select case(intertype)
  case(0)
    interpolation => standard_interpolation
  case(1)
    interpolation => linear_interpolation
  case DEFAULT
    call stp('Unrecognised option for intertype.')
  end select
 
  select case(ptype) 
  case(0) 
    vmicro = Cao(rho, F_Nm_Nm, iso, interpolation)
  case DEFAULT
    call stp('Unrecognized ptype option.')
  end select

  call stop_timer(T_microscopic_pairing)
 
 end function vmicro

 function Cao(rho, F_Nm_Nm, iso, interpolation, debug) result (vp)
  !-----------------------------------------------------------------------------
  ! Deduce the microscopic pairing strength (vp) from the density (rho).
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Input:
  !   rho         : density D_I_I, including BOTH isospins
  !   iso         : isospin component to return (1 or 2)
  ! F_Nm_Nm       : potential associated with D_Nm_Nm, for calculating
  !                 the position-dependent effective mass.
  !  interpolation: routine to deal with the interpolation to densities
  !                 that are not strictly symmetric or neutron matter.
  !  debug : if .true., print a ton of debugging output to stdout.
  ! Output:
  !   vp  : relevant pairing strength deduced, vp(r) [position-dependent!]
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in)  :: rho(mv,4), F_Nm_Nm(mv,4)
  integer, intent(in)        :: iso
  integer                    :: i
  real(KIND=dp)              :: vp(mv), kf0(mv), kfp(mv), kfn(mv), eta(mv)
  real(KIND=dp)              :: deltann(mv), deltanp(mv), deltans(mv)
  real(KIND=dp)              :: x(mv), mu(mv), effm(mv)
  real(KIND=dp)              :: integral(mv)
  real(KIND=dp), allocatable :: Delta(:)
  
  ! I originally coded this routine as taking a procedure as input. 
  ! Turns out that IFORT puts out catastrophic errors at some points...
  procedure(inter_abstract)  :: interpolation
 
  logical, optional, intent(in) :: debug
  logical                       :: debugflag

  1 format (' ----- Debugging function Cao ---------')
  2 format (' ISO on input =  ', i3)
  
  debugflag = .false.
  if(present(debug)) then
    if(debug) then
      debugflag = .true.
    endif
  endif


  if(debugflag) then
    ! Opening the debugging file
    print 1
    print 2, iso

    if(iso .eq. 2) then
      open(10,file='Cao.debug.p.out')
    else
      open(10,file='Cao.debug.n.out')
    endif
  endif
  

  vp = 0.0d0

  ! Fermi wavelengths
  kf0=(3.d0/2.0d0*pi**2*    rho(:,3))**(1./3.) ! Isoscalar density
  kfn=(3.d0      *pi**2*    rho(:,1))**(1./3.) ! Neutron density
  kfp=(3.d0      *pi**2*    rho(:,2))**(1./3.) ! Proton  density

  ! Asymmetry \eta
  do i=1,mv
    if(abs(rho(i,3)).gt.1d-12) then
      eta(i) = rho(i,4)/rho(i,3)
    else
      eta(i) = 0.0d0
    endif
  enddo
  
  deltann = Cao_delta(kfn, 1) ! pairing gap in pure neutron matter
  deltanp = Cao_delta(kfp, 2) !                pure proton  matter
  deltans = Cao_delta(kf0, 3) !                symmetric    matter

  ! Calculate the effective masses
  effm = hbm(iso) + F_Nm_Nm(:,iso) 
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  ! Calculating of gaps for each nucleon species using the interpolation 
  ! routine selected
  delta = interpolation(deltans, deltann, deltanp, eta, iso)  

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Calculation of the Fermi energies using the position-dependent 
  ! effective masses
  !
  ! mu_q = hbar^2/2M* k_f^2
  ! hbar^2/2M^* = hbar^2/2M + F_Nm_Nm(r)
  ! 
  select case(iso)
  case(1)
    mu = effm * kfn**2
  case(2)
    mu = effm * kfp**2
  end select   
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! calculation of the pairing strengths, inverted from the gaps. 
  ! The original way would be that of Eq. (1) of
  !     S. Goriely, N. Chamel and N. Pearson, PRL 102, 152503 (2009).
  ! which involves a numerical integral. 
  !
  ! The last development of the BSk-family is somewhat simpler: Eqs. (7-8) from 
  !   S. Goriely, N. Chamel and J. M. Pearson, PRC 93, 034337 (2016).
  ! where the integral is approximated as
  !
  ! I_q = \sqrt{\mu_q} [2 \log (2 \mu_q/\Delta_q) + \Lambda (\epsilon_l/\mu_q)]
  !
  ! where \mu_q is the INM approximation for the Fermi energy of species q, 
  ! \Delta_q is the gap for species q and epsilon_l is the pairing cutoff. 
  ! The function Lambda is 
  ! \Lambda(x) = log(16*x) + 2 * sqrt(1 + x) - 2 log(1 + sqrt{1 + x}) - 4.

  x = pairingcut(iso)/mu
  do i=1,mv
    if(delta(i) .gt. 0) then
      if(rho(i,iso) .gt. 1d-15) then
        integral(i) = sqrt(mu(i))* (2.d0*dlog(2.d0*mu(i)/delta(i))+Lambda(x(i)))
      else
        integral(i) = 2 * sqrt(pairingcut(iso))
      endif
    else
      integral(i) = 1d99
    endif
  enddo
  
!  if(iso .eq.1) then
!    open( unit=47 )
!    do i=1,nx
!      write(47,'(i3,1p,20e11.3)') iso, rho(i,1), mu(i), 6.5, eta(i), kf0(i), kfn(i), &
!      & delta(i), deltans(i), deltann(i), lambda(x(i)), integral(i)
!    enddo
!    write(47, '()')
!  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Final results for the pairing strengths
  vp = - (8.*pi**2)  /integral*(effm)**1.5d0 

  if(debugflag) then
    do i=1,mv
      write(10, fmt='(3f8.3, 13es25.12E3)') &
        &                        meshgrid(i,1), meshgrid(i,2), meshgrid(i,3),  &
        &                        eta(i), kf0(i), kfn(i), kfp(i), deltann(i),   &
        &                        deltanp(i), deltans(i), delta(i), mu(i), vp(i)
    enddo

    close(10)
  endif
 end function Cao

!-------------------------------------------------------------------------------
! Obtain an estimate for the microscopic gap by interpolating between the 
! microscopic gap in symmetric matter and the gap in pure matter. 
! The following functions are all different ways of doing the interpolation:
! 
! - standard_interpolation: used for the older BSk-models, like HFB31
!                           but this seems to be problematic, as we discovered
!                           while fitting BSkG3. NOT RECOMMENDED.
!
! - linear_interpolation  : simpler recipe to correct for the deficiencies of
!                           standard_inerpolation
!
 function standard_interpolation(deltans, deltann, deltanp, eta, iso) &
 &    result(Delta)
  !-----------------------------------------------------------------------------
  !  
  !     Delta_n = Delta_sym. (1 - abs(eta)) + eta (eta + 1)/2 Delta_{n, pure}
  !     Delta_p = Delta_sym. (1 - abs(eta)) + eta (eta - 1)/2 Delta_{p, pure}
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  !
  ! Input:
  !   deltans :  the gap in symmetric infinite matter
  !   deltann :  the gap in pure neutron matter
  !   deltanp :  the gap in pure proton matter
  !   eta     :  local asymmetry
  !   iso     :  for 1 (2) the code calculated the neutron (proton) gap by
  !              interpolation.
  ! Output:
  !   delta   : interpolated gap
  !
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in)  :: deltans(:), deltann(:), deltanp(:), eta(:)
  integer, intent(in)        :: iso 
  real(KIND=dp), allocatable :: Delta(:)
 
  Delta = Deltans * (1-abs(eta))
  select case(iso)
  case(1)
    ! neutrons
    Delta = Delta + eta * (eta + 1.0d0)/2.0d0 * deltann
  case(2)
    ! protons
    Delta = Delta + eta * (eta - 1.0d0)/2.0d0 * deltanp
  case DEFAULT
    call stp('Unrecognised input value for iso in standard_interpolation.')
  end select

 end function standard_interpolation
 
 function linear_interpolation(deltans, deltann, deltanp, eta, iso) &
 &    result(Delta)
  !-----------------------------------------------------------------------------
  !     Delta_n = Delta_sym. (1 - |eta|) + |eta| Delta_{n, pure}
  !     Delta_p = Delta_sym. (1 - |eta|) + |eta| Delta_{p, pure}
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  !
  ! Input:
  !   deltans :  the gap in symmetric infinite matter
  !   deltann :  the gap in pure neutron matter
  !   deltanp :  the gap in pure proton matter
  !   eta     :  local asymmetry
  !   iso     :  for 1 (2) the function calculates the neutron (proton) gap by
  !              interpolation.
  ! Output:
  !   delta   : interpolated gap
  !
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in)  :: deltans(:), deltann(:), deltanp(:), eta(:)
  integer, intent(in)        :: iso 
  real(KIND=dp), allocatable :: Delta(:)
 
  Delta = Deltans * (1-abs(eta))
  select case(iso)
  case(1)
    ! neutrons
    Delta = Delta + abs(eta) * deltann
  case(2)
    ! protons
    Delta = Delta + abs(eta) * deltanp
  case DEFAULT
    call stp('Unrecognised input value for iso in standard_interpolation.')
  end select

 end function linear_interpolation
!-------------------------------------------------------------------------------

 pure function Lambda(x) result(l)
  !-----------------------------------------------------------------------------
  ! Numerical function for the approximation of an integral in subroutine Cao.
  !
  ! \Lambda(x) = log(16*x) + 2 * sqrt(1 + x) - 2 log(1 + sqrt{1 + x}) - 4.
  !
  ! Input : 
  !    x: real argument
  !
  ! Output:
  !    L: Lambda(x)
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in) :: x
  real(KIND=dp), allocatable:: l
  
  l=dlog(16.d0*x)+2.d0*dsqrt(1.d0+x)- 2d0*dlog(1.d0+dsqrt(1.d0+x))-4.d0
  
 end function Lambda
 
 pure function Cao_delta(kf, iso) result(delta)
  !-----------------------------------------------------------------------------
  ! Calculation of gaps in infinite neutron, proton and symmetric matter.
  ! All expressions taken from the Brussels axial HFB code.
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in) :: kf(mv)
  integer, intent(in)       :: iso
  real(KIND=dp)             :: delta(mv), xkfint, xkfmax
  integer                   :: i
  
  select case(iso)
  case(1,2)
    ! Neutron matter and proton matter have identical gaps for this prescription
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xkfmax= 1.382d0
    xkfint= 1.25
    delta = 3.37968*kf**2/(kf**2+0.556092**2)*(kf-1.38236)**2 &
    &    /((kf-1.38236)**2+0.327517**2)

    do i=1, mv
      if (kf(i).gt.xkfint) delta(i)=0.39609*dexp(-(kf(i)-xkfint)/0.1)
    enddo

  case(3)
    ! Symmetric matter

    xkfmax=1.314d0
    xkfint=1.12
    delta =11.5586*kf**2*(kf-1.3142)**2/(kf**2+0.489932**2)/ & 
    &       ((kf-1.3142)**2+0.906146**2)
    do i=1,mv
      if (kf(i).ge.xkfint) delta(i)=0.42605d0*exp(-(kf(i)-xkfint)/0.1d0)
    enddo
  end select
    
 end function Cao_delta
 
end module pairing_strengths

