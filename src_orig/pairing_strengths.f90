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
 ! Hephaestos keywords:
 !
 !==============================================================================
 
 use geninfo
 use parameterization
 use pairingcutoffs

 implicit none
 
contains 

 function vmicro(rho, F_Nm_Nm, iso, ptype) 
  !-----------------------------------------------------------------------------
  ! Select the right routine for calculation the (position-dependent)
  ! microscopic pairing strength from among possible options. 
  ! Available right now:
  !
  !   ptype   subroutine
  !     0     Cao
  !
  ! Input:
  !   rho     : density   
  !   iso     : isospin (1 or 2 for neutrons or protons) 
  !   F_NM_NM : potential associated with D_Nm_Nm, for calculating
  !             the position-dependent effective mass.
  !   ptype   : select the type of microscopic pairing strength
  ! Output:
  !   vmicro: deduced pairing strength for both isospin species
  !-----------------------------------------------------------------------------
 
  integer, intent(in)       :: ptype, iso
  real(KIND=dp), intent(in) :: rho(mv,4), F_Nm_Nm(mv,4)
  real(KIND=dp)             :: vmicro(mv)
  
  call start_timer(T_microscopic_pairing)
 
  select case(ptype) 
  case(0) 
    vmicro = Cao(rho, F_Nm_Nm, iso)
  case DEFAULT
    print *, 'Unrecognized micptype value.'
    print *, 'Accepted values:'
    print *, '  0 : subroutine Cao'
    stop
  end select

  call stop_timer(T_microscopic_pairing)
 
 end function vmicro

 function Cao(rho, F_Nm_Nm, iso) result (vp)
  !-----------------------------------------------------------------------------
  ! Deduce the microscopic pairing strength (vp) from the density (rho).
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Input:
  !   rho  : density D_I_I, including BOTH isospins
  !   iso  : isospin component to return (1 or 2)
  ! F_Nm_Nm: potential associated with D_Nm_Nm, for calculating
  !          the position-dependent effective mass.
  ! Output:
  !   vp  : relevant pairing strength deduced, vp(r) [position-dependent!]
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in) :: rho(mv,4), F_Nm_Nm(mv,4)
  integer, intent(in)       :: iso
  integer                   :: i
  real(KIND=dp)             :: vp(mv), kf0(mv), kfp(mv), kfn(mv), eta(mv)
  real(KIND=dp)             :: deltann(mv), deltanp(mv), deltans(mv)
  real(KIND=dp)             :: delta(mv), x(mv), mu(mv), effm(mv)
  real(KIND=dp)             :: integral(mv)
 
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
  ! Calculating of gaps for each nucleon species
  delta =(1.-abs(eta))*deltans 
  select case(iso)
  case(1)
    delta = delta + eta * (eta + 1.0d0)/2.0d0 * deltann
  case(2)
    delta = delta + eta * (eta - 1.0d0)/2.0d0 * deltanp
  end select
  
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

 end function Cao

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

