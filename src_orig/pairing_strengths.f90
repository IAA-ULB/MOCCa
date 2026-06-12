!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
! Written mainly by W. Ryssens & M. Bender
!
! Opensource software distributed under the GNU AGPLv3 licence, see the
!  LICENCE file in the root of this project.
!===============================================================================
module pairing_strengths
!==============================================================================
!
! This module deals with the determination of so-called microscopic pairing
! strengths as discussed in detail in
!    N. Chamel, PRC 82, 014313 (2010),
! and first used in a BSk model in
!    N. Chamel et al., NPA 812, 72 (2008)
! and originally (though seemingly not completely) formulated first in
!    T. Duguet et al., PRC 69, 054317 (2004).
!
! Large parts of the coding in this module reflect the variety of choices that
! can be made for the determination of the gaps. One can
!   (a) choose the recipe for the pairing gaps that should be reproduced
!   (b) choose an interpolation recipe for extending said gaps to INM at
!       arbitrary symmetry
!   (c) choose an integration strategy for calculating a rather complicated
!       integral that occurs in the formulas.
!
!==============================================================================
! Hephaestos keywords: [NONE at the moment]
!
!==============================================================================
 use geninfo
 use parameterization

 implicit none

 abstract interface 
    ! Interface for the different microscopic pairing gaps included
    function delta_abstract(kf, iso) result(delta)
      import :: mv, dp 
      real(KIND=dp), intent(in)  :: kf(mv)
      integer, intent(in)        :: iso
      real(KIND=dp)              :: delta(mv)
    end function delta_abstract
 end interface

 abstract interface 
      ! Interface for the different ways to interpolate the micros310copic pairing gaps
    function inter_abstract(delta_function, kfn, kfp, kf0, eta, iso) result(Delta)
      ! import statement to make this interface aware of the one above
      import                             :: delta_abstract, mv, dp
      procedure(delta_abstract), pointer :: delta_function
      real(KIND=dp), intent(in)          :: kfn(mv), kfp(mv), kf0(mv), eta(mv)
      integer, intent(in)                :: iso 
      real(KIND=dp)                      :: Delta(mv) 
    end function inter_abstract
 end interface

 abstract interface
    ! Interface for the different ways to calculate the integral in the inversion
    ! of the pairing strengths. 
    function integr_abstract(rho,mu,delta,u,cut) result(integral)
      import                    :: mv, dp
      real(KIND=dp), intent(in) :: rho(mv), mu(mv), delta(mv), u(mv), cut
      real(KIND=dp)             :: integral(mv)
    end function integr_abstract
 end interface

 interface vmicro
    ! Interface to traffic the work towards routines dealing with real or complex numbers.
    module procedure vmicro_real
    module procedure vmicro_complex
 end interface vmicro

 !------------------------------------------------------------------------------
 ! Global storage for the microscopically-derived pairing strengths.
 ! Since these numbers require some numerical integration to obtain and feature
 ! in quite a lot of places, it is worthwhile to calculate these values only
 ! once per iteration and store them.
 !------------------------------------------------------------------------------
 logical                    :: vmicro_stored(2) = .false.
 real(KIND=dp), allocatable :: vmicro_storage(:,:)
 !------------------------------------------------------------------------------

contains 

 subroutine print_micro_pairing_info(ptype, interpolationtype, integrationtype)
  !-----------------------------------------------------------------------------
  ! Print some information on the type of microscopic pairing type detected.
  !
  ! Input:
  ! ------- 
  !   ptype                  : selection of gap to get to
  !   interpolationtype      : selection of INM interpolation routine
  !   integrationtype        : selection of the integration routine
  !
  !   ptype   type of gap     
  !   -----   -----------     
  !     0     Cao             
  !
  !
  !             interpolation routine         Interpolation approach
  ! ---------   ----------------------       ---------------------------------
  !  0          standard_interpolation      Goriely et al. PRL 102, 152503 (2009)
  !  1          linear interpolation        D=(1-|eta|)D_sym + |eta|D_{q,pure}
  !  2          weak coupling interpolation [To be published]
  !
  !
  !             integration routine            Integration
  ! ---------   --------------------------    ---------------------------------
  !  0          weak_coupling_integration     analytical approximation
  !  1          tanh_sinh_integration         numerical integration
  !
  ! Note: the analytical formula employed for integrationtype=0 is not suited
  !       for N2LO forms.
  !
  !-----------------------------------------------------------------------------
  integer, intent(in) :: ptype, interpolationtype, integrationtype
  
  1 format (' Microscopic treatment of the pairing active')
  2 format ('      Vmic prescription   :  ', a50)
  3 format ('      INM interpolation   :  ', a50)
  4 format ('      Integral calculation:  ', a50)

  print 1
  select case(ptype)
  case(0)
    print 2, 'Cao et al., PRC 74 064301 (2006)'
  case DEFAULT
    call stp('PTYPE not recognized in pring_micro_pairing_info.')
  end select 

  select case(interpolationtype)
  case(0)
    print 3, ' "Standard interpolation" from N. Chamel et al., PRC 80, 065804 (2009).'
    print *, '  ATTENTION: this INM interpolation is NOT recommended. '
  case (1)
    print 3, ' Linear interpolation '
    print *, '  Delta_q = (1 - |eta|) Delta_sym + |eta| Delta_{q,pure}'
  case (2)
    print 3, ' Weak coupling interpolation'
    print *, '      Delta_n = Delta_NM(k_Fn)*[Delta_SM(k_F)/Delta_NM(k_F)]^{1-eta}'
    print *, '      Delta_p = Delta_NM(k_Fp)*[Delta_SM(k_F)/Delta_NM(k_F)]^{1+eta}'
  case DEFAULT
    call stp('interpolationtype not recognized in print_micro_pairing_info.')
  end select 

  select case(integrationtype)
  case(0)
    print 4, ' Analytical weak coupling approximation from N. Chamel, PRC 82, 014313 (2010).'
  case (1)
    print 4, ' Numerical integration by means of tanh-sinh quadrature.'
  case DEFAULT
    call stp('integrationtype not recognized in print_micro_pairing_info.')
  end select

 end subroutine print_micro_pairing_info

 function vmicro_complex(rho, U2, U4, iso, ptype, interpolationtype, integrationtype) result(vp)
  !-----------------------------------------------------------------------------
  ! Calculate a microscopically motivated (position-dependent) pairing strength
  !
  ! Input:
  !   rho               : density
  !   iso               : isospin (1 or 2 for neutrons or protons)
  !   U2          : potential multiplying k^2 in homogeneous INM
  !   U4          : potential multiplying k^4 in homogeneous INM
  !   ptype             : select the prescription for microscopic pairing strength
  !                      (0) gaps from BHF calculations by Cao et al.
  !   interpolationtype : select the prescription for INM matter interpolation
  !                      (0) "standard" interpolation of N. Chamel et al.
  !                      (1) linear interpolation of original BSkG3
  !                      (2) weak coupling interpolation as proposed by
  !                          N. Shchechilin.
  !   integrationtype   : select the type of integration to employ when
  !                     determining the microscopic pairing strengths
  !                      (0) analytical result for the weak coupling
  !                          approximation as proposed by N. Chamel.
  !                      (1) direct numerical integration through tanh-sinh
  !                          techniques.
  !
  ! Output:
  !   vp: deduced pairing strength for both species
  !-----------------------------------------------------------------------------
  integer, intent(in)          :: ptype, iso, interpolationtype, integrationtype
  complex(KIND=dp), intent(in) :: rho(mv,4)
  real(KIND=dp), intent(in)    :: U2(mv,4), U4(mv,4)
  real(KIND=dp)                :: vp(mv)

  vp = 0.0d0
  !call stp('VMICRO does not know how to handle complex densities yet.')
 end function vmicro_complex

 function vmicro_real(rho, U2, U4, iso, ptype, interpolationtype, integrationtype) &
 & result(vmicro)
  !-----------------------------------------------------------------------------
  ! Calculate a microscopically motivated (position-dependent) pairing strength
  !
  ! Technical note: both U2 and U4 are explicit inputs to this routine, because
  ! these functions are housed in the functional.f90 module; as the calculation
  ! of this INM potentials depends on the specific EDF employed in the
  ! calculation, this is the most logical place for them.
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  ! Input:
  !   rho               : density
  !   iso               : isospin (1 or 2 for neutrons or protons)
  !   U2          : potential multiplying k^2 in homogeneous INM
  !   U4          : potential multiplying k^4 in homogeneous INM
  !   ptype             : select the prescription for microscopic pairing strength
  !                      (0) gaps from BHF calculations by Cao et al.
  !   interpolationtype : select the prescription for INM matter interpolation
  !                      (0) "standard" interpolation of N. Chamel et al.
  !                      (1) linear interpolation of original BSkG3
  !                      (2) weak coupling interpolation as proposed by
  !                          N. Shchechilin.
  !   integrationtype   : select the type of integration to employ when
  !                     determining the microscopic pairing strengths
  !                      (0) analytical result for the weak coupling
  !                          approximation as proposed by N. Chamel.
  !                      (1) direct numerical integration through tanh-sinh
  !                          techniques.
  !
  ! Output:
  !   vmicro: deduced pairing strength for the requested species
  !-----------------------------------------------------------------------------
 
  integer, intent(in)       :: ptype, iso, interpolationtype, integrationtype
  real(KIND=dp), intent(in) :: rho(mv,4), U2(mv,4), U4(mv,4)
  real(KIND=dp)             :: vmicro(mv)

  procedure(inter_abstract), pointer  :: interpolation
  procedure(integr_abstract), pointer :: integration

  if(.not.allocated(vmicro_storage)) then
    vmicro_stored = .false.
    allocate(vmicro_storage(mv,2))
  endif

  ! Select the right type of interpolation
  select case(interpolationtype)
  case(0)
    interpolation => standard_interpolation
  case(1)
    interpolation => linear_interpolation
  case(2)
    interpolation => weak_coupling_interpolation
  case DEFAULT
    call stp('Unrecognised option for interpolationtype.')
  end select

  ! Select the right type of integration
  select case(integrationtype)
  case(0)
    integration => weak_coupling_integration
  case(1)
    integration => tanh_sinh_integration
  case DEFAULT
    call stp('Unrecognised option for integrationtype.')
  end select

  if(.not. vmicro_stored(iso)) then
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Perform the calculation of vmicro for this isospin
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call start_timer(T_microscopic_pairing)
    select case(ptype)
    case(0)
     vmicro_storage(:,iso) = Cao(rho, U2, U4, iso, interpolation,integration,.false.)
    case DEFAULT
     call stp('Unrecognized ptype option.')
    end select
    ! Signal to future calls that we have the numbers already
    vmicro_stored(iso) = .true.
    call stop_timer(T_microscopic_pairing)
  endif
  vmicro = vmicro_storage(:,iso)

 end function vmicro_real

 function Cao(rho, U2, U4, iso, interpolation, integration, debug) result (vp)
  !-----------------------------------------------------------------------------
  ! Deduce the microscopic pairing strength (vp) from the density (rho) at all
  ! mesh points, using the routine interpolation to obtain results away from
  ! pure matter and symmetric matter.
  !
  ! Note: this routine assumes that the single-particle hamiltonian in INM
  !       takes the following form:
  !
  !              epsilon_q = U_{0,q} + (hbar^2/2*m_q + U_{2,q}) k^2
  !                                  +                 U_{4,q}  k^4
  !
  !       in a notation close to that introduced by N. Chamel.
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  !
  ! Input:
  !   rho         : density D_I_I, including BOTH isospins
  !   iso         : isospin component to return (1 or 2)
  !   U2          : potential multiplying k^2 in homogeneous INM
  !   U4          : potential multiplying k^4 in homogeneous INM
  !  interpolation: routine to deal with the interpolation to densities
  !                 that are not strictly symmetric or neutron matter.
  !  integration  : routine to calculate the integral involved in the
  !                 calculatoin of the pairing strengths
  !  debug : if .true., print a ton of debugging output to stdout.
  !
  ! Output:
  !   vp  : relevant pairing strength deduced, vp(r) [position-dependent!]
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in)  :: rho(mv,4), U2(mv,4), U4(mv,4)
  integer, intent(in)        :: iso
  procedure(delta_abstract), pointer :: delta_function
  integer                    :: i
  real(KIND=dp)              :: vp(mv), kf0(mv), kfp(mv), kfn(mv), eta(mv)
  real(KIND=dp)              :: mu(mv), effm(mv), Delta(mv), u(mv), integral(mv)

  ! I originally coded this routine as taking a procedure as input. 
  ! Turns out that IFORT puts out catastrophic errors at some points...
  procedure(inter_abstract), pointer   :: interpolation
  procedure(integr_abstract), pointer  :: integration
 
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
  kf0=(3.d0/2.0d0*pi**2*    rho(:,3))**(1.0d0/3.0d0) ! Isoscalar density
  kfn=(3.d0      *pi**2*    rho(:,1))**(1.0d0/3.0d0) ! Neutron density
  kfp=(3.d0      *pi**2*    rho(:,2))**(1.0d0/3.0d0) ! Proton  density

  ! Asymmetry \eta
  do i=1,mv
    if(abs(rho(i,3)).gt.1d-12) then
      eta(i) = rho(i,4)/rho(i,3)
    else
      eta(i) = 0.0d0
    endif
  enddo

  ! Calculate the (local) effective mass
  effm = hbm(iso) + U2(:,iso)
  ! .... and the ratio between U2 and U4
  u = 4 * U4(:,iso)/(effm**2)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Calculating of pairing gaps for each nucleon species using the interpolation 
  ! routine selected
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  delta_function => cao_delta 
  delta = interpolation(delta_function, kfn, kfp, kf0, eta, iso)
  ! This routine seems to need to take a procedure POINTER; cray compilers 
  ! segfault if it is not a pointer.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Calculation of the Fermi energies using the position-dependent 
  ! effective masses
  !
  ! mu_q = hbar^2/2M* k_f^2 + U4 k_f^4
  ! hbar^2/2M^* = hbar^2/2M + U2
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  select case(iso)
  case(1)
    mu = effm * kfn**2 + U4(:,iso) * kfn**4
  case(2)
    mu = effm * kfp**2 + U4(:,iso) * kfp**4
  end select   
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! calculation of the complicated integral
  integral = integration(rho(:,iso), mu, delta, u, pairingcut(iso))
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Final results for the pairing strengths:
  !  see S. Goriely, N. Chamel and N. Pearson, PRL 102, 152503 (2009).
  vp = - (8.0d0*pi**2)  /integral*(effm)**1.5d0
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(debugflag) then
    do i=1,mv
      write(10, fmt='(3f8.3, 16es30.16E3)') &
        &           meshgrid(i,1), meshgrid(i,2), meshgrid(i,3),         &
        &           rho(i,1), rho(i,2), eta(i),                  &
        &           kf0(i), kfn(i), kfp(i), delta(i), mu(i),     &
        &           U2(i,iso), U4(i,iso), u(i), integral(i), vp(i)
    enddo
    close(10)
  endif
 end function Cao

 function integrand(xi, mu, delta, u) result(I)
    !----------------------------------------------------------------------------
    ! Integrand of the complicated integral
    !
    !                  sqrt (xi)
    !  I (xi) = --------------------------- sqrt(2(1+uxi)(1+sqrt(1+uxi)))^-1
    !           sqrt((xi - mu)^2 + Delta^2)
    !
    ! Input:
    !    xi   : integration variable
    !    mu   : reduced Fermi energy
    !    delta: targetted pairing gap
    !    u    : 4 * U4/U2**2, see above
    !----------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: xi, mu, delta, u
    real(KIND=dp)             :: I, Eqp, N2LOfac

    Eqp    = sqrt((xi-mu)**2 + delta**2)
    N2LOfac= 2.0d0/( ( 1 + u * xi) * (1 + sqrt(1+ u * xi)))
    N2LOfac= sqrt(N2LOfac)
    I   = sqrt(xi)/Eqp * N2LOfac

 end function integrand

 function tanh_sinh(mu, delta, u , a, b, h, tol) result(I)
    !----------------------------------------------------------------------------
    ! Calculate the following integral numerically:
    !
    !    I = \int_{0}^{\mu + \epsilon_lambda} integrand(xi, mu, delta, u) dxi
    !
    ! where the function integrand is defined above. This is not trivial, because
    ! the integrand is a very peaked function at \mu if \Delta is small, which
    ! typically happens at very low or very high densities. This routine employs
    ! a tanh-sinh quadrature (also known as a double exponential formula): the
    ! main interest of this technique is its robustness with respect to
    ! integrable singularities on the borders of the integration range.
    !
    ! This routines strategy is first to split the integral
    !
    !   I =  \int_{0}^{\mu} integrand(xi, mu, delta) dxi
    !     +  \int_{mu}^{\mu+\epsilon_lambda} integrand(xi, mu, delta) dxi
    !
    ! and apply a tanh-sinh change of variables to each:
    !
    !    xi = tanh(pi/2 * sinh(t))
    !
    ! which can changes
    !
    !     int_{-1, 1} f(x) dx = \int_{-\infty}^{+\infty} f(x(t)) dx/dt dt
    !
    ! Of course, one needs first a linear transformation to recast each original
    ! integral as one restricted to [-1,1]. The advantage of the r.h.s. is that
    ! the extremely quickly decreasing dx/dt means we don't need collocation
    ! points at very large values of t and kills any contribution from the
    ! points near \mu.
    !
    ! In practice, we take symmetric collocation points at
    !
    !      t_i = (i + 0.5) * h    and t_-i =-(i + 0.5) * h
    !
    ! for i = 1, N. N is adaptively determined by simply adding points until the
    ! contribution to the integral becomes small enough.
    !
    ! This routine is somewhat stupid in that it uses a fixed stepsize h and offers
    ! no error control in practice; this could be improved easily through refinement
    ! with additional points.
    !
    ! This implementation is strongly inspired by the one in the mpmath python library.
    ! For more information see:
    !   - Numerical Recipes in C++, section 4.4 in the second edition.
    !   - https://en.wikipedia.org/wiki/Tanh-sinh_quadrature
    !   - Bailey, David H, "Tanh-Sinh High-Precision Quadrature". (2006).
    !     https://www.davidhbailey.com/dhbpapers/dhb-tanh-sinh.pdf
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !
    !  mu    : reduced Fermi energy             | arguments of the function integrand
    !  delta : targetted pairing gap            | defined above
    !  u     : ratio of NLO and N2LO potentials |
    !  a, b  : limits of the integration interval
    !  h     : step size of the collocation points
    !  tol   : tolerance determining when to stop adding points
    !
    ! Output:
    !
    !   I    : value of the integral
    !
    !----------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: a, b, h, tol, mu, delta, u
    real(KIND=dp)             :: I, C, D, fxm, fxp, t, t0, x0, w, w0, xm, xp
    integer                   :: k

    I = 0
    k = -1
    C = (b-a)/2 ;  D = (b+a)/2
    fxm = 10 ; fxp = 10
    w  = 10
    t0 = h/2

    do while (abs(fxm*w*h) .gt. tol .or. abs(fxp*w*h) .gt. tol)
        k = k +1
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Add two integration points at +(t0 + k* h) and -((t0 + k* h)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        t  = t0 + k * h
        ! this is the coordinate and weight for integration on the interval [-1,1]
        x0 = tanh(pi/2 * sinh(t))
        w0 = pi/2 * cosh(t) / cosh(pi/2 * sinh(t))**2

        ! We transform x and w to the [a,b] interval with a linear transformation
        w   = C*w0
        xp  = C*x0 + D      ; xm  =-C*x0 + D

        ! Function evaluation for both points
        fxp = integrand(xp, mu, delta, u)
        fxm = integrand(xm, mu, delta, u)

        ! ... and we add both points to the integral with the appropriate weight
        I   = I + (fxp+fxm)*w*h
    enddo
 end function tanh_sinh

 function tanh_sinh_integration(rho, mu, delta, u, cut) result (integral)
   !----------------------------------------------------------------------------
   ! Use tanh-sinh quadrature to evaluate the complicated integral of the
   ! function integrand defined above with the strategy described for the
   ! routine tanh_sinh.
   !
   ! The settings h = 0.2 MeV and tol=2d-4 have been manually optimised to
   ! result in a relative error below 1-e3 for rho in [0, 0.4] fm^{-3} when
   ! using the gaps of Cao et al. as target with a feasible number of points,
   ! typically about 50. Of course, this means a superior accuracy (up to 1e-7)
   ! for the densities encountered in finite nuclei.
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Input:
   !   rho      : density of the nucleon species under consideration
   !   mu       : reduced Fermi energy of the species under consideration
   !   delta    : targetted pairing gap
   !   u        : ratio of NLO and N2LO potentials
   !   cut      : value of the pairing cutoff
   !
   ! Output :
   !   integral : numerical value of the integral
   !----------------------------------------------------------------------------
   real(KIND=dp), intent(in) :: rho(mv), mu(mv), delta(mv), u(mv), cut
   integer                   :: i
   real(KIND=dp)             :: integral(mv)

   do i=1, mv
    if(delta(i) .gt. 0) then
      ! numerical safeguard for underflowing \Delta
      if(rho(i) .gt. 1d-15) then
        ! Numerical safeguard for very low or negative density
        integral(i) =&
        &       tanh_sinh(mu(i), delta(i), u(i), 0.0d0, mu(i)    , 0.2d0, 2d-4)  &
        &     + tanh_sinh(mu(i), delta(i), u(i), mu(i), mu(i)+cut, 0.2d0, 2d-4)
      else
        ! Analytical limit of rho and delta tending to zero
        integral(i) = 2 * sqrt(cut)
      endif
    else
        integral(i) = 1d99
    endif
  enddo

 end function tanh_sinh_integration

 function weak_coupling_integration(rho,mu, delta, u, cut) result (integral)
  !----------------------------------------------------------------------------
  ! Employ an analytical formula to approximate the following integral
  !
  !    I = \int_{0}^{\mu + \epsilon_lambda} integrand(xi, mu, delta, u=0) dxi
  !
  ! where the function integrand is defined above. The analytical approximation
  ! is Eqs. (7-8) from
  !
  !   S. Goriely, N. Chamel and J. M. Pearson, PRC 93, 034337 (2016).
  !
  ! and reads
  !
  !  I_q = \sqrt{\mu_q} [2 \log (2 \mu_q/\Delta_q) + \Lambda (\epsilon_l/\mu_q)]
  !
  ! where \mu_q is the INM approximation for the Fermi energy of species q,
  ! \Delta_q is the gap for species q and epsilon_l is the pairing cutoff.
  ! The function Lambda is
  !   \Lambda(x) = log(16*x) + 2 * sqrt(1 + x) - 2 log(1 + sqrt{1 + x}) - 4.
  !
  ! Note that this approximation has so far not been generalized to N2LO EDFs.
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Input:
  !   rho      : density
  !   mu       : reduced Fermi lveel
  !   delta    : targetted pairing gaps
  !   cut      : pairing cutoff
  !
  ! Output : 
  !   integral : valueof the integral
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  real(KIND=dp), intent(in) :: rho(mv), mu(mv), delta(mv), cut, u(mv)
  real(KIND=dp)             :: integral(mv),x(mv)
  integer :: i

  if(any(abs(u) .gt. 1.0d-16)) then
    call stp('One should not combine integrationtype=0 with an N2LO EDF.')
  endif

  x = cut/mu
  do i=1,mv
    if(delta(i) .gt. 0) then
      if(rho(i) .gt. 1d-15) then
        ! Numerical safeguard for very low or negative density
        integral(i) = sqrt(mu(i))* (2.d0*dlog(2.d0*mu(i)/delta(i))+Lambda(x(i)))
      else
        ! Analytical limit of rho and delta tending to zero
        integral(i) = 2 * sqrt(cut)
      endif
    else
      integral(i) = 1d99
    endif
  enddo

 end function weak_coupling_integration

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
!                           standard_interpolation.
!
! - weak_coupling_interpolation: interpolation scheme designed by N. Shchechilin 
!                                to match better the asymmetry-dependence 
!                                predicted in: 
!                                 S.S. Zhan et al., PRC 81(4):044313. (2010)
 function standard_interpolation(delta_function, kfn, kfp, kf0, eta, iso) &
 &                               result(Delta)
  !-----------------------------------------------------------------------------
  !
  !  Delta_n = Delta_SM(k_F) (1 - abs(eta)) + eta (eta + 1)/2 Delta_NM(k_Fn)
  !  Delta_p = Delta_SM(k_F) (1 - abs(eta)) + eta (eta - 1)/2 Delta_NM(k_Fp)
  !
  ! S. Goriely et al. PRL 102, 152503 (2009)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  !
  ! Input:
  !   delta_function :  procedure pointer to the routine giving the gaps as 
  !                     a function of Fermi wavenumber
  !   kfn            :  Fermi wavenumber for neutrons
  !   kfp            :  Fermi wavenumber for protons
  !   kf0            :  Fermi wavenumber for total density
  !   eta            :  local asymmetry
  !   iso            :  for 1 (2) the function calculates the neutron (proton) 
  !                     gap by interpolation.
  ! Output:
  !   delta   : interpolated gap
  !
  !-----------------------------------------------------------------------------
  procedure(delta_abstract), pointer :: delta_function
  real(KIND=dp), intent(in)          :: kfn(mv), kfp(mv), kf0(mv), eta(mv)
  integer, intent(in)                :: iso 
  real(KIND=dp)                      :: Delta(mv) 
  real(KIND=dp)                      :: deltann(mv), deltanp(mv), deltans(mv)
  integer                            :: i
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Calculate all pairing gaps that occur in the interpolation formula
  deltann = delta_function(kfn, 1) ! pairing gap in pure neutron matter at k_Fn
  deltanp = delta_function(kfp, 2) !                pure proton  matter at k_Fp
  deltans = delta_function(kf0, 3) !                symmetric    matter at k_F
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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

  ! Safeguard: if this results in negative pairing gaps, just make pairing vanish
  do i=1,mv
    if(Delta(i) .lt. 1e-14) Delta(i) = 0.0d0
  enddo

 end function standard_interpolation
 
 function linear_interpolation(delta_function, kfn, kfp, kf0, eta, iso) &
 &    result(Delta)
  !-----------------------------------------------------------------------------
  !     Delta_n = Delta_SM(k_F) (1 - |eta|) + |eta| Delta_NM(k_Fn)
  !     Delta_p = Delta_SM(k_F) (1 - |eta|) + |eta| Delta_NM(k_Fp)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  !
  ! Input:
  !   delta_function :  procedure pointer to the routine giving the gaps as 
  !                     a function of Fermi wavenumber
  !   kfn            :  Fermi wavenumber for neutrons
  !   kfp            :  Fermi wavenumber for protons
  !   kf0            :  Fermi wavenumber for total density
  !   eta            :  local asymmetry
  !   iso            :  for 1 (2) the function calculates the neutron (proton) 
  !                     gap by interpolation.
  ! Output:
  !   delta   : interpolated gap
  !
  !-----------------------------------------------------------------------------
  procedure(delta_abstract), pointer:: delta_function
  real(KIND=dp), intent(in)         :: kfn(mv), kfp(mv), kf0(mv), eta(mv)
  integer, intent(in)               :: iso 
  real(KIND=dp)                     :: Delta(mv) 
  real(KIND=dp)                     :: deltann(mv), deltanp(mv), deltans(mv)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Calculate all pairing gaps that occur in the interpolation formula
  deltann = delta_function(kfn, 1) ! pairing gap in pure neutron matter at k_Fn
  deltanp = delta_function(kfp, 2) !                pure proton  matter at k_Fp
  deltans = delta_function(kf0, 3) !                symmetric    matter at k_F
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  Delta = Deltans * (1-abs(eta))
  select case(iso)
  case(1)
    ! neutrons
    Delta = Delta + abs(eta) * deltann
  case(2)
    ! protons
    Delta = Delta + abs(eta) * deltanp
  case DEFAULT
    call stp('Unrecognised input value for iso in linear_interpolation.')
  end select

 end function linear_interpolation
 
 function weak_coupling_interpolation(delta_function, kfn, kfp, kf0, eta, iso) &
 &    result(Delta)
  !-----------------------------------------------------------------------------
  !   Delta_n = Delta_NM(k_Fn)*[Delta_SM(k_F)/Delta_NM(k_F)]^{1-eta} 
  !   Delta_p = Delta_NM(k_Fp)*[Delta_SM(k_F)/Delta_NM(k_F)]^{1+eta}
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Input:
  !   delta_function :  procedure pointer to the routine giving the gaps as 
  !                     a function of Fermi wavenumber
  !   kfn            :  Fermi wavenumber for neutrons
  !   kfp            :  Fermi wavenumber for protons
  !   kf0            :  Fermi wavenumber for total density
  !   eta            :  local asymmetry
  !   iso            :  for 1 (2) the function calculates the neutron (proton) 
  !                     gap by interpolation.
  ! Output:
  !   delta   : interpolated gap
  !-----------------------------------------------------------------------------
  procedure(delta_abstract), pointer :: delta_function
  real(KIND=dp), intent(in)          :: kfn(mv), kfp(mv), kf0(mv), eta(mv)
  integer, intent(in)                :: iso 
  real(KIND=dp)                      :: Delta(mv), fac(mv)
  real(KIND=dp)                      :: deltaq(mv), deltaN0(mv), deltaS0(mv)
  integer                            :: i

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Calculate gaps that occur in the interpolation formula for both cases
  deltaN0 = delta_function(kf0, 1) ! pairing gap in pure neutron matter at k_F0
  deltaS0 = delta_function(kf0, 3) !                symmetric    matter at k_F0
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  do i=1,mv
    ! Here we silently assume that (if delta_N0 is very small =>  delta_S0
    ! is even smaller) to avoid numerical issues. This is certainly not 
    ! universally true, but (almost everywhere) valid for the Cao gaps at least.
    if(abs(deltaN0(i)) .lt. 1e-10 ) then
      fac(i) = 0.0d0
    else
      fac(i) = deltaS0(i)/ deltaN0(i)
    endif
  enddo
 
  select case(iso)
  case(1)
    ! neutrons
    deltaq = delta_function(kfn, 1)
    Delta = deltaq * fac **(1-eta)
  case(2)
    ! protons
    deltaq = delta_function(kfp, 2)
    Delta = deltaq * fac **(1+eta)
  case DEFAULT
    call stp('Unrecognised input value for iso in weak_coupling_interpolation.')
  end select

 end function weak_coupling_interpolation
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
  ! Calculation of gaps in infinite neutron, proton and symmetric matter through
  ! functions fitted to data from 
  !     L. G. Cao et al., Phys. Rev C., 74(6):06430 (2006).
  !
  ! All expressions coded were taken from the Brussels axial HFB code.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Input:
  !   kf : Fermi wavelength
  !   iso: 1 => neutrons
  !        2 => protons
  !        3 => symmetric matter
  !
  ! Output:
  !  delta : gap for the given iso option at kf
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in) :: kf(mv)
  integer, intent(in)       :: iso
  real(KIND=dp)             :: delta(mv), xkfint, xkfmax
  integer                   :: i
  
  select case(iso)
  case(1,2)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Symmetric matter
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

