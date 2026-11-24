!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
!    Copyright (C) 2026 W. Ryssens and M. Bender
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Affero General Public License as published
!    by the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!===============================================================================
module folding
 !==============================================================================
 ! Module that is useful for folding various functions on the mesh. 
 !==============================================================================
 
 use geninfo

 implicit none

 !------------------------------------------------------------------------------
 ! Form factor matrices (their mathematematical form depend on the form factor parameterization chosen)
 !, to be used when folding of the nucleon densities to 
 ! obtain the charge densities are required. 
 ! Notes:
 !  - Neutrons and protons get separate definitions as folding is not 
 !    necessarily active for both species at once.
 !  - The fourth index ranges over the different gaussian form factors 
 !    to be applied. 
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable :: form_fact_x_neutron(:,:,:,:)
 real(KIND=dp), allocatable :: form_fact_y_neutron(:,:,:,:)
 real(KIND=dp), allocatable :: form_fact_z_neutron(:,:,:,:)

 real(KIND=dp), allocatable :: form_fact_x_proton(:,:,:,:)
 real(KIND=dp), allocatable :: form_fact_y_proton(:,:,:,:)
 real(KIND=dp), allocatable :: form_fact_z_proton(:,:,:,:)

 interface fold_form_factor
      module procedure fold_form_factor_real
      module procedure fold_form_factor_complex
 end interface

 interface fold_form_factor_reverse
      module procedure fold_form_factor_reverse_real
      module procedure fold_form_factor_reverse_complex
 end interface


 ! Friedrich-Walcher 2003 parameters 
 ! a11, a21, ab, Qb and sb are in GeV-based units 
 integer, parameter :: fw_proton = 1
 integer, parameter :: fw_neutron = 2

 real(kind=dp), parameter :: fw_p_a10 =  1.041_dp
 real(kind=dp), parameter :: fw_p_a11 =  0.765_dp
 real(kind=dp), parameter :: fw_p_a20 = -0.041_dp
 real(kind=dp), parameter :: fw_p_a21 =  6.2_dp
 real(kind=dp), parameter :: fw_p_ab  = -0.23_dp
 real(kind=dp), parameter :: fw_p_Qb  =  0.07_dp
 real(kind=dp), parameter :: fw_p_sb  =  0.27_dp

 real(kind=dp), parameter :: fw_n_a10 =  1.04_dp
 real(kind=dp), parameter :: fw_n_a11 =  1.73_dp
 real(kind=dp), parameter :: fw_n_a20 = -1.04_dp
 real(kind=dp), parameter :: fw_n_a21 =  1.54_dp
 real(kind=dp), parameter :: fw_n_ab  =  0.23_dp
 real(kind=dp), parameter :: fw_n_Qb  =  0.29_dp
 real(kind=dp), parameter :: fw_n_sb  =  0.20_dp


contains


 pure function gaussian(r1,r2, r0) result(G)
  !-----------------------------------------------------------------------------
  ! Returns the 1D Gaussian function given by
  ! 
  !     G(|r1 - r2|,n) = exp[-(r/r0)**2]
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in) :: r1, r2, r0
  real(KIND=dp) :: dr, G
  
  dr = abs(r1-r2)/(r0)
  G = 1.0/(r0 * sqrt(pi)) * exp(-dr**2)
  return
 end function gaussian

subroutine gauss_1D(G, mesh, m, r0, p)
    !---------------------------------------------------------------------------
    ! Function that constructs a matrix to fold in 1-D with a Gaussian of
    ! parameter r0 and symmetry sign p. 
    !
    ! Input:
    !   mesh : mesh coordinates along this direction
    !   m    : number of mesh points along this direction
    !   p    : symmetry sign for reflection along this direction.
    !   r0   : width of the Gaussian
    !
    ! Output:
    !   G  :  Constructed folding matrix
    !
    !          G(i,j) =     (r0 * sqrt(pi))**(-1) * exp[-(|+r_i - r_j|/r0)**2]
    !                 + p * (r0 * sqrt(pi))**(-1) * exp[-(|-r_i - r_j|/r0)**2]
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! A few notes:
    !
    !  *) After construction, we normalize the folding matrix to offset 
    !     numerical errors due to the mesh discretisation. In this way, we 
    !     can guarantee that 
    !
    !         int dx int dx' G(x,x') f(x') = sum_ij G(i,j) f(j) = int dx f(x)  
    !                                                           = sum_i  f(j)
    !     i.e. that we don't change integrals on the mesh.
    ! 
    !  *) If the relevant Cartesian axis is completely stored in the code, 
    !     meaning that there is no reflection symmetry, it is up to the user
    !     to pass in p = 0, removing trivially the contribution of any 
    !     reflection symmetry along the axis.
    !---------------------------------------------------------------------------

    real(KIND=dp)             :: G(m,m)
    integer, intent(in)       ::  m, p
    real(KIND=dp), intent(in) :: r0, mesh(m)
    
    integer :: i,j, ind
    
    ! Elements actually represented on the mesh
    do i=1,m
        do j=1,m          
            G(i,j) = Gaussian(mesh(i), mesh(j), r0)
        enddo
    enddo
    !---------------------------------------------------------------------------
    ! Elements to be gotten by symmetry.
    do i=1,m
        do j=1,m          
            G(i,j) = G(i,j) + p*Gaussian(-mesh(i), mesh(j), r0)
        enddo
    enddo
    !---------------------------------------------------------------------------
    !NS: for periodic boundary conditions add contrubution from 2 (symmetric)
    !neighboors. Should be enough for realistic box sizes due to rapid fall down
    !of the exponent.
#if(USE_Periodic==1) 
    do i=1,m
        do j=1,m          
            G(i,j) = G(i,j) + Gaussian(mesh(i)-(1+p)*m*dx, mesh(j), r0)        &
            &      +  Gaussian((1-2*p)*mesh(i)+(1+p)*m*dx, mesh(j), r0)
        enddo
    enddo
#endif
    !---------------------------------------------------------------------------
    ! Normalize, to avoid the numerical errors due to the mesh discretization.
    ! Technical note: we normalize all columns with the norm of ONE PARTICULAR
    !                 column, chosen "sufficiently far away" from the boundary
    !                 of the mesh. If we would normalize G for j = m, on the 
    !                 boundary, we would divide by too small a number, as the 
    !                 Gaussian should extend BEYOND the mesh. 
    !                 Naively, we could choose ind = 1 for this, 
    !                 but this is ON the boundary of the mesh when this axis
    !                 is not reduced through a conserved symmetry. For 
    !                 reasonable meshes and reasonable folding sizes, m/2+1
    !                 is several points away from either boundary of the mesh.
    ind = m/2 + 1 
    
    G(:,:) = G(:,:)/(sum(G(:,ind)*dx))
 end subroutine gauss_1D
 

 !Trying to define the galster and dipole form factors in r space


 !---------------------------------------------------------------------------
 ! Proton electric form factor: dipole spatial distribution
 !
 ! Input:
 !   r1 and r2 with r2-r1 : radial distance in fm
 !
 ! Output:
 !   rho_p : spatial density in fm**(-3) = 
 !---------------------------------------------------------------------------

pure function dipole(r1, r2, beta_GeV2) result(D)

  real(kind=dp), intent(in) :: r1, r2, beta_GeV2
  real(kind=dp) :: D
  real(kind=dp) :: r, m_D_fm

  real(kind=dp), parameter :: hbarc = 0.1973269804_dp

  ! Distance along the folding axis (in fm)
  r = abs(r1-r2)

  ! Convert the dipole parameter beta from GeV^2 to fm^(-2)
  ! The dipole density in 3D is given by:
  !
  ! rho_D(r) = m_D^3/(8*pi) * exp(-m_D*r)
  !
  ! with m_D = sqrt(beta)/hbarc
  m_D_fm = sqrt(beta_GeV2)/hbarc

  !---------------------------------------------------------------------------
  ! Change of coordinates and transverse integration
  !
  ! The original 3D dipole density depends on the radial distance:
  !
  !   r_3D = sqrt(x^2 + y^2 + z^2)
  !
  ! To construct a 1D folding kernel along the x-axis, we integrate the density
  ! over the two transverse coordinates y and z:
  !
  !   K_D(x) = integral dy dz rho_D(sqrt(x^2+y^2+z^2))
  !
  ! Using cylindrical coordinates in the transverse plane:
  !
  !   y = s*cos(phi), z = s*sin(phi)
  !   dy dz = s*ds*dphi
  !
  ! This gives:
  !
  !   K_D(x) = 2*pi * integral_0^infinity s ds
  !             * m_D^3/(8*pi)
  !             * exp(-m_D*sqrt(x^2+s^2))
  !
  ! After performing the integral analytically:
  !
  !   K_D(x) = m_D^2/4 * (1+m_D*|x|) * exp(-m_D*|x|)
  !
  ! This is the transverse marginal of the 3D dipole density. It is normalized
  ! to unity over the x coordinate and is the kernel used for 1D folding.
  !---------------------------------------------------------------------------

  D = m_D_fm**2 / 4.0_dp * (1.0_dp + m_D_fm*r) &
      & * exp(-m_D_fm*r)

end function dipole

subroutine dipole_1D(D, mesh, m, beta_GeV2, p)

    !---------------------------------------------------------------------------
    ! Constructs a 1D foldinD matrix usinD the dipole form factor kernel
    !
    ! Input:
    !   mesh       : mesh coordinates along this direction
    !   m          : number of mesh points
    !   beta_GeV2 : dipole parameter in GeV^2
    !   p          : symmetry sign for reflection along this direction
    !
    ! Output:
    !   G          : constructed folding matrix
    !
    ! The 1D kernel is the transverse marginal of the 3D dipole density:
    !
    ! K_D(x) = m_D^2/4 * (1 + m_D*|x|) * exp(-m_D*|x|)
    !
    ! where m_D = sqrt(beta)/hbarc.
    !---------------------------------------------------------------------------

    real(kind=dp)             :: D(m,m)
    integer, intent(in)       :: m, p
    real(kind=dp), intent(in) :: beta_GeV2, mesh(m)

    integer :: i, j, ind

    ! Elements actually represented on the mesh
    do i=1,m
        do j=1,m
            D(i,j) = dipole(mesh(i), mesh(j), beta_GeV2)
        enddo
    enddo

    !---------------------------------------------------------------------------
    ! Elements obtained by symmetry
    do i=1,m
        do j=1,m
            D(i,j) = D(i,j) + p*dipole(-mesh(i), mesh(j), beta_GeV2)
        enddo
    enddo

    !---------------------------------------------------------------------------
    ! Periodic boundary conditions
#if(USE_Periodic==1)
    do i=1,m
        do j=1,m
            D(i,j) = D(i,j) &
            & + dipole(mesh(i)-(1+p)*m*dx, mesh(j), beta_GeV2) &
            & + dipole((1-2*p)*mesh(i)+(1+p)*m*dx, mesh(j), beta_GeV2)
        enddo
    enddo
#endif

    !---------------------------------------------------------------------------
    ! Normalize to compensate for mesh discretization errors
    ind = m/2 + 1

    D(:,:) = D(:,:) / sum(D(:,ind)*dx)

end subroutine dipole_1D

 !---------------------------------------------------------------------------
 ! Neutron electric form factor: generalized Galster spatial distribution
 !
 ! Input:
 !   r : radial distance in fm
 !
 ! Output:
 !   rho_n : spatial density in fm**(-3)
 !---------------------------------------------------------------------------



pure function galster(r1, r2) result(G)

  real(kind=dp), intent(in) :: r1, r2
  real(kind=dp) :: G

  real(kind=dp) :: radius
  real(kind=dp) :: alpha, beta
  real(kind=dp) :: mass_alpha, mass_beta
  real(kind=dp) :: coeff_c, coeff_e

  real(kind=dp), parameter :: hbarc = 0.1973269804_dp
  real(kind=dp), parameter :: mn = 0.939565421_dp
  real(kind=dp), parameter :: lambda2 = 0.71_dp

  real(kind=dp), parameter :: A_g = 1.73_dp
  real(kind=dp), parameter :: B_g = 5.6_dp

  ! Distance along the folding axis (in fm)
  radius = abs(r1-r2)

  ! Momentum scales, expressed in fm^(-2)
  alpha = 4.0_dp*mn**2 / (B_g*hbarc**2)
  beta  = lambda2 / hbarc**2

  mass_alpha = sqrt(alpha)
  mass_beta  = sqrt(beta)

  ! Coefficients from the partial fraction decomposition
  coeff_c = -alpha*beta**2 / (beta-alpha)**2
  coeff_e = beta**3 / (beta-alpha)

  !---------------------------------------------------------------------------
  ! Change of coordinates and transverse integration
  !
  ! The original 3D Galster distribution depends on:
  !
  !   r_3D = sqrt(x^2 + y^2 + z^2)
  !
  ! To construct a 1D folding kernel along the x-axis, we integrate the density
  ! over the two transverse coordinates y and z:
  !
  !   K_n(x) = integral dy dz rho_n(sqrt(x^2+y^2+z^2))
  !
  ! Using cylindrical coordinates in the transverse plane:
  !
  !   y = s*cos(phi), z = s*sin(phi)
  !   dy dz = s*ds*dphi
  !
  ! This gives:
  !
  !   K_n(x) = 2*pi * integral_0^infinity s ds rho_n(sqrt(x^2+s^2))
  !
  ! The transverse integrals can be evaluated analytically:
  !
  !   integral_0^infinity s*exp(-m*sqrt(x^2+s^2))
  !                    /sqrt(x^2+s^2) ds = exp(-m*|x|)/m
  !
  !   integral_0^infinity s*exp(-m*sqrt(x^2+s^2)) ds
  !                  = exp(-m*|x|)*(|x|/m + 1/m^2)
  !
  ! The resulting 1D Galster kernel is:
  !
  !   K_n(x) = (A_g/B_g) * [
  !       C/(2*m_alpha) * exp(-m_alpha*|x|)
  !      -C/(2*m_beta)  * exp(-m_beta*|x|)
  !      +E/(4*m_beta^3) * (1+m_beta*|x|) * exp(-m_beta*|x|)
  !   ]
  !
  ! No division by radius is present in the final expression, therefore the
  ! numerical cancellation occurring in the original 3D expression at r=0
  ! is avoided.
  !
  !
  !
  ! Probleme pour la normalisation a 0 ...
  !---------------------------------------------------------------------------

  G = (A_g/B_g) * ( &
       coeff_c/(2.0_dp*mass_alpha) * exp(-mass_alpha*radius) &
       - coeff_c/(2.0_dp*mass_beta) * exp(-mass_beta*radius) &
       + coeff_e/(4.0_dp*mass_beta**3) &
         * (1.0_dp+mass_beta*radius) * exp(-mass_beta*radius) )


end function galster

subroutine galster_1D(G, mesh, m, p)

    real(kind=dp)             :: G(m,m)
    integer, intent(in)       :: m, p
    real(kind=dp), intent(in) :: mesh(m)

    integer :: i, j
    integer :: log_unit, ios, ind
    real(kind=dp) :: galster_integral

    ! Construct the kernel
    do i=1,m
        do j=1,m
            G(i,j) = galster(mesh(i), mesh(j))
        enddo
    enddo

    ! Symmetry contribution
    do i=1,m
        do j=1,m
            G(i,j) = G(i,j) &
                & + p*galster(-mesh(i), mesh(j))
        enddo
    enddo

#if(USE_Periodic==1)
    ! Periodic images
    do i=1,m
        do j=1,m
            G(i,j) = G(i,j) &
                & + galster(mesh(i)-(1+p)*m*dx, mesh(j)) &
                & + galster((1-2*p)*mesh(i)+(1+p)*m*dx, mesh(j))
        enddo
    enddo
#endif

    ! Numerical integral of the central column
    ind = m/2 + 1
    galster_integral = sum(G(:,ind))*dx

    ! Write diagnostic to log file
    open(newunit=log_unit, &
         file='/home/schroeder/Desktop/Codes/CoulombCorrections/galster_errors.log', &
         status='unknown', position='append', action='write', iostat=ios)

    if (ios == 0) then
        write(log_unit,'(A,I6,A,I6,A,ES16.8)') &
            'm = ', m, ', p = ', p, &
            ', Galster integral = ', galster_integral
        close(log_unit)
    endif

end subroutine galster_1D


!-----------------------------------------------------------------------------
 ! Exact 3D Breit-frame density.
 ! This is useful for validation and plotting; the Cartesian folding matrices
 ! use fw(), the transverse marginal, rather than evaluating fw at |x|.
 ! Output units: fm^-3.
 !-----------------------------------------------------------------------------
 pure function fw(r_input, nucleon) result(rho)
   real(kind=dp), intent(in) :: r_input
   integer, intent(in) :: nucleon
   real(kind=dp) :: rho

   real(kind=dp), parameter :: hbarc_fw = 0.1973269804_dp
   real(kind=dp), parameter :: pi = 3.1415926535897932384626433832795_dp
   real(kind=dp) :: a10, a11, a20, a21, ab, Qb, sb
   real(kind=dp) :: m1, m2, a, qb_fm, s, r, r_safe
   real(kind=dp) :: dip, Rterm, Tterm, bump, bump0, R0

   select case(nucleon)
   case(fw_proton)
     a10 = fw_p_a10
     a11 = fw_p_a11
     a20 = fw_p_a20
     a21 = fw_p_a21
     ab  = fw_p_ab
     Qb  = fw_p_Qb
     sb  = fw_p_sb
   case(fw_neutron)
     a10 = fw_n_a10
     a11 = fw_n_a11
     a20 = fw_n_a20
     a21 = fw_n_a21
     ab  = fw_n_ab
     Qb  = fw_n_Qb
     sb  = fw_n_sb
   case default
     rho = 0.0_dp
     return
   end select

   r = abs(r_input)

  
   m1 = a11/hbarc_fw**2
   m2 = a21/hbarc_fw**2
   a  = ab*hbarc_fw**2
   qb_fm = Qb/hbarc_fw
   s  = sb/hbarc_fw

   dip = sqrt(pi/8.0_dp) * ( &
       a10*m1**1.5_dp*exp(-sqrt(m1)*r) &
     + a20*m2**1.5_dp*exp(-sqrt(m2)*r) )

   Rterm = qb_fm**3 + 3.0_dp*qb_fm*s**2*(1.0_dp-s**2*r**2)
   Tterm = 3.0_dp*qb_fm**2 + 3.0_dp*s**2 - s**4*r**2

   if (r < 1.0e-8_dp) then
     R0 = qb_fm**3 + 3.0_dp*qb_fm*s**2
     bump0 = 2.0_dp*a*s * &
       (R0*qb_fm + s**2*(3.0_dp*qb_fm**2+3.0_dp*s**2))
     bump = bump0
   else
     r_safe = r
     bump = (2.0_dp*a*s/r_safe) &
       * exp(-0.5_dp*s**2*r_safe**2) &
       * ( Rterm*sin(qb_fm*r_safe) &
           + s**2*r_safe*Tterm*cos(qb_fm*r_safe) )
   endif

   rho = (dip+bump)/(2.0_dp*pi)**1.5_dp
 end function fw


 !-----------------------------------------------------------------------------
 ! Build the FW folding matrix. 
 !-----------------------------------------------------------------------------
subroutine fw_1D(Fmat, mesh, m, p, nucleon)

    integer, intent(in) :: m, p, nucleon
    real(kind=dp), intent(in) :: mesh(m)
    real(kind=dp), intent(out) :: Fmat(m,m)

    integer :: i, j, ind
    real(kind=dp) :: fw_integral

    ! Construct the kernel matrix
    do i=1,m
        do j=1,m
            Fmat(i,j) = fw(mesh(i), mesh(j), nucleon)
        enddo
    enddo

    ! Reflection symmetry
    do i=1,m
        do j=1,m
            Fmat(i,j) = Fmat(i,j) &
                + p*fw(-mesh(i), mesh(j), nucleon)
        enddo
    enddo

#if(USE_Periodic==1)
    ! Periodic images: same convention as dipole_1D
    do i=1,m
        do j=1,m
            Fmat(i,j) = Fmat(i,j) &
                + fw(mesh(i)-(1+p)*m*dx, mesh(j), nucleon) &
                + fw((1-2*p)*mesh(i)+(1+p)*m*dx, mesh(j), nucleon)
        enddo
    enddo
#endif

    ! Numerical integral of the central column
    ind = m/2 + 1
    fw_integral = sum(Fmat(:,ind))*dx

    ! Proton: compensate for discretization errors
    if (nucleon == fw_proton) then

        if (abs(fw_integral) < 1.0e-12_dp) then
            error stop "FW proton: central-column integral is zero"

        endif

        Fmat(:,:) = Fmat(:,:) / fw_integral

    endif

    ! Neutron: do NOT normalize.
    ! Its physical integral is zero, so division is not meaningful.

end subroutine fw_1D


 pure function gamma_function(r1,r2, r0,n,L) result(G)
  !-----------------------------------------------------------------------------
  ! Returns the 1D Gamma function given by
  ! 
  !     G(|r1 - r2|) = exp[-(pi*r0/L*(n+1/2))**2]*cos(2*pi/L*(n+1/2)*(r1-r2))
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in) :: r1, r2, r0, L
  integer, intent(in)       :: n
  real(KIND=dp)             :: exponent, cos_arg, G
  
  exponent = pi * r0 / L * (n + 1.0_dp/2.0_dp)
  cos_arg  = 2.0_dp * pi / L * (n + 1.0_dp/2.0_dp) * (r1 - r2)
  G = exp(-exponent**2) * cos(cos_arg)
  return
 end function gamma_function
 
 










 subroutine gamma_1D(G, mesh, m, r0, p)
    !---------------------------------------------------------------------------
    ! Function that constructs a matrix to fold in 1-D with the Gamma matrix
    ! coming from the exact interpolation of the folding of the density matrix
    ! and a Gaussian
    !
    ! Input:
    !   mesh : mesh coordinates along this direction
    !   m    : number of mesh points along this direction
    !   p    : symmetry sign for reflection along this direction.
    !   r0   : width of the Gaussian
    !
    ! Output:
    !   G  :  Constructed folding matrix
    !
    !          G(i,j) =     sum_n gamma_n( r_i - r_j)
    !                 + p * sum_n gamma_n(-r_i - r_j)
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! A few notes:
    !
    !  *) Contrary to gauss_1D, this function must NOT be normalized.
    ! 
    !  *) If the relevant Cartesian axis is completely stored in the code, 
    !     meaning that there is no reflection symmetry, it is up to the user
    !     to pass in p = 0, removing trivially the contribution of any 
    !     reflection symmetry along the axis.
    !---------------------------------------------------------------------------

    real(KIND=dp)             :: G(m,m), L 
    integer, intent(in)       ::  m, p
    real(KIND=dp), intent(in) :: r0, mesh(m)
    integer :: i,j, N_total, n 

    ! Intializing to 0 the G array (Neccesary because of the sum over n)
    G = 0.0_dp

    ! The summation over N in the Gamma matrices goes over n \in  [0,N-1], with
    ! N being the total number of sample points of HALF of a dimension of the cell
    if(abs(p).eq.1) then
        N_total = m
    else
        N_total = m/2
    endif 
    
    !Getting the total length of the box
    L = 2*N_total*dx

    ! Elements actually represented on the mesh
    do i=1,m
        do j=1,m          
            do n=0,N_total-1
                G(i,j) = G(i,j)+gamma_function(mesh(i), mesh(j), r0, n,L)*2_dp/L
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    ! Elements to be gotten by symmetry.
    do i=1,m
        do j=1,m          
            do n=0,N_total-1
                G(i,j) = G(i,j) + p*gamma_function(-mesh(i), mesh(j), r0,n,L)*2_dp/L
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    !NS: for periodic boundary conditions add contribution from 2 (symmetric)
    !neighboors. Should be enough for realistic box sizes due to rapid fall down
    !of the exponent.
#if(USE_Periodic==1) 
    do i=1,m
        do j=1,m          
            do n=0, N_total-1
                G(i,j) = G(i,j) + gamma_function(mesh(i)-(1+p)*m*dx, mesh(j), r0,n,L)*2_dp/L    &
                &      +  gamma_function((1-2*p)*mesh(i)+(1+p)*m*dx, mesh(j), r0,n,L)*2_dp/L
            enddo
        enddo
    enddo
#endif
 end subroutine gamma_1D

 subroutine fold_form_factor_real(f, folded, mx, my, mz, sx, sy, sz) 
  !--------------------------------------------------------------------------
  ! Perform the folding of a function f with neutron and proton components, 
  ! with the relevant form factors for neutrons and protons; summing 
  ! the contributions of protons and neutrons in the final result.
  !
  ! Important: 
  ! - this function assumes that the gaussian folding matrices
  !   have correctly been initialized beforehand, using the
  !   construct_folding_matrices subroutine.
  ! - if no folding is needed, the output folded is not allocated
  ! 
  ! Input:
  !  f       : the function to be folded, with a neutron and proton component
  ! mx/my/mz : the number of mesh points along the X/Y/Z direction
  ! sx/sy/sz : the symmetry sign for folding along the X/Y/Z direction
  !
  ! Output:
  ! Folded : the folded function
  !          with summed contributions of neutrons and protons,
  !---------------------------------------------------------------------------
  real(KIND=dp), intent(in)               :: f(:,:,:,:)
  integer, intent(in)                     :: mx,my,mz
  integer, intent(in)                     :: sx,sy,sz
  real(KIND=dp), intent(out), allocatable :: folded(:,:,:)
  real(KIND=dp)                           :: fn(mx,my,mz), fp(mx,my,mz)

  if(allocated(form_fact_x_neutron) .or. allocated(form_fact_x_proton)) then
    allocate(folded(mx,my,mz))
    folded = 0.0d0
  else
    return ! early return if no folding is needed
  endif 

  fn = fold_form_factor_species(f(:,:,:,1),      &
  &                             form_fact_x_neutron, & 
  &                             form_fact_y_neutron, &
  &                             form_fact_z_neutron, &
  &                             mx, my, mz, sx, sy, sz)

  fp = fold_form_factor_species(f(:,:,:,2),     &
  &                             form_fact_x_proton, & 
  &                             form_fact_y_proton, &
  &                             form_fact_z_proton, &
  &                             mx, my, mz, sx, sy, sz)

  ! Sum the contributions
  folded = fn + fp
 end subroutine fold_form_factor_real

 subroutine fold_form_factor_complex(f, folded, mx, my, mz, sx, sy, sz) 
  !--------------------------------------------------------------------------
  ! Perform the folding of a function f with neutron and proton components, 
  ! with the relevant form factors for neutrons and protons; summing 
  ! the contributions of protons and neutrons in the final result.
  !
  ! Important: 
  ! - this function assumes that the gaussian folding matrices
  !   have correctly been initialized beforehand, using the
  !   construct_folding_matrices subroutine.
  ! - if no folding is needed, the output folded is not allocated
  ! 
  ! Input:
  !  f       : the function to be folded, with a neutron and proton component
  ! mx/my/mz : the number of mesh points along the X/Y/Z direction
  ! sx/sy/sz : the symmetry sign for folding along the X/Y/Z direction
  !
  ! Output:
  ! Folded : the folded function
  !          with summed contributions of neutrons and protons,
  !---------------------------------------------------------------------------
  complex(KIND=dp), intent(in)               :: f(:,:,:,:)
  integer, intent(in)                        :: mx,my,mz
  integer, intent(in)                        :: sx,sy,sz
  complex(KIND=dp), intent(out), allocatable :: folded(:,:,:)

  real(KIND=dp)                              :: f_re(mx,my,mz,2)
  real(KIND=dp)                              :: f_im(mx,my,mz,2)
  real(KIND=dp), allocatable                 :: folded_re(:,:,:)
  real(KIND=dp), allocatable                 :: folded_im(:,:,:)

  if(allocated(form_fact_x_neutron) .or. allocated(form_fact_x_proton)) then
    allocate(folded(mx,my,mz))
    folded = 0.0d0
  else
    return ! early return if no folding is needed
  endif 

  ! Split the complex function into real and imaginary parts
  f_re = DBLE( f)
  f_im = AIMAG(f)

  ! Fold each part separately
  call fold_form_factor_real(f_re, folded_re, mx, my, mz, sx, sy, sz)
  call fold_form_factor_real(f_im, folded_im, mx, my, mz, sx, sy, sz)

  ! recombine complex numbers
  folded = CMPLX(folded_re, folded_im, dp)

 end subroutine fold_form_factor_complex

 subroutine fold_form_factor_reverse_real(f, folded, mx, my, mz, sx, sy, sz) 
  !--------------------------------------------------------------------------
  ! Perform the folding of an (isoscalar) function f with the relevant form 
  ! factors for neutrons and protons, returning both components separately.
  !
  ! Important: 
  ! - this function assumes that the gaussian folding matrices
  !   have correctly been initialized beforehand, using the
  !   construct_folding_matrices subroutine.
  ! - if no folding is needed, the output folded is not allocated
  ! 
  ! Input:
  !  f       : the function to be folded
  ! mx/my/mz : the number of mesh points along the X/Y/Z direction
  ! sx/sy/sz : the symmetry sign for folding along the X/Y/Z direction
  !
  ! Output:
  ! Folded : the folded function, differentiating between protons and neutrons
  !---------------------------------------------------------------------------
  real(KIND=dp), intent(in)         :: f(:,:,:)
  integer, intent(in)               :: mx,my,mz
  integer, intent(in)               :: sx,sy,sz
  real(KIND=dp), intent(out), allocatable :: folded(:,:,:,:)

  if(allocated(form_fact_x_neutron) .or. allocated(form_fact_x_proton)) then
    allocate(folded(mx,my,mz,2))
    folded = 0.0d0
  else
    return ! early return if no folding is needed
  endif 
  

  folded(:,:,:,1) = fold_form_factor_species(f,  &
  &                             form_fact_x_neutron, & 
  &                             form_fact_y_neutron, &
  &                             form_fact_z_neutron, &
  &                             mx, my, mz, sx, sy, sz)

  folded(:,:,:,2) = fold_form_factor_species(f, &
  &                             form_fact_x_proton, & 
  &                             form_fact_y_proton, &
  &                             form_fact_z_proton, &
  &                             mx, my, mz, sx, sy, sz)

 end subroutine fold_form_factor_reverse_real
 
 subroutine fold_form_factor_reverse_complex(f, folded, mx, my, mz, sx, sy, sz) 
  !--------------------------------------------------------------------------
  ! Perform the folding of an (isoscalar) function f with the relevant form 
  ! factors for neutrons and protons, returning both components separately.
  !
  ! Important: 
  ! - this function assumes that the gaussian folding matrices
  !   have correctly been initialized beforehand, using the
  !   construct_folding_matrices subroutine.
  ! - if no folding is needed, the output folded is not allocated
  ! 
  ! Input:
  !  f       : the function to be folded
  ! mx/my/mz : the number of mesh points along the X/Y/Z direction
  ! sx/sy/sz : the symmetry sign for folding along the X/Y/Z direction
  !
  ! Output:
  ! Folded : the folded function, differentiating between protons and neutrons
  !---------------------------------------------------------------------------
  complex(KIND=dp), intent(in)               :: f(:,:,:)
  integer, intent(in)                        :: mx,my,mz
  integer, intent(in)                        :: sx,sy,sz
  complex(KIND=dp), intent(out), allocatable :: folded(:,:,:,:)

  real(KIND=dp)                              :: f_re(mx,my,mz)
  real(KIND=dp)                              :: f_im(mx,my,mz)
  real(KIND=dp), allocatable                 :: folded_re(:,:,:,:)
  real(KIND=dp), allocatable                 :: folded_im(:,:,:,:)

  if(allocated(form_fact_x_neutron) .or. allocated(form_fact_x_proton)) then
    allocate(folded(mx,my,mz,2))
    folded = 0.0d0
  else
    return ! early return if no folding is needed
  endif 

  ! Separate the complex function into real and imaginary parts
  f_re = DBLE( f)
  f_im = AIMAG(f)

  ! ... and perform the folding for each part individually
  call fold_form_factor_reverse_real(f_re, folded_re, mx, my, mz, sx, sy, sz)
  call fold_form_factor_reverse_real(f_im, folded_im, mx, my, mz, sx, sy, sz)

  ! recombine complex numbers
  folded = CMPLX(folded_re, folded_im, dp)

 end subroutine fold_form_factor_reverse_complex

 function fold_form_factor_species(f, FFx,FFy,FFz, mx, my, mz, sx, sy, sz) result(Folded)
  !----------------------------------------------------------------------------------
  ! Fold a function f with possibly multiple form Gaussian form factors.
  ! 
  ! Input:
  !   f       : the function to be folded, defined on a (mx,my,mz) mesh
  !   Gx/y/z  : the folding matrices, defined on a (mx,my,mz) mesh
  !   mx/my/mz: the number of mesh points along the X/Y/Z direction
  !   sx/sy/sz: the symmetry sign for folding along the X/Y/Z direction
  !
  ! Output:
  !   folded : the folded function, defined on a (mx,my,mz) mesh
  !-----------------------------------------------------------------------------------
  real(KIND=dp), intent(in)              :: f(mx,my,mz)
  real(KIND=dp), intent(in), allocatable :: FFx(:,:,:,:), FFy(:,:,:,:), FFz(:,:,:,:)
  integer, intent(in)                    :: mx,my,mz
  integer, intent(in)                    :: sx,sy,sz
  real(KIND=dp)                          :: folded(mx,my,mz)

   folded = 0.0d0 
   if(allocated(FFx)) then
    ! Fold with a Gaussian with positive sign 
    folded =   folded + fold_one_form_factor(f, &
    &                                        FFx(:,:,:,1), &
    &                                        FFy(:,:,:,1), &
    &                                        FFz(:,:,:,1), & 
    &                                        mx, my, mz, sx, sy, sz)
    if(size(FFx,4) == 2) then 
      ! ... and with a Gaussian with negative sign if needed!
      folded = folded - fold_one_form_factor(f, &
      &                                      FFx(:,:,:,2), &
      &                                      FFy(:,:,:,2), &
      &                                      FFz(:,:,:,2), & 
      &                                      mx, my, mz, sx, sy, sz)
    endif
  endif

 end function fold_form_factor_species

 function fold_one_form_factor(f,FFx,FFy,FFz, mx, my, mz, sx, sy, sz) result(Folded)
  !-----------------------------------------------------------------------------
  !(Emilien : up to date)
  ! Returns the folded function 
  !
  !
  !  Folded(x2,y2,z2) = int dx1 dy1 dz1 f(x1,y1,z1) 
  !                                Gx(|x1-x2|)Gx(|y1-y2|)Gx(|z1-z2|)
  !
  ! Where the folding functions Gx/Gy/Gz should be constructed beforehand and
  ! passed in. 
  !
  ! Input:
  !   f     : the function to be folded, defined on a (mx,my,mz) mesh
  !   FFx/y/z : the folding matrices, defined on a (mx,my,mz) mesh
  !   mx/my/mz : the number of mesh points along the X/Y/Z direction
  !   sx/sy/sz : the symmetry sign for folding along the X/Y/Z direction
  ! 
  ! Output:
  !   Folded : the folded function, defined on a (mx,my,mz) mesh
  !-----------------------------------------------------------------------------   
  real(KIND=dp), intent(in)         :: f(:,:,:)
  real(KIND=dp), intent(in)         :: FFx(:,:,:), FFy(:,:,:), FFz(:,:,:)
  integer, intent(in)               :: mx,my,mz
  integer, intent(in)               :: sx,sy,sz
  real(KIND=dp)                     :: folded(mx,my,mz)
  integer                           :: i,j,k, index
  !-----------------------------------------------------------------------------
  ! First convolute/fold along the X-direction
  index = 1 + (sx + 1)/2 ! index = 1 for sx = -1/0, index = 2 for sx = +1
  do k=1,mz
    do j=1,my
      folded(:,j,k) = matmul(FFx(:,:,index), f(1:mx,j,k))*dx
    enddo
  enddo
  ! Then the Y-direction
  index = 1 + (sy + 1)/2 ! index = 1 for sy = -1/0, index = 2 for sy = +1
  do k=1,mz
    do i=1,mx
      folded(i,:,k) = matmul(FFy(:,:,index), folded(i,1:my,k))*dx
    enddo
  enddo
  ! Then the Z-direction
  index = 1 + (sz + 1)/2 ! index = 1 for sz = -1/0, index = 2 for sz = +1
  do j=1,my
    do i=1,mx
      folded(i,j,:) = matmul(FFz(:,:,index), folded(i,j,1:mz))*dx
    enddo
  enddo
  !-----------------------------------------------------------------------------
 end function fold_one_form_factor

  subroutine construct_folding_matrices(proton_size, neutron_size, hocomform, & 
 &                             exact_coulomb_folding, nucleon_form_factor, hbm,FFxn,FFyn,FFzn, FFxp, FFyp, FFzp)
    !---------------------------------------------------------------------------
    ! Construct the matrices for Gaussian folding.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input :
    !  protonsize, neutronsize : sizes of the Gaussians for folding
    !  hoconform  : whether to apply the harmonic-oscillator correction
    !  exact_coulomb_folding: whether to apply the integral interpolation procedure
    !  hbm        : hbar^2/m for use in the harmonic-oscillator correction
    !
    ! Output:
    !   FFx/y/zn/p : Form factors for folding the density
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: neutron_size(2), proton_size(2), hbm(2)
    logical, intent(in)        :: hocomform, exact_coulomb_folding
    real(KIND=dp), allocatable, intent(out) :: FFxn(:,:,:,:), FFyn(:,:,:,:), FFzn(:,:,:,:)
    real(KIND=dp), allocatable, intent(out) :: FFxp(:,:,:,:), FFyp(:,:,:,:), FFzp(:,:,:,:)
    real(KIND=dp)              :: rplus_n, rplus_p, rmin_n, rmin_p
    real(KIND=dp)              :: hbom, mhb, B
    integer                    :: n_gauss_n, n_gauss_p, sign, index
    real(kind=dp), parameter :: beta_dipole = 0.71_dp
    integer, intent(in) :: nucleon_form_factor
!For now, I will hard code the parameters of my new parameterization in here, before I can discuss with Wouter the good strategy
!maybe the most elegant solution would be to use instantiation insted of  cases, so the variables are managed in a better way

    ! The determination of folding parameters from the parameterization input 
    ! for neutrons and protons is not the same; see documentation.
    rplus_n = sqrt(neutron_size(1))
    rmin_n  = sqrt(neutron_size(2))

    rplus_p = proton_size(1) * sqrt(2.0d0/3.0d0)
    rmin_p  = proton_size(2) * sqrt(2.0d0/3.0d0)


  





    ! Memory management
  
    if(allocated(FFxn)) deallocate(FFxn)
    if(allocated(FFyn)) deallocate(FFyn)
    if(allocated(FFzn)) deallocate(FFzn)
    if(allocated(FFxp)) deallocate(FFxp)
    if(allocated(FFyp)) deallocate(FFyp)
    if(allocated(FFzp)) deallocate(FFzp)

    ! Allocate the folding matrices
    n_gauss_n = 0
    if(rplus_n .gt. 0.0d0) n_gauss_n = 1
    if(rmin_n  .gt. 0.0d0) n_gauss_n = 2

    n_gauss_p = 0
    if(rplus_p .gt. 0.0d0) n_gauss_p = 1
    if(rmin_p .gt. 0.0d0)  n_gauss_p = 2


    select case(nucleon_form_factor)

    case(0)
    ! correspond to the original code, where the form factor is an addition of gaussian    

    if(n_gauss_n .ne. 0) then
      allocate(FFxn(nx,nx,3,n_gauss_n), FFyn(ny,ny,3,n_gauss_n), FFzn(nz,nz,3,n_gauss_n))
      FFxn = 0.0_dp
      FFyn = 0.0_dp
      FFzn = 0.0_dp
    endif
    if(n_gauss_p .ne. 0) then 
      allocate(FFxp(nx,nx,3,n_gauss_p), FFyp(ny,ny,3,n_gauss_p), FFzp(nz,nz,3,n_gauss_p))
      FFxp = 0.0_dp
      FFyp = 0.0_dp
      FFzp = 0.0_dp
    endif
    ! early return if possible
    if(n_gauss_n .eq. 0 .and. n_gauss_p .eq. 0) return
    !---------------------------------------------------------------------------
    ! Harmonic-oscillator correction
    if(hocomform) then
        ! hbar x omega
        hbom  = 41.0d0 * (neutrons + protons)**(-1.0d0/3.0d0)
        ! 2m/hbar^2
        mhb = 2.0d0/(1.0d0/hbm(1)+1.0d0/hbm(2))
        ! B^{-1} = hbar * omega/m * A = 1/2 * A * hbar omega * 2m/hbar^2
        B = sqrt( 1.0d0/( 0.5d0 * hbom/mhb  * (neutrons + protons)))

        if(rplus_n.ne.0.0d0) then
            rplus_n = sqrt(rplus_n**2 - B**2)
        endif
        if(rmin_n.ne.0.0d0) then
            rmin_n   = sqrt(rmin_n**2 - B**2)
        endif

        if(rplus_p.ne.0.0d0) then
            rplus_p = sqrt(rplus_p**2 - B**2)
        endif
        if(rmin_p.ne.0.0d0) then
            rmin_p   = sqrt(rmin_p**2 - B**2)
        endif
    endif

    if (exact_coulomb_folding) then
    	do sign = -1,+1, 2
    	  index = 1 + (sign + 1)/2 ! index = 1 for sign = -1, index = 2 for sign = +1

    	  ! HACK: the multiplication with reduX/Y/Z ensures that the breaking of 
    	  !       a reflection symmetry automatically leads to Gaussian matrices 
    	  !       being constructed with no reflection symmetry.

    	  ! Neutrons
    	  if(rplus_n .ne. 0.0_dp) then
    	    call gamma_1D(FFxn(:,:,index,1), meshx, nx, rplus_n, sign * reduX)
    	    call gamma_1D(FFyn(:,:,index,1), meshy, ny, rplus_n, sign * reduY)
    	    call gamma_1D(FFzn(:,:,index,1), meshz, nz, rplus_n, sign * reduZ)
    	  endif
    	  if(rmin_n .ne. 0.0_dp) then
    	    call gamma_1D(FFxn(:,:,index,2), meshx, nx, rmin_n,  sign * reduX)
    	    call gamma_1D(FFyn(:,:,index,2), meshy, ny, rmin_n,  sign * reduY)
    	    call gamma_1D(FFzn(:,:,index,2), meshz, nz, rmin_n,  sign * reduZ)
    	  endif
    	  ! Protons
    	  if(rplus_p .ne. 0.0_dp) then
    	    call gamma_1D(FFxp(:,:,index,1), meshx, nx, rplus_p, sign * reduX)
    	    call gamma_1D(FFyp(:,:,index,1), meshy, ny, rplus_p, sign * reduY)
    	    call gamma_1D(FFzp(:,:,index,1), meshz, nz, rplus_p, sign * reduZ)
    	  endif
    	  if(rmin_p .ne. 0.0_dp) then
    	    call gamma_1D(FFxp(:,:,index,2), meshx, nx, rmin_p,  sign * reduX)
    	    call gamma_1D(FFyp(:,:,index,2), meshy, ny, rmin_p,  sign * reduY)
    	    call gamma_1D(FFzp(:,:,index,2), meshz, nz, rmin_p,  sign * reduZ)
    	  endif
    	enddo
    else
    	do sign = -1,+1, 2
    	  index = 1 + (sign + 1)/2 ! index = 1 for sign = -1, index = 2 for sign = +1

    	  ! HACK: the multiplication with reduX/Y/Z ensures that the breaking of 
    	  !       a reflection symmetry automatically leads to Gaussian matrices 
    	  !       being constructed with no reflection symmetry.

    	  ! Neutrons
    	  if(rplus_n .ne. 0.0_dp) then
    	    call gauss_1D(FFxn(:,:,index,1), meshx, nx, rplus_n, sign * reduX)
    	    call gauss_1D(FFyn(:,:,index,1), meshy, ny, rplus_n, sign * reduY)
    	    call gauss_1D(FFzn(:,:,index,1), meshz, nz, rplus_n, sign * reduZ)
    	  endif
    	  if(rmin_n .ne. 0.0_dp) then
    	    call gauss_1D(FFxn(:,:,index,2), meshx, nx, rmin_n,  sign * reduX)
    	    call gauss_1D(FFyn(:,:,index,2), meshy, ny, rmin_n,  sign * reduY)
    	    call gauss_1D(FFzn(:,:,index,2), meshz, nz, rmin_n,  sign * reduZ)
    	  endif
    	  ! Protons
    	  if(rplus_p .ne. 0.0_dp) then
    	    call gauss_1D(FFxp(:,:,index,1), meshx, nx, rplus_p, sign * reduX)
    	    call gauss_1D(FFyp(:,:,index,1), meshy, ny, rplus_p, sign * reduY)
    	    call gauss_1D(FFzp(:,:,index,1), meshz, nz, rplus_p, sign * reduZ)
    	  endif
    	  if(rmin_p .ne. 0.0_dp) then
    	    call gauss_1D(FFxp(:,:,index,2), meshx, nx, rmin_p,  sign * reduX)
    	    call gauss_1D(FFyp(:,:,index,2), meshy, ny, rmin_p,  sign * reduY)
    	    call gauss_1D(FFzp(:,:,index,2), meshz, nz, rmin_p,  sign * reduZ)
    	  endif
    	enddo
    endif


    case(1)




    ! early return if possible
    !I use the same logic to let the code know if a folding is necessary. Depending on the nucleon sizes
    !But I just need it to be non zero, probably there is an architectural choice to make here
    ! because its just a hack based on the previous code and is not relevant to the new form factors
    if(n_gauss_n .eq. 0 .and. n_gauss_p .eq. 0) return

    !---------------------------------------------------------------------------
    ! Harmonic-oscillator correction
    if(hocomform) then
        ! hbar x omega
        hbom  = 41.0d0 * (neutrons + protons)**(-1.0d0/3.0d0)
        ! 2m/hbar^2
        mhb = 2.0d0/(1.0d0/hbm(1)+1.0d0/hbm(2))
        ! B^{-1} = hbar * omega/m * A = 1/2 * A * hbar omega * 2m/hbar^2
        B = sqrt( 1.0d0/( 0.5d0 * hbom/mhb  * (neutrons + protons)))

        if(rplus_n.ne.0.0d0) then
            rplus_n = sqrt(rplus_n**2 - B**2)
        endif
        if(rmin_n.ne.0.0d0) then
            rmin_n   = sqrt(rmin_n**2 - B**2)
        endif

        if(rplus_p.ne.0.0d0) then
            rplus_p = sqrt(rplus_p**2 - B**2)
        endif
        if(rmin_p.ne.0.0d0) then
            rmin_p   = sqrt(rmin_p**2 - B**2)
        endif
    endif


    !here the last dimension is a relicate of the previous strategy two construct the neutron form factor
    ! I ll see if I can get rid of it after
    allocate(FFxn(nx,nx,3,1), FFyn(ny,ny,3,1), FFzn(nz,nz,3,1))
    FFxn = 0.0_dp
    FFyn = 0.0_dp
    FFzn = 0.0_dp

    allocate(FFxp(nx,nx,3,1), FFyp(ny,ny,3,1), FFzp(nz,nz,3,1))
    FFxp = 0.0_dp
    FFyp = 0.0_dp
    FFzp = 0.0_dp

    do sign = -1,+1, 2
      index = 1 + (sign + 1)/2 ! index = 1 for sign = -1, index = 2 for sign = +1

      ! HACK: the multiplication with reduX/Y/Z ensures that the breaking of 
      !       a reflection symmetry automatically leads to Gaussian matrices 
      !       being constructed with no reflection symmetry.

      ! Neutrons

      call galster_1D(FFxn(:,:,index,1), meshx, nx, sign * reduX)
      call galster_1D(FFyn(:,:,index,1), meshy, ny,  sign * reduY)
      call galster_1D(FFzn(:,:,index,1), meshz, nz,  sign * reduZ)

      ! Protons

      call dipole_1D(FFxp(:,:,index,1), meshx, nx, beta_dipole, sign * reduX)
      call dipole_1D(FFyp(:,:,index,1), meshy, ny, beta_dipole, sign * reduY)
      call dipole_1D(FFzp(:,:,index,1), meshz, nz, beta_dipole, sign * reduZ)

    enddo


    case(2)
    ! Friedrich-Walcher 2003 for both proton and neutron.
    ! One component is retained to match the existing matrix interface.

    allocate(FFxn(nx,nx,3,1), FFyn(ny,ny,3,1), FFzn(nz,nz,3,1))
    FFxn = 0.0_dp
    FFyn = 0.0_dp
    FFzn = 0.0_dp

    allocate(FFxp(nx,nx,3,1), FFyp(ny,ny,3,1), FFzp(nz,nz,3,1))
    FFxp = 0.0_dp
    FFyp = 0.0_dp
    FFzp = 0.0_dp

    do sign = -1,+1,2
      index = 1 + (sign + 1)/2

      ! Neutron: FW neutron coefficients.
      call fw_1D(FFxn(:,:,index,1), meshx, nx, sign*reduX, fw_neutron)
      call fw_1D(FFyn(:,:,index,1), meshy, ny, sign*reduY, fw_neutron)
      call fw_1D(FFzn(:,:,index,1), meshz, nz, sign*reduZ, fw_neutron)

      ! Proton: FW proton coefficients.
      call fw_1D(FFxp(:,:,index,1), meshx, nx, sign*reduX, fw_proton)
      call fw_1D(FFyp(:,:,index,1), meshy, ny, sign*reduY, fw_proton)
      call fw_1D(FFzp(:,:,index,1), meshz, nz, sign*reduZ, fw_proton)
    enddo

   case default
      error stop "Unknown nucleon form factor selection"

   end select

end subroutine construct_folding_matrices