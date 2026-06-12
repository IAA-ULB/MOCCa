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
module folding
 !==============================================================================
 ! Module that is useful for folding various functions on the mesh. 
 !==============================================================================
 
 use geninfo

 implicit none

 !------------------------------------------------------------------------------
 ! Gaussian matrices, to be used when folding of the nucleon densities to 
 ! obtain the charge densities are required. 
 ! Notes:
 !  - Neutrons and protons get separate definitions as folding is not 
 !    necessarily active for both species at once.
 !  - The fourth index ranges over the different gaussian form factors 
 !    to be applied. 
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable :: gauss_x_neutron(:,:,:,:)
 real(KIND=dp), allocatable :: gauss_y_neutron(:,:,:,:)
 real(KIND=dp), allocatable :: gauss_z_neutron(:,:,:,:)

 real(KIND=dp), allocatable :: gauss_x_proton(:,:,:,:)
 real(KIND=dp), allocatable :: gauss_y_proton(:,:,:,:)
 real(KIND=dp), allocatable :: gauss_z_proton(:,:,:,:)

 interface fold_form_factor
      module procedure fold_form_factor_real
      module procedure fold_form_factor_complex
 end interface

 interface fold_form_factor_reverse
      module procedure fold_form_factor_reverse_real
      module procedure fold_form_factor_reverse_complex
 end interface

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

  if(allocated(gauss_x_neutron) .or. allocated(gauss_x_proton)) then
    allocate(folded(mx,my,mz))
    folded = 0.0d0
  else
    return ! early return if no folding is needed
  endif 

  fn = fold_form_factor_species(f(:,:,:,1),      &
  &                             gauss_x_neutron, & 
  &                             gauss_y_neutron, &
  &                             gauss_z_neutron, &
  &                             mx, my, mz, sx, sy, sz)

  fp = fold_form_factor_species(f(:,:,:,2),     &
  &                             gauss_x_proton, & 
  &                             gauss_y_proton, &
  &                             gauss_z_proton, &
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

  if(allocated(gauss_x_neutron) .or. allocated(gauss_x_proton)) then
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

  if(allocated(gauss_x_neutron) .or. allocated(gauss_x_proton)) then
    allocate(folded(mx,my,mz,2))
    folded = 0.0d0
  else
    return ! early return if no folding is needed
  endif 
  

  folded(:,:,:,1) = fold_form_factor_species(f,  &
  &                             gauss_x_neutron, & 
  &                             gauss_y_neutron, &
  &                             gauss_z_neutron, &
  &                             mx, my, mz, sx, sy, sz)

  folded(:,:,:,2) = fold_form_factor_species(f, &
  &                             gauss_x_proton, & 
  &                             gauss_y_proton, &
  &                             gauss_z_proton, &
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

  if(allocated(gauss_x_neutron) .or. allocated(gauss_x_proton)) then
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

 function fold_form_factor_species(f, Gx,Gy,Gz, mx, my, mz, sx, sy, sz) result(Folded)
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
  real(KIND=dp), intent(in), allocatable :: Gx(:,:,:,:), Gy(:,:,:,:), Gz(:,:,:,:)
  integer, intent(in)                    :: mx,my,mz
  integer, intent(in)                    :: sx,sy,sz
  real(KIND=dp)                          :: folded(mx,my,mz)

   folded = 0.0d0 
   if(allocated(Gx)) then
    ! Fold with a Gaussian with positive sign 
    folded =   folded + fold_one_gaussian(f, &
    &                                        Gx(:,:,:,1), &
    &                                        Gy(:,:,:,1), &
    &                                        Gz(:,:,:,1), & 
    &                                        mx, my, mz, sx, sy, sz)
    if(size(Gx,4) == 2) then 
      ! ... and with a Gaussian with negative sign if needed!
      folded = folded - fold_one_gaussian(f, &
      &                                      Gx(:,:,:,2), &
      &                                      Gy(:,:,:,2), &
      &                                      Gz(:,:,:,2), & 
      &                                      mx, my, mz, sx, sy, sz)
    endif
  endif

 end function fold_form_factor_species

 function fold_one_gaussian(f,Gx,Gy,Gz, mx, my, mz, sx, sy, sz) result(Folded)
  !-----------------------------------------------------------------------------
  ! Returns the folded function 
  !
  !  Folded(x2,y2,z2) = int dx1 dy1 dz1 f(x1,y1,z1) 
  !                                Gx(|x1-x2|)Gx(|y1-y2|)Gx(|z1-z2|)
  !
  ! Where the folding functions Gx/Gy/Gz should be constructed beforehand and
  ! passed in. 
  !
  ! Input:
  !   f     : the function to be folded, defined on a (mx,my,mz) mesh
  !   Gx/y/z : the folding matrices, defined on a (mx,my,mz) mesh
  !   mx/my/mz : the number of mesh points along the X/Y/Z direction
  !   sx/sy/sz : the symmetry sign for folding along the X/Y/Z direction
  ! 
  ! Output:
  !   Folded : the folded function, defined on a (mx,my,mz) mesh
  !-----------------------------------------------------------------------------   
  real(KIND=dp), intent(in)         :: f(:,:,:)
  real(KIND=dp), intent(in)         :: Gx(:,:,:), Gy(:,:,:), Gz(:,:,:)
  integer, intent(in)               :: mx,my,mz
  integer, intent(in)               :: sx,sy,sz
  real(KIND=dp)                     :: folded(mx,my,mz)
  integer                           :: i,j,k, index
  !-----------------------------------------------------------------------------
  ! First convolute/fold along the X-direction
  index = 1 + (sx + 1)/2 ! index = 1 for sx = -1/0, index = 2 for sx = +1
  do k=1,mz
    do j=1,my
      folded(:,j,k) = matmul(Gx(:,:,index), f(1:mx,j,k))*dx
    enddo
  enddo
  ! Then the Y-direction
  index = 1 + (sy + 1)/2 ! index = 1 for sy = -1/0, index = 2 for sy = +1
  do k=1,mz
    do i=1,mx
      folded(i,:,k) = matmul(Gy(:,:,index), folded(i,1:my,k))*dx
    enddo
  enddo
  ! Then the Z-direction
  index = 1 + (sz + 1)/2 ! index = 1 for sz = -1/0, index = 2 for sz = +1
  do j=1,my
    do i=1,mx
      folded(i,j,:) = matmul(Gz(:,:,index), folded(i,j,1:mz))*dx
    enddo
  enddo
  !-----------------------------------------------------------------------------
 end function fold_one_gaussian

 subroutine construct_folding_matrices(proton_size, neutron_size, hocomform, & 
 &                             exact_coulomb_folding, hbm,Gxn,Gyn,Gzn, Gxp, Gyp, Gzp)
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
    !   Gx/y/zn/p : Gaussian factors for folding the density
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: neutron_size(2), proton_size(2), hbm(2)
    logical, intent(in)        :: hocomform, exact_coulomb_folding
    real(KIND=dp), allocatable, intent(out) :: Gxn(:,:,:,:), Gyn(:,:,:,:), Gzn(:,:,:,:)
    real(KIND=dp), allocatable, intent(out) :: Gxp(:,:,:,:), Gyp(:,:,:,:), Gzp(:,:,:,:)
    real(KIND=dp)              :: rplus_n, rplus_p, rmin_n, rmin_p
    real(KIND=dp)              :: hbom, mhb, B
    integer                    :: n_gauss_n, n_gauss_p, sign, index



    ! The determination of folding parameters from the parameterization input 
    ! for neutrons and protons is not the same; see documentation.
    rplus_n = sqrt(neutron_size(1))
    rmin_n  = sqrt(neutron_size(2))

    rplus_p = proton_size(1) * sqrt(2.0/3.0)
    rmin_p  = proton_size(2) * sqrt(2.0/3.0)

    ! Memory management 
    if(allocated(Gxn)) deallocate(Gxn)
    if(allocated(Gyn)) deallocate(Gyn)
    if(allocated(Gzn)) deallocate(Gzn)
    if(allocated(Gxp)) deallocate(Gxp)
    if(allocated(Gyp)) deallocate(Gyp)
    if(allocated(Gzp)) deallocate(Gzp)

    ! Allocate the folding matrices
    n_gauss_n = 0
    if(rplus_n .gt. 0.0d0) n_gauss_n = 1
    if(rmin_n  .gt. 0.0d0) n_gauss_n = 2

    n_gauss_p = 0
    if(rplus_p .gt. 0.0d0) n_gauss_p = 1
    if(rmin_p .gt. 0.0d0)  n_gauss_p = 2

    if(n_gauss_n .ne. 0) then
      allocate(Gxn(nx,nx,3,n_gauss_n), Gyn(ny,ny,3,n_gauss_n), Gzn(nz,nz,3,n_gauss_n))
      Gxn = 0.0_dp
      Gyn = 0.0_dp
      Gzn = 0.0_dp
    endif
    if(n_gauss_p .ne. 0) then 
      allocate(Gxp(nx,nx,3,n_gauss_p), Gyp(ny,ny,3,n_gauss_p), Gzp(nz,nz,3,n_gauss_p))
      Gxp = 0.0_dp
      Gyp = 0.0_dp
      Gzp = 0.0_dp
    endif
    ! early return if possible
    if(n_gauss_n .eq. 0 .and. n_gauss_p .eq. 0) return
    !---------------------------------------------------------------------------
    ! Harmonic-oscillator correction
    if(hocomform) then
        ! hbar x omega
        hbom  = 41.0 * (neutrons + protons)**(-1.0/3.0)
        ! 2m/hbar^2
        mhb = 2.0/(1.0/hbm(1)+1.0/hbm(2))
        ! B^{-1} = hbar * omega/m * A = 1/2 * A * hbar omega * 2m/hbar^2
        B = sqrt( 1.0/( 0.5 * hbom/mhb  * (neutrons + protons)))

        if(rplus_n.ne.0.0) then
            rplus_n = sqrt(rplus_n**2 - B**2)
        endif
        if(rmin_n.ne.0.0) then
            rmin_n   = sqrt(rmin_n**2 - B**2)
        endif

        if(rplus_p.ne.0.0) then
            rplus_p = sqrt(rplus_p**2 - B**2)
        endif
        if(rmin_p.ne.0.0) then
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
    	    call gamma_1D(Gxn(:,:,index,1), meshx, nx, rplus_n, sign * reduX)
    	    call gamma_1D(Gyn(:,:,index,1), meshy, ny, rplus_n, sign * reduY)
    	    call gamma_1D(Gzn(:,:,index,1), meshz, nz, rplus_n, sign * reduZ)
    	  endif
    	  if(rmin_n .ne. 0.0_dp) then
    	    call gamma_1D(Gxn(:,:,index,2), meshx, nx, rmin_n,  sign * reduX)
    	    call gamma_1D(Gyn(:,:,index,2), meshy, ny, rmin_n,  sign * reduY)
    	    call gamma_1D(Gzn(:,:,index,2), meshz, nz, rmin_n,  sign * reduZ)
    	  endif
    	  ! Protons
    	  if(rplus_p .ne. 0.0_dp) then
    	    call gamma_1D(Gxp(:,:,index,1), meshx, nx, rplus_p, sign * reduX)
    	    call gamma_1D(Gyp(:,:,index,1), meshy, ny, rplus_p, sign * reduY)
    	    call gamma_1D(Gzp(:,:,index,1), meshz, nz, rplus_p, sign * reduZ)
    	  endif
    	  if(rmin_p .ne. 0.0_dp) then
    	    call gamma_1D(Gxp(:,:,index,2), meshx, nx, rmin_p,  sign * reduX)
    	    call gamma_1D(Gyp(:,:,index,2), meshy, ny, rmin_p,  sign * reduY)
    	    call gamma_1D(Gzp(:,:,index,2), meshz, nz, rmin_p,  sign * reduZ)
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
    	    call gauss_1D(Gxn(:,:,index,1), meshx, nx, rplus_n, sign * reduX)
    	    call gauss_1D(Gyn(:,:,index,1), meshy, ny, rplus_n, sign * reduY)
    	    call gauss_1D(Gzn(:,:,index,1), meshz, nz, rplus_n, sign * reduZ)
    	  endif
    	  if(rmin_n .ne. 0.0_dp) then
    	    call gauss_1D(Gxn(:,:,index,2), meshx, nx, rmin_n,  sign * reduX)
    	    call gauss_1D(Gyn(:,:,index,2), meshy, ny, rmin_n,  sign * reduY)
    	    call gauss_1D(Gzn(:,:,index,2), meshz, nz, rmin_n,  sign * reduZ)
    	  endif
    	  ! Protons
    	  if(rplus_p .ne. 0.0_dp) then
    	    call gauss_1D(Gxp(:,:,index,1), meshx, nx, rplus_p, sign * reduX)
    	    call gauss_1D(Gyp(:,:,index,1), meshy, ny, rplus_p, sign * reduY)
    	    call gauss_1D(Gzp(:,:,index,1), meshz, nz, rplus_p, sign * reduZ)
    	  endif
    	  if(rmin_p .ne. 0.0_dp) then
    	    call gauss_1D(Gxp(:,:,index,2), meshx, nx, rmin_p,  sign * reduX)
    	    call gauss_1D(Gyp(:,:,index,2), meshy, ny, rmin_p,  sign * reduY)
    	    call gauss_1D(Gzp(:,:,index,2), meshz, nz, rmin_p,  sign * reduZ)
    	  endif
    	enddo
    endif
 end subroutine construct_folding_matrices


end module folding
