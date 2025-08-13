module Coulombmod
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
 ! Module that solves the Coulomb problem for the charge density (and is able
 ! to construct said charge density.)
 !==============================================================================
 ! Hephaestos keywords
 !
 !   REDUX  : $REDUX
 !   REDUY  : $REDUY
 !   REDUZ  : $REDUZ
 !
 !   FULLX  : $FULLX
 !   FULLY  : $FULLY
 !   FULLZ  : $FULLZ
 !
 !==============================================================================
 !
 ! Reminder about Coultreatment
 !
 !   (0) => No Coulomb
 !   (1) => Ordinary: direct and exchange
 !   (2) => Only direct
 !==============================================================================

 use geninfo
 use densities
 use parameterization
 use timing

 implicit none

 public
 !------------------------------------------------------------------------------
 !Precision required of the Coulomb Solvers; gets dynamically set as a function
 ! of dx down below
 real(KIND=dp)              :: Prec
 !------------------------------------------------------------------------------
 ! Number of boundary conditions to put on all sides of the box.
#if(USE_Periodic==1)
 integer :: BC = 0
#else
 integer :: BC = 2
#endif
 !------------------------------------------------------------------------------
 ! Arrays containing,
 ! 1) the values of the spherical harmonics on the extended mesh and
 ! 2) the value of the radial coordinate r on the extended mesh.
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable, target :: spherical_harmonics_coulomb(:,:,:,:,:,:)
 real(KIND=dp), allocatable, target :: r(:,:,:)
 !------------------------------------------------------------------------------
 ! Coordinates of the mesh in the enlarged coulomb box.
 real(KIND=dp), allocatable :: coulmeshx(:), coulmeshy(:), coulmeshz(:)
 real(KIND=dp), allocatable :: coulgrid(:,:)
 !------------------------------------------------------------------------------
 ! Maximum l of the multipole moments to use in the boundary conditions
 ! Currently hardcoded at 8: does not cost anything CPU-time wise and
 ! has been shown to be sufficient in MOCCa.
 integer, parameter :: max_moment_coulomb=8
 !------------------------------------------------------------------------------
 ! Offsets for the Coulomb box.
 integer :: coul_offset_x, coul_offset_y, coul_offset_z
 !------------------------------------------------------------------------------
 ! Flag indicating whether or not the direct and exchange Coulomb potentials
 ! were read from file. This is by default .false.; it should only happen when
 ! reading .pot files.
 logical :: Coulomb_read_from_file = .false.

 !------------------------------------------------------------------------------
 ! Coefficients of the Coulomb laplacian
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable :: CoulCoefs(:)
 real(KIND=dp),parameter,dimension(3) :: CoulCoefs_3=(/ &
 &  1.0_dp,-2.0_dp,1.0_dp /)
 real(KIND=dp),parameter,dimension(5) :: CoulCoefs_5=(/ &
 & -1.0_dp/12.0_dp,4.0_dp/3.0_dp,-5.0_dp/2.0_dp,4.0_dp/3.0_dp, -1.0_dp/12.0_dp/)
 real(KIND=dp),parameter,dimension(7) :: CoulCoefs_7 =(/ &
 &     1.0_dp/90.0_dp, -3.0_dp/20.0_dp, 3.0_dp/2.0_dp, -49.0_dp/18.0_dp,       &
 &     3.0_dp/2.0_dp,  -3.0_dp/20.0_dp, 1.0_dp/90.0_dp/)
 real(KIND=dp),parameter,dimension(9) :: CoulCoefs_9 =(/ &
 & -9.0_dp/8064.0_dp, 128.0_dp/8064.0_dp, -1008.0_dp/8064.0_dp, 1.0_dp,     &
 & -14350.0_dp/8064.0_dp, 1.0_dp, -1008.0_dp/8064.0_dp, 128.0_dp/8064.0_dp, &
 & -9.0_dp/8064.0_dp /)

 interface SolveCoulomb_worker
  module procedure solve_coulomb_worker_real
  module procedure solve_coulomb_worker_complex
 end interface

contains

 subroutine SolveCoulomb(R,F,sx,sy,sz)
    !---------------------------------------------------------------------------
    ! Master routine to solve the Coulomb problem for a given source-density.
    !
    ! Input :
    !     R       : density vector with a precalculated charge density
    !     sx/sy/sz: reflection symmetries of the charge density
    !               these are explicit inputs to accomodate perturbed mean-field
    !               densities whose charge density need not have the same
    !               symmetry properties
    ! Output:
    !     F : potential vector; only the Coulomb fields are modified
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Note: this routine does not calculate the charge density based on the
    !       the proton and neutron densities. This is assumed to have been
    !       done before.
    !---------------------------------------------------------------------------
    type(DensityVector), intent(in)      :: R
    type(PotentialVector), intent(inout) :: F
    integer, intent(in)                  :: sx, sy, sz

    ! Initialize all of the arrays.
    call setupcoulomb(F)
    ! Solve Poissons equation for the charge density
    call SolveCoulomb_worker(R%chargedensity, F%CoulombPotential, F%ExchangePotential, &
    &                        sx, sy, sz)
    ! Perform a folding of the potentials if needed
    call Obtain_folded_potentials(F)

 end subroutine SolveCoulomb

 subroutine solve_coulomb_worker_real(chargedensity, CoulombPotential, ExchangePotential, &
 &                                    sx, sy, sz)
    !-------------------------------------------------------------------------------------
    ! TODO: DOcument
    !
    !-------------------------------------------------------------------------------------
    real(KIND=dp), allocatable, intent(in)    :: Chargedensity(:,:,:)
    real(KIND=dp), allocatable, intent(inout) :: CoulombPotential(:,:,:)
    real(KIND=dp), allocatable, intent(inout) :: ExchangePotential(:,:,:)
    integer, intent(in)                       :: sx, sy, sz

    call SolveCoulomb_solver(chargedensity, CoulombPotential, ExchangePotential,sx,sy,sz)

 end subroutine solve_coulomb_worker_real

  subroutine solve_coulomb_worker_complex(chargedensity, CoulombPotential, ExchangePotential, &
 &                               sx, sy, sz)
    !-------------------------------------------------------------------------------------
    ! TODO: DOcument
    !
    !-------------------------------------------------------------------------------------
    complex(KIND=dp), allocatable, intent(in)    :: Chargedensity(:,:,:)
    complex(KIND=dp), allocatable, intent(inout) :: CoulombPotential(:,:,:)
    complex(KIND=dp), allocatable, intent(inout) :: ExchangePotential(:,:,:)
    integer, intent(in)                          :: sx, sy, sz

    real(KIND=dp), allocatable :: Re_CD(:,:,:), Im_CD(:,:,:)
    real(KIND=dp), allocatable :: Re_CP(:,:,:), Im_CP(:,:,:)
    real(KIND=dp), allocatable :: Re_EX(:,:,:), Im_EX(:,:,:)

    Re_CD = DBLE(ChargeDensity)    ; Im_CD = IMAG(ChargeDensity)
    Re_CP = DBLE(CoulombPotential) ; Im_CP = IMAG(CoulombPotential)
    Re_EX = DBLE(ExchangePotential); Im_EX = IMAG(ExchangePotential)

    ! Solve the real part ...
    print *, 'SOLVING REAL PART'
    call SolveCoulomb_solver(Re_CD, Re_CP, Re_Ex, sx, sy, sz)
    ! ... and the imaginary part of the Coulomb equation
    print *, 'SOLVING IMAGINARY PART'
    call SolveCoulomb_solver(Im_CD, Im_CP, Im_Ex, sx, sy, sz)

    ! ... and sum the results
    CoulombPotential  = CMPLX(Re_CD, Im_CD)
    ExchangePotential = CMPLX(Re_EX, Im_EX)

 end subroutine solve_coulomb_worker_complex

 subroutine SolveCoulomb_solver(chargedensity, CoulombPotential, ExchangePotential, &
 &                               sx, sy, sz)
    !----------------------------------------------------------------------------------------
    ! TODO: document this worker routine
    !
    !
    !
    !
    !----------------------------------------------------------------------------------------
    real(KIND=dp), allocatable, intent(in)    :: Chargedensity(:,:,:)
    real(KIND=dp), allocatable, intent(inout) :: CoulombPotential(:,:,:)
    real(KIND=dp), allocatable, intent(inout) :: ExchangePotential(:,:,:)
    integer, intent(in)                       :: sx, sy, sz

    real(KIND=dp), allocatable      :: source(:,:,:)
    integer                         :: i,j,k,ii

    call start_timer(T_coulomb)

    if(.not.allocated(CoulCoefs)) then
     ! Determine the coefficients of the finite difference scheme
     select case(coulorder)
       case(1)
        BC = 1
        allocate(CoulCoefs(3)) ; CoulCoefs = CoulCoefs_3
       case(2)
        BC = 2
        allocate(CoulCoefs(5)) ; CoulCoefs = CoulCoefs_5
       case(3)
        BC = 3
        allocate(CoulCoefs(7)) ; CoulCoefs = CoulCoefs_7
       case(4)
        BC = 4
        allocate(CoulCoefs(9)) ; CoulCoefs = CoulCoefs_9
        CoulCoefs = Coulcoefs * (8064.0_dp)/(5040.0_dp) !Legacy from CR8
       case DEFAULT
        call stp('This order for the Coulomb discretisation is not supported.')
     end select
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! If periodic boundary conditions are active, we need no additional points
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#if(USE_Periodic==1)
     BC=0
#endif
    endif

    ! Determine the offsets of the original mesh inside the larger Coulomb mesh
    coul_offset_x = BC ; coul_offset_Y = BC ; coul_offset_z = BC

    ! If any given axis is not actually represented in computer memory, the
    ! offset of the mesh in that direction is zero.
$REDUX  coul_offset_x = 0
$REDUY  coul_offset_y = 0
$REDUZ  coul_offset_z = 0

    if(.not.allocated(Source)) then
        allocate(Source(nx+BC+coul_offset_x, &
        &               ny+BC+coul_offset_y, &
        &               nz+BC+coul_offset_z))
        Source = 0.0_dp
    endif

    !---------------------------------------------------------------------------
    ! Early return if possible
    if(coultreatment.eq.0) then
       call stop_timer(T_coulomb)
       return
    endif
    !---------------------------------------------------------------------------
    ! Set up the source term: - 4 * pi * charge_density
    ! Note that this is set up in the middle of the box, i.e. no source density
    ! at the edges of the Coulomb box
    Source = 0.0_dp
    do k=1,nz
      do j=1,ny
        do i=1,nx
          Source(i+coul_offset_x,j+coul_offset_y,k+coul_offset_z) = &
          &                                  -4*pi*e2*Chargedensity(i,j,k)
        enddo
      enddo
    enddo
    !---------------------------------------------------------------------------
    ! Set the boundary condition if dealing with non-periodic boundary conditions
    ! In the peridic case, these are automatically taken care of
#if(USE_Periodic == 0)
    call coulomb_bound(Source, sx, sy, sz, CoulombPotential)
#endif
    !---------------------------------------------------------------------------
    ! Solve for the direct coulomb potential
    ! Note that the symmetry properties (+1,+1,+1) are never changed:
    ! Hephaestos modifies directly the Coulomb_Laplacian routine when necessary
    call ConjugGrad (CoulombPotential,Source, sx,sy,sz,100,.false.,prec)

    if(Coultreatment.eq.1) then
      !-------------------------------------------------------------------------
      ! Exchange potential in Slater approximation
      ExchangePotential =  &
      &           -(3.0/pi)**(1.0/3.0_dp)*e2*(ChargeDensity**(1.0_dp/3.0_dp))
    else
      !-------------------------------------------------------------------------
      ! No Coulomb Exchange
      ExchangePotential = 0.0
    endif

    call stop_timer(T_coulomb)
    deallocate(source)
 end subroutine SolveCoulomb_solver

 subroutine Obtain_folded_potentials(F)
    !---------------------------------------------------------------------------
    ! Obtain the folded Coulomb potentials (direct and exchange) if needed.
    ! It is a separate routine from SolveCoulomb because it should also be
    ! callable from the routine to read potentials.
    !---------------------------------------------------------------------------
    type(PotentialVector), intent(inout) :: F

    !if(any(protonsize .ne. 0.0_dp) .or. any(neutronsize.ne.0.0_dp)) then
    !  if(nucleonsize_selfconsistent) then
    !     F%FoldedCoul    =FoldCoulombPotential(F%CoulombPotential(             &
    !     &                                    coul_offset_x+1:coul_offset_x+nx,&
    !     &                                    coul_offset_y+1:coul_offset_y+ny,&
    !     &                                    coul_offset_z+1:coul_offset_z+nz))
!
!         if(coultreatment .eq. 1) then
!           F%FoldedExchange=FoldCoulombPotential(F%ExchangePotential(          &
!           &                                  coul_offset_x+1:coul_offset_x+nx,&
!           &                                  coul_offset_y+1:coul_offset_y+ny,&
!           &                                  coul_offset_z+1:coul_offset_z+nz))
!         else
!           if(.not. allocated(F%FoldedExchange)) then
!            allocate (F%FoldedExchange(nx,ny,nz,2))
!           endif
!           F%foldedexchange = 0.0d0
!         endif
!      endif
!    endif

 end subroutine Obtain_folded_potentials

 function FoldCoulombPotential(pot) result(Folded)
    !---------------------------------------------------------------------------
    ! Obtain the folded Coulomb potential, for use in the single-particle
    ! hamiltonian when finite size corrections are included selfconsistently.
    !---------------------------------------------------------------------------
    use Folding

    real(KIND=dp), intent(in)  :: pot(:,:,:)
    real(KIND=dp), allocatable :: Folded(:,:,:,:)

    allocate(folded(nx,ny,nz,2)) ; folded = 0.0

    if(protonsize(1).gt.0.0) then
        ! Fold the potential with the Gaussian of positive sign for protons
        folded(:,:,:,2) = folded(:,:,:,2) + FoldGaussian( pot,GaussX(:,:,1,2), &
        &                                                     GaussY(:,:,1,2), &
        &                                                     GaussZ(:,:,1,2), &
        &                                                            nx, ny, nz)
    endif
    if(protonsize(2).gt.0.0) then
        ! Fold the potential with the Gaussian of negative sign for protons
        folded(:,:,:,2) = folded(:,:,:,2) - FoldGaussian( pot,GaussX(:,:,2,2), &
        &                                                     GaussY(:,:,2,2), &
        &                                                     GaussZ(:,:,2,2), &
        &                                                            nx, ny, nz)
    endif

    !---------------------------------------------------------------------------
    ! Note that, if the neutron charge form factor is included, they feel a
    ! Coulomb potential as well!
    if(all(neutronsize.eq.0.0) .or. neutroncoulomberror) return
    if(neutronsize(1).gt.0.0) then
        ! Fold the potential with the Gaussian of positive sign for neutrons
        folded(:,:,:,1) = folded(:,:,:,1) + FoldGaussian( pot,GaussX(:,:,1,1), &
        &                                                     GaussY(:,:,1,1), &
        &                                                     GaussZ(:,:,1,1), &
        &                                                            nx, ny, nz)
    endif
    if(neutronsize(2).gt.0.0) then
        ! Fold the potential with the Gaussian of negative sign for neutrons
        folded(:,:,:,1) = folded(:,:,:,1) - FoldGaussian( pot,GaussX(:,:,2,1), &
        &                                                     GaussY(:,:,2,1), &
        &                                                     GaussZ(:,:,2,1), &
        &                                                            nx, ny, nz)
    endif

 end function FoldCoulombPotential

 subroutine SetupCoulomb(F)
    !---------------------------------------------------------------------------
    ! Initialize the entire module and the fields in the potentialvector
    !
    ! Input:
    !    F : potentialvector to be initialized
    !---------------------------------------------------------------------------
    use sphericalharmonics, only: generate_spherical_harmonics
    use moments,            only: QuantisationAxis, SecondaryAxis
    use folding
    use vectors

    type(PotentialVector), intent(inout) :: F
    integer       :: i,j,k, ox, oy, oz

    ! Determine the offsets of the original mesh inside the larger Coulomb mesh
    coul_offset_x = BC ; coul_offset_Y = BC ; coul_offset_z = BC

    ! If any given axis is not represented, the offset of the mesh in that
    ! direction is zero.
$REDUX  coul_offset_x = 0
$REDUY  coul_offset_y = 0
$REDUZ  coul_offset_z = 0

    ox = nx+BC+coul_offset_x
    oy = ny+BC+coul_offset_y
    oz = nz+BC+coul_offset_z

    !---------------------------------------------------------------------------
    ! Allocate the CoulombPotential array on the full Coulomb mesh
    if(.not.allocated(F%CoulombPotential)) then
      allocate(F%CoulombPotential(ox,oy,oz))
      F%CoulombPotential = 0.0_dp
      ! The exchange potential is only defined on the original mesh
      allocate(F%ExchangePotential(nx,ny,nz))
      F%ExchangePotential = 0.0_dp
    endif
    !---------------------------------------------------------------------------
    ! Precision desired of the Coulomb solver
    Prec = 1.d-12/(dx**3)

    !---------------------------------------------------------------------------
    ! Set-up the value of r on the Coulomb mesh.
#if(USE_Periodic == 0)
    if(.not. allocated(r)) then
      call inimesh(coulmeshx,coulmeshy,coulmeshz,nx+BC+coul_offset_x, &
      &                                          ny+BC+coul_offset_y, &
      &                                          nz+BC+coul_offset_z, &
      &                                          coulgrid,0.0d0,0.0d0,0.0d0)

      allocate(r(ox,oy,oz))  ;  r = 0.0_dp
      allocate(spherical_harmonics_coulomb(ox,oy,oz,0:max_moment_coulomb,0:max_moment_coulomb,2))
      spherical_harmonics_coulomb = 0.0_dp

      do k=1,oz
        do j=1,oy
          do i=1,ox
            r(i,j,k) = sqrt(coulmeshx(i)**2 + coulmeshy(j)**2 + coulmeshz(k)**2)
          enddo
        enddo
      enddo

      call generate_spherical_harmonics(max_moment_coulomb,ox,oy,oz,                         &
      &                          coulmeshx,coulmeshy, coulmeshz,spherical_harmonics_coulomb, &
      &                          QuantisationAxis,SecondaryAxis)
    endif
#endif
    !---------------------------------------------------------------------------
    ! If we account for the finite extent of the charge of the nucleus, then
    ! we need to fold densities and potentials with gaussians. This sets up the
    ! required matrices.
    !
    ! Note: this little piece of code is duplicated, since in different
    !       runmodes of the code different Coulomb routines get called in
    !       different order; this makes sure we get no segfaults.
    !
    !---------------------------------------------------------------------------
    if(any(protonsize .ne. 0.0_dp) .or. any(neutronsize.ne.0.0_dp)) then
        if(.not.allocated(Gaussx)) then
            allocate(Gaussx(nx,nx,2,2), Gaussy(ny,ny,2,2), Gaussz(nz,nz,2,2))
            Gaussx = 0.0 ;  Gaussy = 0.0 ; Gaussz = 0.0
        endif
        !-----------------------------------------------------------------------
        ! Construct Gauss matrices
        call ConstructFoldingMatrices(Gaussx,Gaussy,Gaussz,sx_rho, sy_rho, sz_rho)
    endif

 end subroutine SetupCoulomb

 subroutine coulomb_bound(source, sx, sy, sz, coulomb_potential)
    !-----------------------------------------------------------------------------
    ! Sets correct boundary conditions of the Coulomb potential at the most extreme
    ! points of the Coulomb mesh.
    !
    ! Far from the source, the potential due to a charge at r' experienced at r 
    ! (|r'| > |r|) can be written as an expansion in terms of solid harmonics: 
    !
    !   (|r -r'|)^-1 = \sum_{l} \sum_{m=-l}^{l} (-1)^m I^{-m}_l(r) R^{m}_l(r') 
    ! 
    ! where the notation is that of Wikipedia; i.e.
    !
    ! - the irregular solid harmonics I^{m}_l(r) are defined as
    !         I^{m}_l(r)  = N Y^m_l(\theta,\phi) / r^{l+1} 
    ! - the regular solid harmonics R^{m}_l(r') are defined as
    !         R^{m}_l(r') = N Y^m_l(\theta',\phi') r^{l}
    ! - the spherical harmonics Y^m_l(\theta,\phi)
    ! - N is a normalisation constant 
    !      N = \sqrt{4 \pi /(2l+1)}
    !
    ! For the Coulomb potential of a charge distribution \rho_c, this means 
    ! 
    !  V(r) = \int d^3r' (|r -r'|)^-1 \rho_c(r')
    !       = \sum_l \sum_{m=-l}^{l} (-1)^m I^{-m}_l(r) Q_lm
    !
    ! with Q_lm = \int d^3r' R^{m}_l(r') \rho_c(r').
    !
    ! Because Y^{-m}_l = (-1)^m [Y^m_l(\theta,\phi)]^* and - as a consequence -
    ! Q_l(-m) = (-1)^m Q_lm^*, we can rewrite the potential as
    !
    !  V(r)= \sum_l I^{0}_0(r) Q_00 + 
    !        \sum_{m=1}^{l} (-1)^m I^{-m}_l(r) Q_lm + (-1)^(-m) I^{m}_l(r) Q_l(-m)
    !      = \sum_l I^{0}_0(r) Q_00 + 
    !        \sum_{m=1}^{l} (-1)^m I^{-m}_l(r) Q_lm + (-1)^(-m) I^{-m,*}_l(r) Q^*_lm
    !      = \sum_l I^{0}_0(r) Q_00 + 
    !        \sum_{m=1}^{l}  I^{+m,*}_l(r) Q_lm +  I^{+m}_l(r) Q^*_lm
    !      = \sum_l I^{0}_0(r) Q_00 + 2 \sum_{m=1}^{l} \Re [ I^{+m}_l(r) Q_lm ] 
    !      = \sum_l I^{0}_0(r) Q_00 
    !              + 2 \sum_{m=1}^{l} \Re [ I^{+m}_l(r) ] \Re [ Q_lm ] 
    !              - 2 \sum_{m=1}^{l} \Im [ I^{+m}_l(r) ] \Im [ Q_lm ] 
    ! 
    ! which is a useable expression for a code that deals with the real and
    ! imaginary parts of solid harmonics.
    !
    ! Slightly easier to code is the expression
    !
    !      =  \sum_{m=0}^{l} (2 - \delta_{m 0}) \Re [ I^{+m}_l(r) ] \Re [ Q_lm ] 
    !       - \sum_{m=1}^{l} (2 - \delta_{m 0}) \Im [ I^{+m}_l(r) ] \Im [ Q_lm ] 
    !
    ! where I used that Im (Y_00) = 0.
    !
    ! TODO: document normalisation!
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   source    : source charge density to compute boundary conditions for
    !               (the factor 4*pi*e2 should be included!)
    !   sx/sy/sz  : symmetry properties of the source charge density
    ! Output:
    !   coulomb_potential : the potential with boundary conditions applied.
    !---------------------------------------------------------------------------
    use folding
    use vectors
    use sphericalharmonics, only : generate_spherical_harmonics
    use moments,            only : QuantisationAxis, SecondaryAxis
    use moments,            only : figure_out_multipole_moments

    real(KIND=dp), intent(in)    :: source(:,:,:)
    integer, INTENT(IN)          :: sx, sy, sz
    real(KIND=dp), intent(inout) :: coulomb_potential(:,:,:)

    integer                   :: moment_list(0:max_moment_coulomb,0:max_moment_coulomb,0:1)
    real(KIND=dp)             :: fac, Qlm
    integer                   :: i,j,k,l,m, im, ox, oy, oz
    logical                   ::  condition

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !  First, set the potential to zero on the edges of the Coulomb mesh
    ox = nx+BC+coul_offset_x
    oy = ny+BC+coul_offset_y
    oz = nz+BC+coul_offset_z

    do k=1,oz
      do j=1,oy
        do i=1,ox

          condition = .false.
          if(i.gt.nx) condition =.true.
          if(j.gt.ny) condition =.true.

$REDUZ     if(k.gt.nz) condition =.true.

$FULLZ     if(k.le.BC)    condition =.true.
$FULLZ     if(k.gt.nz+BC) condition =.true.

          if(condition) then
              coulomb_potential(i,j,k) = 0.0d0
          endif
        enddo
      enddo
    enddo
 
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Figure out what multipole moments of the source density will not vanish.
    !  moment_list(l,m,k) = 0 => restricted by symmetry, 
    !              | | |           should not be calculated
    !              | | |
    !  moment_list(l,m,k) = 1 => should be calculated
    !              | | |
    !              | | -> k : real (0) or imaginary(1) part
    !              | ---> m : second characteristic number
    !              -----> l : first characteristic number
    !                      
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    moment_list = figure_out_multipole_moments(sx,sy,sz,max_moment_coulomb,&
    &                                        quantisationaxis,secondaryaxis) 

    do im = 0,1 
      do l=0, max_moment_coulomb
        do m=0, l
          Qlm = 0

          ! If the multipole moment is restricted by symmetry, we move on
          if(moment_list(l,m,im).eq.0) cycle 

          ! calculate the multipole moment Q_lm of the source density
          do k=1,oz
            do j=1,oy
              do i=1,ox
                Qlm = Qlm - Source(i,j,k) * spherical_harmonics_coulomb(i,j,k,l,m,Im+1)
              enddo
            enddo
          enddo

          Qlm = Qlm * dv/(2*l+1) ! volume element and normalisation
          print *, 'Multipole moment in Coulomb', l, m, im, ' = ', Qlm, sx, sy, sz

          ! The prefactor +/-(2 - \delta_{m 0})
          fac = 1
          if(m  .ne. 0) fac =    2 ! Real parts of Q_lm and Q_l(-m) are identical
          if(im .eq. 1) fac = -fac ! Imaginary parts obtain a minus sign 

          do k=1,oz
            do j=1,oy
              do i=1,ox

                condition = .false.
$REDUX          if(i.gt.nx) condition =.true.
$REDUY          if(j.gt.ny) condition =.true.
$REDUZ          if(k.gt.nz) condition =.true.

$FULLX          if(i.le.BC)    condition =.true.
$FULLX          if(i.gt.nx+BC) condition =.true.

$FULLY          if(j.le.BC)    condition =.true.
$FULLY          if(j.gt.ny+BC) condition =.true.

$FULLZ          if(k.le.BC)    condition =.true.
$FULLZ          if(k.gt.nz+BC) condition =.true.
                if(condition) then
                  coulomb_potential(i,j,k) = coulomb_potential(i,j,k) +          &
                  &       fac*Qlm*spherical_harmonics_coulomb(i,j,k,l,m,Im+1)/(r(i,j,k)**(2*l+1))
                endif
              enddo
            enddo
          enddo

        enddo
      enddo
    enddo
 end subroutine coulomb_bound

 function CoulombEnergy_direct(R, F) result(CEnergy)
    !---------------------------------------------------------------------------
    ! Calculate the (direct) electrostatic energy of the system.
    !
    ! Input:
    !     R: Densityvector (which contains the charge density)
    !     F: Potentialvector (which contains the Coulomp potential)
    !---------------------------------------------------------------------------
    use vectors

    type(DensityVector), intent(in)   :: R
    type(PotentialVector), intent(in) :: F
    real(KIND=dp) :: CEnergy
    integer       :: i,j,k, ox, oy, oz

    ox = coul_offset_x ; oy = coul_offset_y ; oz = coul_offset_z

    CEnergy = 0.0_dp
    do k=1,nz
        do j=1,ny
            do i=1,nx
                CEnergy = CEnergy + R%chargedensity(i,j,k) *   &
                &                   F%CoulombPotential(i+ox,j+oy,k+oz)
            enddo
        enddo
    enddo
    CEnergy = CEnergy * dv * 0.5_dp
 end function CoulombEnergy_Direct

 function CoulombEnergy_Exchange(R) result(CEnergy)
    !---------------------------------------------------------------------------
    ! Calculate the (exchange) electrostatic energy of the system in the Slater
    ! approximation.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !    R :  a set of mean-field densities
    !---------------------------------------------------------------------------
    use vectors

    type(DensityVector), intent(in) ::R
    real(KIND=dp) :: factor, Cenergy

    Cenergy = 0.0_dp
    if(coultreatment.ne.1) return

    factor  = -0.75_dp*(3/pi)**(1/3._dp)*e2*dv
    Cenergy = factor*sum(R%chargedensity**(4.0/3.0))

 end function CoulombEnergy_Exchange

 subroutine ConjugGrad (Solution,SourceTerm,sx,sy,sz,MaxIteration,iprint,Precis)
    !---------------------------------------------------------------------------
    ! This subroutine solves the Coulomb problem. The technique is identical to
    ! the ones employed in EV8 and CR8 and is a straight-forward conjugate
    ! gradients algorithm.
    !
    !---------------------------------------------------------------------------
    ! Actual Algorithm
    !
    !  Define the following quantities (indices k refer to iteration k):
    !            "Proton density "        \rho
    !            "Coulomb Potential"      \phi_k
    !            "Actual Residue"         r_k = \rho - \Delta \phi
    !            "Iterative Residue"      p_k
    !            "Conjugate coefficient"  a_k = (r_k^T r_k)/(p_k^T \Delta p_k^T)
    !            "Evolving Coefficient"   b_k = (r_(k+1)^T r_(k+1))/(r_k^T r_k)
    !  Note:
    !
    !  The actual algorithm goes as follows:
    !    1) From a starting potential CoulombPotential (0 if this is the first
    !       time this algorithm is invoked),
    !       compute:
    !                r_0 = -e/(epsilon)*Rho - \Delta \Phi
    !                p_0 = r_0
    ! |--2)  Then compute
    ! |              a_k       = (r_k^T r_k)/(p_k^T \Delta p_k)
    ! |              phi_(k+1) = phi_k + a_k p_k
    ! |              r_(k+1)   = r_k - a_k \Delta p_k
    ! |  3) Check if sum(|r_k|) is smaller than the desired precision.
    ! |      If so, phi_(k+1) is the desired potential.
    ! |      If not compute:
    ! |              b_k       = (r_(k+1)^T r_(k+1))/(r_k^T r_k)
    ! |              p_(k+1)   = r_(k+1) + b_k p_k
    ! |
    ! |--4) Go to 2).
    !
    !---------------------------------------------------------------------------

    1 format ("Coulomb Calculation converged after", i4, " Iterations.")
    2 format ("Coulomb Calculation did not converge after", i4,                &
    &         " Iterations. Residual:", e15.8)

    integer,intent(in)                :: MaxIteration, sx, sy,sz
    real(KIND=dp), intent(in)         :: SourceTerm(:,:,:)
    real(KIND=dp), intent(inout)      :: Solution(:,:,:)
    real(KIND=dp), intent(in),optional:: Precis
    logical, intent(in)               :: iprint

    integer       :: iteration
    real(KIND=dp) :: p_k(nx+BC+coul_offset_x, ny+BC+coul_offset_y, nz+BC+coul_offset_z)
    real(KIND=dp) :: Temp(nx+BC+coul_offset_x, ny+BC+coul_offset_y, nz+BC+coul_offset_z)
    real(KIND=dp) :: Residual(nx+BC+coul_offset_x, ny+BC+coul_offset_y, nz+BC+coul_offset_z)
    real(KIND=dp) :: PoissonNorm, Integral, a_k, c_k
    real(KIND=dp) :: NewPoissonNorm

    Residual = - coulomb_laplacian(Solution,sx,sy,sz)
    Residual = Residual + SourceTerm
    !---------------------------------------------------------------------------
    !The variable p_k is the conjugate direction. It starts out equal to our
    ! initial Residual.
    p_k=Residual
    !Poissonnorm is r_(k+1)^T r_(k+1)
    PoissonNorm = sum(Residual**2)

    ! safeguard against the possiblity of their being no charge at all
    if(PoissonNorm .eq. 0.0) return

    do iteration = 1,MaxIteration
      !Applying lagrangian to p_k
      Temp = coulomb_laplacian(p_k,sx,sy,sz)

      !Integral is p_k^T \Delta p_k^T
      Integral = sum(Residual*Temp)
      a_k = PoissonNorm/Integral

      !Increment The Coulomb Potential
      Solution = Solution + a_k*p_k

      !Increment the Poisson Equation
      Residual = Residual - a_k*Temp
      !NewPoissonNorm = zz2
      NewPoissonNorm = sum(Residual**2)
      if(sum(Residual**2).le.Precis .and. iteration.gt.5) then
        !Check if convergence is reached.
        if(iprint) print 1, Iteration
        exit
      elseif(Iteration.eq.MaxIteration) then
        if(iprint) print 2, Iteration, PoissonNorm
      endif
      !c_k = zz2/zz1
      c_k = NewPoissonNorm/PoissonNorm
      PoissonNorm = NewPoissonNorm
      !Calculate the next conjugate gradient
      p_k = Residual    + c_k*p_k

      ! Diagnostic printing
      print *, 'Coul, it',  iteration, PoissonNorm, sx, sy, sz
    enddo

    return
  end subroutine ConjugGrad

  function coulomb_laplacian(f, sx, sy, sz) result(lf)
    !---------------------------------------------------------------------------
    ! Subroutine applying a finite difference operator (of order two) to the
    ! function f. Note that this function should be defined on the box + BC
    ! points in every direction, as boundary conditions are necessary.
    !
    ! Note that these boundary conditions are not touched by this procedure.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) ::  f(:,:,:)
    real(KIND=dp), allocatable:: lf(:,:,:)
    integer, intent(in)       :: sx, sy, sz
    integer                   :: i,j,k,l, ox, oy, oz, trash

    ! Statement to stop the compiler complaining about unused dummy arguments
    trash = sx ; trash = sy ; trash = sz

    ox = coul_offset_x ; oy = coul_offset_y ; oz = coul_offset_z

    allocate(lf(nx+BC+ox, ny+BC+oy, nz+BC+oz)) ;  lf = 0.0_dp

#if(USE_Periodic==0)
    !---------------------------------------------------------------------------
    ! X-direction
    do k=oz+1,oz+nz
        do j=oy+1,oy+ny
            do i=1+BC,nx+ox
              do l=-BC,+BC
                 lf(i,j,k) = lf(i,j,k) + CoulCoefs(l+BC+1) * f(i+l,j,k)
              enddo
            enddo
        enddo
    enddo

$REDUX    do k=oz+1,oz+nz
$REDUX      do j=oy+1,oy+ny
$REDUX        do i=1,BC
$REDUX          ! Forwards and backwards difference without plane reflection
$REDUX          do l=-i+1,BC
$REDUX            lf(i,j,k) = lf(i,j,k) +      CoulCoefs(l+BC+1)  * f(i+l,j,k)
$REDUX          enddo
$REDUX          ! Backwards difference (with plane reflection)
$REDUX          do l=-BC,-i
$REDUX            lf(i,j,k) = lf(i,j,k) + sx * CoulCoefs(l+BC+1) * f(-i-l+1,j,k)
$REDUX          enddo
$REDUX        enddo
$REDUX      enddo
$REDUX    enddo

    !---------------------------------------------------------------------------
    ! Y-direction
    do k=oz+1,oz+nz
        do j=1+BC,ny+oy
            do i=ox+1,ox+nx
              do l=-BC,+BC
                 lf(i,j,k) = lf(i,j,k) + CoulCoefs(l+BC+1) * f(i,j+l,k)
              enddo
            enddo
        enddo
    enddo

$REDUY    do k=oz+1,oz+nz
$REDUY      do j=1,BC
$REDUY        do i=ox+1,ox+nx
$REDUY          ! Forwards and backwards difference without plane reflection
$REDUY          do l=-j+1,BC
$REDUY            lf(i,j,k) = lf(i,j,k) +      CoulCoefs(l+BC+1)  * f(i,j+l,k)
$REDUY          enddo
$REDUY          ! Backwards difference (with plane reflection)
$REDUY          do l=-BC,-j
$REDUY            lf(i,j,k) = lf(i,j,k) + sy * CoulCoefs(l+BC+1) * f(i,-j-l+1,k)
$REDUY          enddo
$REDUY        enddo
$REDUY      enddo
$REDUY    enddo

    !---------------------------------------------------------------------------
    ! Z-direction
    !---------------------------------------------------------------------------
    do k=1+BC,nz+oz
        do j=oy+1,oy+ny
            do i=ox+1,ox+nx
              do l=-BC,+BC
                 lf(i,j,k) = lf(i,j,k) + CoulCoefs(l+BC+1) * f(i,j,k+l)
              enddo
            enddo
        enddo
    enddo

$REDUZ    do k=1,BC
$REDUZ      do j=oy+1,oy+ny
$REDUZ        do i=ox+1,ox+nx
$REDUZ          ! Forwards and backwards difference without plane reflection
$REDUZ          do l=-k+1,BC
$REDUZ            lf(i,j,k) = lf(i,j,k) +      CoulCoefs(l+BC+1)  * f(i,j,k+l)
$REDUZ          enddo
$REDUZ          ! Backwards difference (with plane reflection)
$REDUZ          do l=-BC,-k
$REDUZ            lf(i,j,k) = lf(i,j,k) + sz * CoulCoefs(l+BC+1) * f(i,j, -k-l+1)
$REDUZ          enddo
$REDUZ        enddo
$REDUZ      enddo
$REDUZ    enddo

#else
    !NS: impose periodic BC of order=coulorder
    !---------------------------------------------------------------------------
    ! X-direction
    do k=1,nz
        do j=1,ny
            do i=1+coulorder,nx-coulorder
              do l=-coulorder,+coulorder
                lf(i,j,k) = lf(i,j,k) + CoulCoefs(l+coulorder+1) * f(i+l,j,k)
              enddo
            enddo
            do i=1,coulorder
              !derivatives without extrapoints from boundary conditions
              do l=-i+1,coulorder
                lf(i,j,k)     =     lf(i,j,k)+CoulCoefs(l+coulorder+1) *f(i+l,j,k)
                lf(nx+1-i,j,k)=lf(nx+1-i,j,k)+CoulCoefs(-l+coulorder+1)*f(nx+1-i-l,j,k)
              enddo
              !Boundary conditions
              do l=-coulorder,-i
$REDUX               lf(i,j,k)     =     lf(i,j,k)+CoulCoefs(l+coulorder+1) *f(1-i-l,j,k)
$REDUX               lf(nx+1-i,j,k)=lf(nx+1-i,j,k)+CoulCoefs(-l+coulorder+1)*f(nx+i+l,j,k)
$REDUX               cycle
                lf(i,j,k)     =     lf(i,j,k)+CoulCoefs(l+coulorder+1) *f(nx+i+l,j,k)
                lf(nx+1-i,j,k)=lf(nx+1-i,j,k)+CoulCoefs(-l+coulorder+1)*f(1-i-l,j,k)
              enddo
            enddo
        enddo
    enddo


    !---------------------------------------------------------------------------
    ! Y-direction
    do k=1,nz
        do i=1,nx
            do j=1+coulorder,ny-coulorder
              do l=-coulorder,+coulorder
                 lf(i,j,k) = lf(i,j,k) + CoulCoefs(l+coulorder+1) * f(i,j+l,k)
              enddo
            enddo
            do j=1,coulorder
              do l=-j+1,coulorder
                lf(i,j,k)     =     lf(i,j,k)+CoulCoefs(l+coulorder+1) *f(i,j+l,k)
                lf(i,ny+1-j,k)=lf(i,ny+1-j,k)+CoulCoefs(-l+coulorder+1)*f(i,ny+1-j-l,k)
              enddo
              do l=-coulorder,-j
$REDUY               lf(i,j,k)     =     lf(i,j,k)+CoulCoefs(l+coulorder+1) *f(i,1-j-l,k)
$REDUY               lf(i,ny+1-j,k)=lf(i,ny+1-j,k)+CoulCoefs(-l+coulorder+1)*f(i,ny+j+l,k)
$REDUY               cycle
                lf(i,j,k)     =     lf(i,j,k)+CoulCoefs(l+coulorder+1) *f(i,ny+j+l,k)
                lf(i,ny+1-j,k)=lf(i,ny+1-j,k)+CoulCoefs(-l+coulorder+1)*f(i,1-j-l,k)
              enddo
            enddo
        enddo
    enddo


    !---------------------------------------------------------------------------
    ! Z-direction
    !---------------------------------------------------------------------------
    do j=1,ny
        do i=1,nx
            do k=1+coulorder,nz-coulorder
              do l=-coulorder,+coulorder
                 lf(i,j,k) = lf(i,j,k) + CoulCoefs(l+coulorder+1) * f(i,j,k+l)
              enddo
            enddo
            do k=1,coulorder
              do l=-k+1,coulorder
                lf(i,j,k)     =     lf(i,j,k)+CoulCoefs(l+coulorder+1) *f(i,j,k+l)
                lf(i,j,nz+1-k)=lf(i,j,nz+1-k)+CoulCoefs(-l+coulorder+1)*f(i,j,nz+1-k-l)
              enddo
              do l=-coulorder,-k
$REDUZ               lf(i,j,k)     =     lf(i,j,k)+CoulCoefs(l+coulorder+1) *f(i,j,1-k-l)
$REDUZ               lf(i,j,nz+1-k)=lf(i,j,nz+1-k)+CoulCoefs(-l+coulorder+1)*f(i,j,nz+k+l)
$REDUZ               cycle
                lf(i,j,k)     =     lf(i,j,k)+CoulCoefs(l+coulorder+1) *f(i,j,nz+k+l)
                lf(i,j,nz+1-k)=lf(i,j,nz+1-k)+CoulCoefs(-l+coulorder+1)*f(i,j,1-k-l)
              enddo
            enddo
        enddo
    enddo
#endif

    !---------------------------------------------------------------------------
    lf = lf/(dx**2)

  end function coulomb_laplacian

  subroutine clean_coulomb()
    use folding 
    
    if(allocated(r))                 deallocate(r)
    if(allocated(Gaussx))            deallocate(Gaussx)
    if(allocated(gaussy))            deallocate(gaussy)
    if(allocated(gaussz))            deallocate(gaussz)
    if(allocated(Coulcoefs))         deallocate(coulcoefs)
    if(allocated(coulmeshx))         deallocate(coulmeshx, coulmeshy, coulmeshz)
  end subroutine clean_coulomb

end module Coulombmod
