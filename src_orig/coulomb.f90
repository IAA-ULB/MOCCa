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
 use moments 
 use parameterization
 use timing 

 implicit none

 public
 
 !------------------------------------------------------------------------------
 ! The array containing the Coulomb Potential in the full Coulomb box
 real(KIND=dp), allocatable :: CoulombPotential(:,:,:)
 real(KIND=dp), allocatable :: ExchangePotential(:,:,:)
 ! Array with the folded Coulomb potential, necessary if we take the 
 ! finite size of the nucleons into account
 real(KIND=dp), allocatable :: FoldedCoul(:,:,:,:), FoldedExchange(:,:,:,:)
 !------------------------------------------------------------------------------
 !Precision required of the Coulomb Solvers
 real(KIND=dp)              :: Prec = 1e-10
 !------------------------------------------------------------------------------
 ! Number of boundary conditions to put on all sides of the box.
 integer :: BC = 2
 !------------------------------------------------------------------------------
 ! Arrays containing,
 ! 1) the values of the spherical harmonics on the extended mesh and
 ! 2) the value of the radial coordinate r on the extended mesh.
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable, target :: SpherHarmCoulomb(:,:,:,:,:,:),r(:,:,:)  
 !------------------------------------------------------------------------------
 ! Coordinates of the mesh in the enlarged coulomb box.
 real(KIND=dp), allocatable :: coulmeshx(:), coulmeshy(:), coulmeshz(:)
 real(KIND=dp), allocatable :: coulgrid(:,:)
 !------------------------------------------------------------------------------
 ! Maximum l of the multipole moments to use in the boundary conditions
 ! Currently hardcoded at 8: does not cost anything CPU-time wise and 
 ! has been shown to be sufficient in MOCCa.
 integer, parameter :: maxm=8
 !------------------------------------------------------------------------------
 ! Gaussian matrices, to be used when folding of the nucleon densities to 
 ! obtain the charge densities are required. The fourth index is the isospin 
 ! index, the third whether it is the Gaussian with positive or negative sign.
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable :: Gaussx(:,:,:,:), Gaussy(:,:,:,:), Gaussz(:,:,:,:)
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

contains

 subroutine SolveCoulomb(R, F, guess)
    !---------------------------------------------------------------------------
    ! Master routine to solve the Coulomb problem for a given source-density.
    !
    ! Input : 
    !     R : density vector with a precalculated charge density
    !     initial_guess :  optional initial guess for the Coulomb potential
    ! Output:
    !     F : potential vector; only the Coulomb fields are modified
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Note: this routine does not calculate the charge density based on the 
    !       the proton and neutron densities. This is assumed to have been 
    !       done before.
    !---------------------------------------------------------------------------

    use Folding    

    type(DensityVector), intent(in)      :: R
    type(PotentialVector), intent(inout) :: F
    real(KIND=dp), intent(in), optional  :: guess(:,:,:)
    real(KIND=dp), allocatable           :: initial_guess(:,:,:)

    real(KIND=dp), allocatable      :: source(:,:,:)
    integer                         :: i,j,k,ii, ox, oy, oz
    
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
     !NS: no offset for pasta
#if(USE_Periodic==1)
     BC=0
#endif
    endif

    ! Determine the offsets of the original mesh inside the larger Coulomb mesh    
    coul_offset_x = BC ; coul_offset_Y = BC ; coul_offset_z = BC
    
    
    ! If any given axis is not represented, the offset of the mesh in that
    ! direction is zero.
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
    ! Initialize all of the arrays.
    if(.not.allocated(SpherHarmCoulomb)) then
        ! Check for all things that should have been setup
        ! Note: checking for the allocation of CoulombPotential would be more
        !       natural, but it is possible that the potentials have been read
        !       from file while nothing else in this module has been set up
        if(present(guess)) then
          initial_guess = guess
        else
          allocate(initial_guess(nx+BC+coul_offset_x, &
          &                      ny+BC+coul_offset_y, & 
          &                      nz+BC+coul_offset_z))
          initial_guess = 0.0d0
        endif
        call setupcoulomb(F, initial_guess)
    endif
      if(coultreatment.eq.0) then
       call stop_timer(T_coulomb)
       return
    endif        
    !---------------------------------------------------------------------------
    ! Set up the source term: - 4 * pi * charge_density
    ! Note that this is set up in the middle of the box, i.e. no source density
    ! at the edges of the Coulomb box
    if(.not.allocated(Source)) then
        allocate(Source(nx+BC+coul_offset_x, &
        &               ny+BC+coul_offset_y, &
        &               nz+BC+coul_offset_z))           
    endif
    
    Source = 0.0_dp
    do k=1,nz
      do j=1,ny
        do i=1,nx
          Source(i+coul_offset_x,j+coul_offset_y,k+coul_offset_z) = &
          &                                  -4*pi*e2*R%Chargedensity(i,j,k)
        enddo
      enddo
    enddo
    !---------------------------------------------------------------------------
    ! Set the boundary conditions.
    call CoulombBound(Source, F)
    
    ! Solve for the direct coulomb potential  
    ! Note that the symmetry properties (+1,+1,+1) are never changed: 
    ! Hephaestos modifies directly the Coulomb_Laplacian routine when necessary 
    call ConjugGrad (F%CoulombPotential,Source, +1,+1,+1,1000,.false.,prec)

    if(Coultreatment.eq.1) then
      !-------------------------------------------------------------------------
      ! Exchange potential in Slater approximation
      if( all(protonsize.eq.0.0) .and. all(neutronsize.eq.0.0) ) then
        do k=1,nz
          do j=1,ny
              do i=1,nx
               ii = meshindex(i,j,k)
               F%ExchangePotential(i,j,k) =                                    & 
                 &   -(3.0/pi)**(1.0/3.0_dp)*e2*(R%D_I_I(ii,2)**(1.0_dp/3.0_dp))   
              enddo
          enddo
        enddo
      else
       ! Incorporating finite size effects
       F%ExchangePotential =  &
       &           -(3.0/pi)**(1.0/3.0_dp)*e2*(R%ChargeDensity**(1.0_dp/3.0_dp)) 
     endif
    else
      !-------------------------------------------------------------------------
      ! No Coulomb Exchange
      F%ExchangePotential = 0.0
    endif  
    
    ! This routine calculates FoldedCoul and FoldedExchange if required
    call Obtain_folded_potentials()

    call stop_timer(T_coulomb)
    deallocate(source)
 end subroutine SolveCoulomb 

 subroutine Obtain_folded_potentials()
    !---------------------------------------------------------------------------
    ! Obtain the folded Coulomb potentials (direct and exchange) if needed.
    ! It is a separate routine from SolveCoulomb because it should also be 
    ! callable from the routine to read potentials.
    !---------------------------------------------------------------------------

    if(any(protonsize .ne. 0.0_dp) .or. any(neutronsize.ne.0.0_dp)) then
      if(nucleonsize_selfconsistent) then
         FoldedCoul    =FoldCoulombPotential(CoulombPotential(                 &
         &                                    coul_offset_x+1:coul_offset_x+nx,&
         &                                    coul_offset_y+1:coul_offset_y+ny,&
         &                                    coul_offset_z+1:coul_offset_z+nz))
         
         if(coultreatment .eq. 1) then
           FoldedExchange=FoldCoulombPotential(ExchangePotential(              &
           &                                  coul_offset_x+1:coul_offset_x+nx,&
           &                                  coul_offset_y+1:coul_offset_y+ny,&
           &                                  coul_offset_z+1:coul_offset_z+nz))
         else
           if(.not. allocated(FoldedExchange)) then
            allocate (FoldedExchange(nx,ny,nz,2))
           endif
           foldedexchange = 0.0d0
         endif
      endif
    endif
 
 end subroutine Obtain_folded_potentials

 subroutine ConstructChargeDensity(R)
    !---------------------------------------------------------------------------
    ! Construct the charge density from the proton and neutron densities,  
    ! using various effective form
    !---------------------------------------------------------------------------

    use Folding

    type(DensityVector),intent(inout) :: R
    !real(KIND=dp), allocatable :: rho_charge(:,:,:)
    real(KIND=dp)              :: temp(nx,ny,nz)
    integer                    :: i,j,k
    real(KIND=dp)              :: rho_el, linX, linY, linZ, volume
    !real(KIND=dp)              :: rrr,xx,yy,zz
    
    call start_timer(T_chargedensity)
    
    ! Deallocation that rho_charge does not have the wrong dimensions
    if(allocated(R%chargedensity)) deallocate(R%chargedensity)
    if(.not.allocated(R%chargedensity)) then
        allocate(R%chargedensity(nx,ny,nz))
    endif
    R%chargedensity = 0.0
    !---------------------------------------------------------------------------
    ! If we account for the finite extent of the charge of the nucleus, then 
    ! we need to fold densities and potentials with gaussians. This sets up the 
    ! required matrices.
    !
    ! Note: this little piece of code is duplicated, since in different 
    !       runmodes of the code different Coulomb routines get called in 
    !       different order; this makes sure we get no segfaults.
    !---------------------------------------------------------------------------
    if(any(protonsize .ne. 0.0_dp) .or. any(neutronsize.ne.0.0_dp)) then
        if(.not.allocated(Gaussx)) then
            allocate(Gaussx(nx,nx,2,2), Gaussy(ny,ny,2,2), Gaussz(nz,nz,2,2)) 
            Gaussx = 0.0 ;  Gaussy = 0.0 ; Gaussz = 0.0
        endif
        !-----------------------------------------------------------------------
        ! Construct Gauss matrices
        call ConstructFoldingMatrices(Gaussx,Gaussy,Gaussz)
    endif

    !---------------------------------------------------------------------------
    ! Proton contributions to the charge density.    
    ! We start from the proton point density
    do k=1,nz
      do j=1,ny
        do i=1,nx
            temp(i,j,k) = R%D_I_I(meshindex(i,j,k),2)
        enddo       
      enddo
    enddo
    
    if(protonsize(1).gt.0.0) then
        ! Fold the source with a Gaussian
        R%chargedensity = R%chargedensity + &
        & FoldGaussian(temp, GaussX(:,:,1,2), GaussY(:,:,1,2), GaussZ(:,:,1,2),&
        &                                                            nx, ny, nz)
    endif
    if(protonsize(2).gt.0.0) then
        ! Fold the source with another Gaussian, this time with minus sign.
        R%chargedensity = R%chargedensity + &
        & FoldGaussian(temp, GaussX(:,:,2,2), GaussY(:,:,2,2), GaussZ(:,:,2,2),&
        &                                                            nx, ny, nz)
    endif

    !do k=1,nz
    !  do j=1,ny
    !    do i=1,nx
    !        xx=dx/2.d0+(i-1)*dx
    !        yy=dx/2.d0+(j-1)*dx
    !        zz=dx/2.d0+(k-1)*dx
    !        rrr=sqrt(xx**2 + yy**2 + zz**2)
    !        temp(i,j,k)=0.070679973d0/(1.d0+exp((rrr-4.d0)/0.3d0))
    !    enddo       
    !  enddo
    !enddo

    if(all(protonsize.eq.0.0)) then
        R%chargedensity = temp
    endif
   
        
    !---------------------------------------------------------------------------
    ! Neutron contributions to the charge density.
    if(all(neutronsize.eq.0.0)) then
    !NS: subtruct electron background
#if(PASTA==1)
      linX=nx
      linY=ny
      linZ=nz
$REDUX     linX=2*nx
$REDUY     linY=2*ny
$REDUZ     linZ=2*nz
      volume=linx*dx*linY*dx*linZ*dx
      rho_el=protons/volume
      rho_charge=rho_charge-rho_el
      !print *,'rho_el=',rho_el
#endif
      call stop_timer(T_chargedensity)
      return
    endif

    do k=1,nz
      do j=1,ny
        do i=1,nx
           temp(i,j,k) = R%D_I_I(meshindex(i,j,k),1)
        enddo
      enddo
    enddo
          
    if(neutronsize(1).gt.0.0) then
        ! Fold the source with a Gaussian
        R%chargedensity = R%chargedensity + & 
        & FoldGaussian(temp, GaussX(:,:,1,1), GaussY(:,:,1,1), GaussZ(:,:,1,1),&
        &                                                            nx, ny, nz)
    endif
    if(neutronsize(2).gt.0.0) then
        ! Fold the source with a Gaussian, minus sign this time
        R%chargedensity = R%chargedensity + & 
        & FoldGaussian(temp, GaussX(:,:,2,1), GaussY(:,:,2,1), GaussZ(:,:,2,1),&
        &                                                            nx, ny, nz)
    endif
    !---------------------------------------------------------------------------
    !NS: subtract electron background
#if(PASTA==1)
    
    linX=nx
    linY=ny
    linZ=nz

$REDUX   linX=2*nx
$REDUY   linY=2*ny
$REDUZ   linZ=2*nz

    volume=linx*dx*linY*dx*linZ*dx
    rho_el=protons/volume
    !print *, protons,sum(rho_charge)*dv
    R%chargedensity = R%chargedensity -rho_el
#endif
    call stop_timer(T_chargedensity)
 end subroutine ConstructChargeDensity

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

 subroutine SetupCoulomb(F,initialguess)
    !---------------------------------------------------------------------------
    ! Initialize the entire module and the fields in the potentialvector
    !
    ! Input:
    !    F : potentialvector to be initialized
    !---------------------------------------------------------------------------
    use sphericalharmonics
    use folding
    use vectors
    
    type(PotentialVector), intent(inout) :: F
    real(KIND=dp), intent(in)            :: initialguess(:,:,:)
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
    allocate(CoulombPotential(ox,oy,oz))
    CoulombPotential = 0.0_dp
    ! The exchange potential is only defined on the original mesh
    allocate(ExchangePotential(nx,ny,nz)) 
    ExchangePotential = 0.0_dp
    !---------------------------------------------------------------------------
    ! Precision desired of the Coulomb solver
    Prec = 1.d-12/(dx**3)

    !---------------------------------------------------------------------------
    ! Set-up the values of r and spherharmcoulomb on the Coulomb mesh.
    call inimesh(coulmeshx,coulmeshy,coulmeshz,nx+BC+coul_offset_x, &
    &                                          ny+BC+coul_offset_y, &
    &                                          nz+BC+coul_offset_z, &
    &                                          coulgrid,0.0d0,0.0d0,0.0d0)

    allocate(r(ox,oy,oz))  ;  r = 0.0_dp
    allocate(SpherHarmCoulomb(ox,oy,oz,0:maxm,0:maxm,2)) 
    SpherHarmCoulomb = 0.0_dp
    
    do k=1,oz
      do j=1,oy
        do i=1,ox
          r(i,j,k) = sqrt(coulmeshx(i)**2 + coulmeshy(j)**2 + coulmeshz(k)**2)
        enddo
      enddo
    enddo
    
    call GenSphericalHarmonics(maxm,ox,oy,oz,                                  &
    &                          coulmeshx,coulmeshy, coulmeshz,SpherHarmCoulomb,&
    &                          QuantisationAxis,SecondaryAxis)
    
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
        call ConstructFoldingMatrices(Gaussx,Gaussy,Gaussz)
    endif
    
 end subroutine SetupCoulomb
    
 subroutine CoulombBound(source, F)
    !---------------------------------------------------------------------------
    ! Calculates the boundary conditions of the Coulomb potential based on the 
    ! multipole moments of the point charge density. 
    !
    ! Input:
    !   source : source charge density to compute boundary conditions for
    ! Output:
    !   F      : potentialvector with boundary conditions applied to the 
    !            CoulombPotential.
    !---------------------------------------------------------------------------
    
    use folding
    use vectors

    real(KIND=dp), intent(in) :: source(:,:,:)
    type(PotentialVector), intent(inout) :: F
    real(KIND=dp)             :: fac
    integer                   :: i,j,k,l,m, im, ox, oy, oz
    real(KIND=dp)             :: Qlm
    type(Moment), pointer     :: Current
    logical                   :: cont, condition

    !---------------------------------------------------------------------------
    ! Calculate the multipole moment expansion of the source term.
    ! We put the boundary condition on every point, and use the potential 
    ! generated this way as an initial guess.
    !---------------------------------------------------------------------------
    ! The source density is expanded into multipole moments for the boundary   
    ! conditions. This is linked to the linked list of multipole moments, 
    ! not because they are calculated with them, but simply to not have
    ! another place in the code where decisions regarding symmetries need
    ! to be chosen.
    !---------------------------------------------------------------------------

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
              F%CoulombPotential(i,j,k) = 0.0d0
          endif
        enddo
      enddo
    enddo
    
    nullify(Current)
    Current => Root
    Cont = .true.
    do while(Cont)
      l = Current%l
      m = Current%m
      Im = 1
      if(Current%Impart) Im = 2

      !------------------------------------------------------------------------
      ! Recalculate the multipole distribution, since source is not
      ! necessarily the point proton distribution.
      !-------------------------------------------------------------------------
      ! The source density is expanded into multipole moments for the boundary   
      ! conditions. This is linked to the linked list of multipole moments, 
      ! not because they are calculated with them, but simply to not have
      ! another place in the code where decisions regarding symmetries need
      ! to be chosen.
      !-------------------------------------------------------------------------
      Qlm = 0
      do k=1,oz
        do j=1,oy
          do i=1,ox
            Qlm = Qlm - Source(i,j,k) * SpherHarmCoulomb(i,j,k,l,m,Im)
          enddo
        enddo
      enddo 
   
      Qlm = Qlm * dv/(2*l+1)    

      
      ! WR 09/02/22 
      ! Bugfix: the code only calculates Q_lm for positive m, but the 
      !         complex conjugate multipole moments Q_l(-m) should contribute
      !         as well. This means that (1) real multipole moments contribute
      !         with a factor 2 and (2) imaginary multipole moments do not 
      !         contribute at all.
      fac = 1
      if(m.ne.0) fac = 2    

      !  Previous implementation based on values of the multipole moments
      !Qlm = e2*Current%Value(2)*(4*pi/(2*l+1)) 
      do k=1,oz
        do j=1,oy
          do i=1,ox

            condition = .false.
 $REDUX     if(i.gt.nx) condition =.true.
 $REDUY     if(j.gt.ny) condition =.true.
 $REDUZ     if(k.gt.nz) condition =.true.
          
 $FULLX     if(i.le.BC)    condition =.true.         
 $FULLX     if(i.gt.nx+BC) condition =.true.         

 $FULLY     if(j.le.BC)    condition =.true.         
 $FULLY     if(j.gt.ny+BC) condition =.true.         
          
 $FULLZ     if(k.le.BC)    condition =.true.         
 $FULLZ     if(k.gt.nz+BC) condition =.true.         
          
            
            if(condition) then
              F%CoulombPotential(i,j,k) = F%CoulombPotential(i,j,k) +          &
              &       fac*Qlm*SpherHarmCoulomb(i,j,k,l,m,Im)/(r(i,j,k)**(2*l+1))
            endif
          enddo
        enddo
      enddo
      !-------------------------------------------------------------------------
      !Transferring to the next moment in the list, until the r**2 is reached or
      ! the highest admissible L.
      if(Current%Next%l .ge. 0 .and. Current%Next%l .le. maxm) then
        Current => Current%Next
      else
      !Signalling that there is no further moment
        Cont=.false.
      endif
    end do
    nullify(current)
 end subroutine CoulombBound

 subroutine ConstructFoldingMatrices(Gx,Gy,Gz)  
    !---------------------------------------------------------------------------
    ! Construct the matrices for Gaussian folding.
    !---------------------------------------------------------------------------
    use folding
    
    real(KIND=dp), intent(out) :: Gx(:,:,:,:), Gy(:,:,:,:), Gz(:,:,:,:)
    real(KIND=dp)              :: rplus(2), rmin(2)
    real(KIND=dp)              :: hbom, mhb, B
    integer                    :: it

    ! The determination from input for neutrons and protons is not the same     
    rplus(1) = sqrt(neutronsize(1))
    rmin(1)  = sqrt(neutronsize(2))

    rplus(2) = protonsize(1) * sqrt(2.0/3.0)
    rmin(2)  = protonsize(2) * sqrt(2.0/3.0)

    !---------------------------------------------------------------------------
    ! Harmonic-oscillator correction
    if(hocomform) then
        ! hbar x omega
        hbom  = 41.0 * (neutrons + protons)**(-1.0/3.0)
        ! 2m/hbar^2 
        mhb = 2.0/(1.0/hbm(1)+1.0/hbm(2)) 
        ! B^{-1} = hbar * omega/m * A = 1/2 * A * hbar omega * 2m/hbar^2 
        B = sqrt( 1.0/( 0.5 * hbom/mhb  * (neutrons + protons)))

        do it=1,2
            if(rplus(it).ne.0.0) then
                rplus(it) = sqrt(rplus(it)**2 - B**2)
            endif
            if(rmin(it).ne.0.0) then
                rmin(it) = sqrt(rmin(it)**2 - B**2)
            endif
        enddo
    endif

    do it=1,2
      if(rplus(it) .ne. 0.0_dp) then 
        call Gauss_1D(Gx(:,:,1,it), meshx, nx, rplus(it), sx_rho)
        call Gauss_1D(Gy(:,:,1,it), meshy, ny, rplus(it), sy_rho)
        call Gauss_1D(Gz(:,:,1,it), meshz, nz, rplus(it), sz_rho)
      endif
      if(rmin(it) .ne. 0.0_dp) then 
        call Gauss_1D(Gx(:,:,2,it), meshx, nx, rmin(it),  sx_rho)
        call Gauss_1D(Gy(:,:,2,it), meshy, ny, rmin(it),  sy_rho)
        call Gauss_1D(Gz(:,:,2,it), meshz, nz, rmin(it),  sz_rho)
      endif
    enddo
    
 end subroutine ConstructFoldingMatrices

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
    ! Calculate the (exchange) electrostatic energy of the system in the 
    ! Slater approximation. Note that rhop is not necessarily the point-proton 
    ! density that is passed in.
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

    Residual = - CoulombLaplacian(Solution,sx,sy,sz)   
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
      Temp = CoulombLaplacian(p_k,sx,sy,sz)

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
      !print *, 'Coul, it',  iteration, PoissonNorm
    enddo
    
    return
  end subroutine ConjugGrad

  function Coulomblaplacian(f, sx, sy, sz) result(lf)
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

  end function Coulomblaplacian

  subroutine clean_coulomb()
    if(allocated(CoulombPotential))  deallocate(CoulombPotential)
    if(allocated(ExchangePotential)) deallocate(ExchangePotential)
    if(allocated(FoldedCoul))        deallocate(FoldedCoul)
    if(allocated(FoldedExchange))    deallocate(FoldedExchange)
    if(allocated(SpherHarmCoulomb))  deallocate(SpherHarmCoulomb)
    if(allocated(r))                 deallocate(r)
    if(allocated(Gaussx))            deallocate(Gaussx)
    if(allocated(gaussy))            deallocate(gaussy)
    if(allocated(gaussz))            deallocate(gaussz)
    if(allocated(Coulcoefs))         deallocate(coulcoefs)
    if(allocated(coulmeshx))         deallocate(coulmeshx, coulmeshy, coulmeshz)

  end subroutine clean_coulomb

end module Coulombmod
