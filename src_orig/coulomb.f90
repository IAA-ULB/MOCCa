module Coulombmod
 !==============================================================================
 !  #######   ##   #    # #####   ##   #      #    #  ####
 !     #     #  #  ##   #   #    #  #  #      #    # #
 !     #    #    # # #  #   #   #    # #      #    #  ####
 !     #    ###### #  # #   #   ###### #      #    #      #
 !     #    #    # #   ##   #   #    # #      #    # #    #
 !     #    #    # #    #   #   #    # ######  ####   ####
 !
 !  Copyright W. Ryssens & M. Bender
 !
 !==============================================================================
 ! Module that solves the Coulomb problem of the proton density. 
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
 use derivatives, only: CoulombLaplacian
 use moments 
 use parameterization
 
 implicit none

 public
 
 !------------------------------------------------------------------------------
 ! The array containing the Coulomb Potential in the original box,
 ! enlarged with the boundary conditions. 
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable :: CoulombPotential(:,:,:)
 real(KIND=dp), allocatable :: ExchangePotential(:)
 real(KIND=dp), allocatable :: chargedensity(:,:,:)
 !------------------------------------------------------------------------------
 !Precision required of the Coulomb Solvers
 real(KIND=dp), public              :: Prec
 !------------------------------------------------------------------------------
 ! Number of boundary conditions to put on all sides of the box.
 ! Currently hard-coded to be 2, corresponding to the Laplacian defined in
 ! derivatives.f90.
 integer, parameter :: BC = 2
 !------------------------------------------------------------------------------
 ! Arrays containing,
 ! 1) the values of the spherical harmonics on the extended mesh and
 ! 2) the value of the radial coordinate r on the extended mesh.
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable, target :: SpherHarmCoulomb(:,:,:,:,:,:),r(:,:,:)
 !------------------------------------------------------------------------------
 ! Maximum l of the multipole moments to use in the boundary conditions
 ! Currently hardcoded at 8: does not cost anything CPU-time wise and 
 ! has been shown to be sufficient in MOCCa.
 integer, parameter :: maxm=8
 !------------------------------------------------------------------------------
 ! Gaussian matrices, to be used when folding is required.
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable :: Gaussx(:,:), Gaussy(:,:), Gaussz(:,:)
   
contains

 subroutine SolveCoulomb(rhop)
    !---------------------------------------------------------------------------
    ! Master routine to solve the Coulomb problem for a given source-density.
    ! Input is the point proton density.
    !---------------------------------------------------------------------------

    use Folding    

    real(KIND=dp), intent(in)  :: rhop(mv)
    real(KIND=dp), allocatable :: source(:,:,:)
    integer                    :: i,j,k
    real(KIND=dp)              :: R, dd(2)
    
    if(.not.allocated(Source)) then
        allocate(Source(nx+2,ny+2,nz+2))           ; Source = 0.0_dp
    endif

    !---------------------------------------------------------------------------
    ! Initialize all of the arrays.
    if(.not.allocated(CoulombPotential)) then
        call setupcoulomb
    endif
    
    if(coultreatment.eq.0) return
    
    do i=1,mv
      ChargeDensity(i,1,1) = rhop(i)
    enddo
    if(protonsize.ne.0.0_dp) then
        ! Fold the source with a Gaussian
        ChargeDensity=FoldGaussian(ChargeDensity, GaussX, GaussY, GaussZ, nx, ny, nz)
    endif
    
    !---------------------------------------------------------------------------
    ! Set up the source term: 
    ! For standard parameterizations it is the simply the proton density with
    ! a prefactor.
    Source = 0.0_dp
    do k=1,nz
        do j=1,ny
            do i=1,nx
                Source(i,j,k) = -4*pi*e2*Chargedensity(i,j,k)
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    ! Set the boundary conditions.
    call CoulombBound(Source)
 
    ! Solve for the direct coulomb potential   
    call ConjugGrad (CoulombPotential,Source,1,1,1,1000,.false.,prec)
 
    if(Coultreatment.eq.1) then
      ! Exchange potential in Slater approximation
      ExchangePotential = -(3.0/pi)**(1.0/3.0_dp)*e2*(rhop**(1.0_dp/3.0_dp))   
    else
      ExchangePotential = 0.0
    endif  
 end subroutine SolveCoulomb 

 subroutine SetupCoulomb
    !---------------------------------------------------------------------------
    ! Initialize the entire module. 
    !---------------------------------------------------------------------------
    use sphericalharmonics
    use folding
    
    integer       :: i,j,k
    
    !---------------------------------------------------------------------------
    ! Allocate the CoulombPotential array (second-order boundary conditions)
    allocate(CoulombPotential(nx+2,ny+2,nz+2)) ; CoulombPotential = 0.0_dp
    allocate(ExchangePotential(mv))            ; ExchangePotential = 0.0_dp
    allocate(ChargeDensity(nx,ny,nz))          ; ChargeDensity     = 0.0_dp
    !---------------------------------------------------------------------------
    ! Precision desired of the Coulomb solver
    Prec = 1.d-9/(dx**3*nx*ny*nz)
    
    !---------------------------------------------------------------------------
    ! Set-up the values of r and spherharmcoulomb on the mesh.
    allocate(r(nx+BC,ny+BC,nz+BC))                                 ;  r = 0.0_dp
    allocate(SpherHarmCoulomb(nx+BC,ny+BC,nz+BC,0:maxm,0:maxm,2)) 
    SpherHarmCoulomb = 0.0_dp
    
    do k=1,nz+BC
      do j=1,ny+BC
        do i=1,nx+BC
          r(i,j,k) = sqrt(coulmeshx(i)**2 + coulmeshx(j)**2 + coulmeshz(k)**2)
        enddo
      enddo
    enddo
    
    call GenSphericalHarmonics(maxm,nx+BC,ny+BC,nz+BC,coulmeshx,coulmeshy,     &
    &                 coulmeshz,SpherHarmCoulomb,QuantisationAxis,SecondaryAxis)
    !---------------------------------------------------------------------------
    ! If the proton has a finite size, we need to fold the density with a
    ! Gaussian. This sets up the necessary matrices.
    !---------------------------------------------------------------------------
    if(protonsize .ne. 0.0_dp) then
        if(.not.allocated(Gaussx)) then
            allocate(Gaussx(nx,nx), Gaussy(ny,ny), Gaussz(nz,nz)) 
            Gaussx = 0.0 ;  Gaussy = 0.0 ; Gaussz = 0.0
        endif
        !-----------------------------------------------------------------------
        ! Construct Gauss matrices
        call ConstructFoldingMatrices(Gaussx,Gaussy,Gaussz)
    endif
    
 end subroutine SetupCoulomb
    
 subroutine CoulombBound(source)
    !---------------------------------------------------------------------------
    ! Calculates the boundary conditions of the Coulomb potential based on the 
    ! multipole moments of the point charge density. 
    !---------------------------------------------------------------------------
    
    use folding

    real(KIND=dp), intent(in) :: source(mv)
    integer                   :: i,j,k,l,m, im
    real(KIND=dp)             :: Qlm
    type(Moment), pointer     :: Current
    logical                   :: cont

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
    CoulombPotential=0
    
    nullify(Current)
    Current => Root
    Cont = .true.
    do while(Cont)
      l = Current%l
      m = Current%m
      Im = 1
      if(Current%Impart) Im = 2

      ! Recalculate the multipole distribution, since source is not
      ! necessarily the point proton distribution.
      Qlm = -sum(Source * SpherHarmCoulomb(1:mv,1,1,l,m,Im)) * dv
      
      !  Previous implementation based on values of the multipole moments
      !Qlm = e2*Current%Value(2)*(4*pi/(2*l+1)) 
      
      do k=1,nz+BC
        do j=1,ny+BC
          do i=1,nx+BC
            if( (i.gt.nx) .or. (j.gt.ny) .or. (k.gt.nz)) then
              CoulombPotential(i,j,k) = CoulombPotential(i,j,k) +             &
              &           Qlm*SpherHarmCoulomb(i,j,k,l,m,Im)/(r(i,j,k)**(2*l+1))
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
    ! Construct the matrices for Gaussian folding, taking into account the 
    ! symmetries of the density. 
    !---------------------------------------------------------------------------
    use folding
    
    real(KIND=dp), intent(out) :: Gx(:,:), Gy(:,:), Gz(:,:)
    real(KIND=dp)              :: r0, D, X, Y,Z
    integer :: linX, linY, linZ, i,j,k
     
    linX = nx
    linY = ny
    linZ = nz

    r0 = protonsize * sqrt(2.0/3.0)
    !---------------------------------------------------------------------------
    ! Elements actually represented on the mesh
    do i=1,nx
        do j=1,nx           
            Gx(i,j) = Gaussian(meshx(i), meshx(j), r0)
        enddo
    enddo
    do i=1,ny
        do j=1,ny          
            Gy(i,j) = Gaussian(meshy(i), meshy(j), r0)
        enddo
    enddo
    do i=1,nz
        do j=1,nz         
            Gz(i,j) = Gaussian(meshz(i), meshz(j), r0)
        enddo
    enddo
    
    !---------------------------------------------------------------------------
    ! Elements to be gotten by symmetry
    do i=1,nx
        do j=1,nx
            Gx(i,j) = Gx(i,j) + Gaussian(-meshx(i), meshx(j),r0)
        enddo
    enddo
    do i=1,ny
        do j=1,ny
            Gy(i,j) = Gy(i,j) + Gaussian(-meshy(i), meshy(j),r0)
        enddo
    enddo
    do i=1,nz
        do j=1,nz
            Gz(i,j) = Gz(i,j) + Gaussian(-meshz(i), meshz(j),r0)
        enddo
    enddo
    !---------------------------------------------------------------------------
    ! Normalize; otherwise the integral of the Gaussian is not necessarily 1.
    Gx = Gx/(sum(Gx(:,1))*dx)
    Gy = Gy/(sum(Gy(:,1))*dx)
    Gz = Gz/(sum(Gz(:,1))*dx)

 end subroutine ConstructFoldingMatrices
 
 function CoulombEnergy_direct(rhop) result(CEnergy)
    !---------------------------------------------------------------------------
    ! Calculate the (direct) electrostatic energy of the system.
    !
    ! Note that rhop is not necessarily the point-proton density that is 
    ! passed in.
    !---------------------------------------------------------------------------
    real(KIND=dp) :: CEnergy
    real(KIND=dp), intent(in) :: rhop(mv)
    integer       :: i,j,k
    
    CEnergy = 0.0_dp
    do k=1,nz
        do j=1,ny
            do i=1,nx
                CEnergy = CEnergy + rhop(i + nx*(j-1) + ny*nx*(k-1)) *         &
                &                   CoulombPotential(i,j,k)
            enddo
        enddo
    enddo
    CEnergy = CEnergy * dv * 0.5_dp

 end function CoulombEnergy_Direct
 
 function CoulombEnergy_Exchange(rhop) result(CEnergy)
    !---------------------------------------------------------------------------
    ! Calculate the (exchange) electrostatic energy of the system in the 
    ! Slater approximation.
    ! Note that rhop is not necessarily the point-proton density that is 
    ! passed in.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: rhop(mv)
    real(KIND=dp) :: factor, Cenergy
    
    Cenergy = 0.0_dp
    if(coultreatment.ne.1) return

    factor  = -0.75_dp*(3/pi)**(1/3._dp)*e2*dv
    Cenergy = factor*sum(rhop**(4.0/3.0))
    
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
    
    integer                    :: iteration
    real(KIND=dp), allocatable :: p_k(:,:,:), Temp(:,:,:)
    real(KIND=dp), allocatable :: Residual(:,:,:)
    real(KIND=dp)              :: PoissonNorm, Integral, a_k, c_k
    real(KIND=dp)              :: NewPoissonNorm

    Residual = - CoulombLaplacian(Solution,sx,sy,sz)
    Residual = Residual + SourceTerm
    !---------------------------------------------------------------------------
    !The variable p_k is the conjugate direction. It starts out equal to our
    ! initial Residual.
    p_k=Residual
    !Poissonnorm is r_(k+1)^T r_(k+1)
    PoissonNorm = sum(Residual**2)

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
      if(sum(Residual**2).le.Precis) then
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
!      print *, 'Coul, it',  iteration, PoissonNorm
    enddo

    return
  end subroutine ConjugGrad

end module Coulombmod
