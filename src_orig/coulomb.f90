module Coulomb
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
 ! Module that solves the Coulomb problem of the proton point density.
 !
 !
 !==============================================================================

 use geninfo
 use densities
 use derivatives, only: CoulombLaplacian
 
 implicit none

 !------------------------------------------------------------------------------
 ! The array containing the Coulomb Potential in the original box,
 ! enlarged with the boundary conditions. 
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable :: CoulombPotential(:,:,:), source(:,:,:)
 real(KIND=dp), allocatable :: ExchangePotential(:)
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
 ! Value of the electron charge, squared
 real(KIND=dp) :: e2 =1.43996446_dp

contains

 subroutine SolveCoulomb(rhop)
    !---------------------------------------------------------------------------
    ! Master routine to solve the Coulomb problem for a given source-density.
    ! Input is the source on the mesh (modulo a factor 4*e2*pi).
    ! 
    ! 1) Puts everything in a box with 2 more points
    ! 2) 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: rhop(mv)
    integer :: i,j,k
    
    !---------------------------------------------------------------------------
    ! Initialize all of the arrays.
    if(.not.allocated(CoulombPotential)) then
        call setupcoulomb
    endif
    !---------------------------------------------------------------------------
    ! Set the boundary conditions.
    call CoulombBound(rhop)
 
    ! Solve for the direct coulomb potential   
    call ConjugGrad (CoulombPotential,Source,1,1,1,1000,.false.,prec)
 
    ! Exchange potential in Slater approximation
    ExchangePotential = -(3.0/pi)**(1.0/3.0_dp)*e2*(rhop**(1.0_dp/3.0_dp))     
 end subroutine SolveCoulomb 

 subroutine SetupCoulomb
    !---------------------------------------------------------------------------
    ! Initialize the entire module. 
    !---------------------------------------------------------------------------
    use sphericalharmonics
    
    integer       :: i,j,k
    real(KIND=dp) :: mesh(3,nx+BC, ny+BC, nz+BC)
    real(KIND=dp), pointer     :: sperharm(:,:,:,:,:,:)
    
    !---------------------------------------------------------------------------
    ! Allocate the CoulombPotential array (second-order boundary conditions)
    allocate(CoulombPotential(nx+2,ny+2,nz+2)) ; CoulombPotential = 0.0_dp
    allocate(Source(nx+2,ny+2,nz+2))           ; Source = 0.0_dp
    allocate(ExchangePotential(mv))            ; ExchangePotential = 0.0_dp
    !---------------------------------------------------------------------------
    ! Precision desired of the Coulomb solver
    Prec = 1.d-9/(dx**3*nx*ny*nz)
    
    !---------------------------------------------------------------------------
    ! Set-up the values of r and spherharmcoulomb on the mesh.
    ! Currently ONLY for EV8-like boxes.
    allocate(r(nx+BC,ny+BC,nz+BC))                                 ;  r = 0.0_dp
    allocate(SpherHarmCoulomb(nx+BC,ny+BC,nz+BC,0:maxm,0:maxm,2)) 
    SpherHarmCoulomb = 0.0_dp
    
    do i=1,nx+BC
          mesh(1,i,:,:) = (1/2.0_dp +(i-1))*dx
    enddo
    do j=1,ny+BC
          mesh(2,:,j,:) = (1/2.0_dp +(j-1))*dx
    enddo
    do k=1,nz+BC
          mesh(3,:,:,k) = (1/2.0_dp +(k-1))*dx
    enddo
    do k=1,nz+BC
      do j=1,ny+BC
        do i=1,nx+BC
          r(i,j,k) = sqrt(sum(mesh(:,i,j,k)**2))
        enddo
      enddo
    enddo
    
    call GenSphericalHarmonics(maxm,nx+BC,ny+BC,nz+BC,mesh,SpherHarmCoulomb,3,1)
    !---------------------------------------------------------------------------    
 end subroutine SetupCoulomb
    
 subroutine CoulombBound(rhop)
    !---------------------------------------------------------------------------
    ! Calculates the boundary conditions of the Coulomb potential based on the 
    ! multipole moments of the point charge density. 
    !
    ! Currently only takes into account the monopole and quadrupole term in EV8
    ! like calculations. 
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in) :: rhop(mv)
    integer                   :: i,j,k,l,m, ex
    real(KIND=dp)             :: factor, Qlm
    
    !---------------------------------------------------------------------------
    ! Set up the source term: 
    ! In mean-field calculations it is the proton density with a prefactor.
    Source = 0.0_dp
    do k=1,nz
        do j=1,ny
            do i=1,nx
                Source(i,j,k) = -  4*pi*e2*rhop(i + (j-1)*nx + (k-1)*ny*nx)
            enddo
        enddo
    enddo

    !---------------------------------------------------------------------------
    ! Calculate the multipole moment expansion of the source term.
    ! We put the boundary condition on every point, and use the potential 
    ! generated this way as an initial guess.
    !---------------------------------------------------------------------------
    ! The selection of multipole moments is currently not at all done according
    ! to the symmetries of the problem, but rather for an EV8 box.
    !
    ! Meaning: even l, even m and only real parts.
    !---------------------------------------------------------------------------
    CoulombPotential=0
    do l=0,maxm,2
      do m=0,l,2
        ! Real part
        Qlm = - sum(Source(1:nx,1:ny,1:nz) * SpherHarmCoulomb(:,:,:,l,m,1))*dv
        ! Imaginary part assumed to be zero for the moment.
        
        
        print *, 'lm', l, m, Qlm/e2/sqrt(4*pi), Qlm/e2/(4*pi)
        
        Qlm = 1.0/(2*l+1) * Qlm
        
        CoulombPotential = CoulombPotential + Qlm*SpherHarmCoulomb(:,:,:,l,m,1)&
        &                                   /(r**(2*l+1))
      enddo
    enddo
   
 end subroutine CoulombBound
 
 function CoulombEnergy_direct(rhop) result(CEnergy)
    !---------------------------------------------------------------------------
    ! Calculate the (direct) electrostatic energy of the system.
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
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: rhop(mv)
    real(KIND=dp) :: factor, Cenergy

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
    real(KIND=dp), allocatable :: p_k(:,:,:), Temp(:,:,:),Temp2(:,:,:)
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
      !print *, iteration, PoissonNorm
    enddo

    return
  end subroutine ConjugGrad
end module Coulomb
