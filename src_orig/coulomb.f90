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
 !------------------------------------------------------------------------------
 !Precision required of the Coulomb Solvers
 real(KIND=dp), public              :: Prec
 
 integer       :: CoulombSolver = 0
    

 !------------------------------------------------------------------------------
 ! Number of boundary conditions to put on all sides of the box.
 integer :: BCNumber = 2

 real(KIND=dp) :: e2 =1.43996446_dp

contains

 subroutine SolveCoulomb
    !
    !
    !
    
    if(.not.allocated(CoulombPotential)) then
        call setupcoulomb
    endif
    
    call CoulombBound
    
     ! Calculate the values of the radial coordinate at every mesh point
    if(CoulombSolver .ne. 0 )then
        call ConjugGrad (CoulombPotential,Source,1,1,1,1000,.false.,prec)
    endif
 end subroutine SolveCoulomb 

 subroutine SetupCoulomb
    !---------------------------------------------------------------------------
    ! Initialize the entire module. 
    !---------------------------------------------------------------------------
    
    ! Allocate the CoulombPotential matrix
    allocate(CoulombPotential(nx+2,ny+2,nz+2))
    CoulombPotential = 0.0_dp
    allocate(Source(nx+2,ny+2,nz+2)) ;
    Source = 0.0_dp
    
    ! Precision desired of the Coulomb solver
    Prec = 1.d-9/(dx**3*nx*ny*nz)
    
 end subroutine SetupCoulomb
    
 subroutine CoulombBound
    !---------------------------------------------------------------------------
    ! Calculates the boundary conditions of the Coulomb potential based on the 
    ! multipole moments of the point charge density. 
    !
    ! Currently only takes into account the monopole and quadrupole term in EV8
    ! like calculations. 
    !---------------------------------------------------------------------------
 
    integer :: i,j,k   
    real(KIND=dp) :: r
    real(KIND=dp) :: Q20, Q22, x,y,z
    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                Source(i,j,k) = -  D_I_I(i + (j-1)*nx + (k-1)*ny*nx,2)
            enddo
        enddo
    enddo
    
    !---------------------------------------------------------------------------
    ! Temporary: compute the quadrupole moments of the density
    !---------------------------------------------------------------------------
    Q20 = 0.0
    Q22 = 0.0
    do k=1,nz
        z = dx/2.0 + (k-1) * dx
        do j=1,ny
            y = dx/2.0 + (j-1) * dx
            do i=1,nx
                x = dx/2.0 + (i-1) * dx
                Q20 = Q20 + D_I_I(i+(j-1)*nx+(k-1)*ny*nx,2)*(2 * z**2 - x**2 - y**2)
                Q22 = Q22 + D_I_I(i+(j-1)*nx+(k-1)*ny*nx,2)*(           x**2 - y**2)
            enddo
        enddo
    enddo
    Q20 = Q20 * dv /( 4 * pi/5.0)
    Q22 = Q22 * dv /( 4 * pi/5.0)
    !print *, Q20, Q22
    !---------------------------------------------------------------------------
    
    do k=1,nz+2
        do j=1,ny+2
            do i=1,nx+2
                if( (i .le. nx) .and. (j .le. ny) .and. (k .le. nz)) cycle
                
                r = dx * sqrt( (i-0.5)**2 + (j-0.5)**2 + (k-0.5)**2)
                ! dx/2.0_dp + dx*(i-1)
                CoulombPotential(i,j,k) = e2 * (   Protons/r) !                   &
                !&                                   + Q20/r**3 + Q22/r**3)
            enddo
        enddo
    enddo
   
 end subroutine CoulombBound
 
 function CalcCoulombEnergy() result(CEnergy)
    
    real(KIND=dp) :: CEnergy
    integer :: i,j,k
    
    CEnergy = 0.0_dp
    do k=1,nz
        do j=1,ny
            do i=1,nx
                CEnergy = CEnergy + D_I_I(i + nx*(j-1) + ny*nx*(k-1),2) *      &
                &                   CoulombPotential(i,j,k)
            enddo
        enddo
    enddo
    do i=1,mv
        
    enddo
    CEnergy = CEnergy * dv * 0.5_dp

 end function CalcCoulombEnergy
 
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
    Residual = Residual + 4 * pi * e2 * SourceTerm
    
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
    enddo

    return
  end subroutine ConjugGrad
end module Coulomb
