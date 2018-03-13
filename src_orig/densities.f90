module densities
!===============================================================================
!  #######   ##   #    # #####   ##   #      #    #  ####
!     #     #  #  ##   #   #    #  #  #      #    # #
!     #    #    # # #  #   #   #    # #      #    #  ####
!     #    ###### #  # #   #   ###### #      #    #      #
!     #    #    # #   ##   #   #    # #      #    # #    #
!     #    #    # #    #   #   #    # ######  ####   ####
!
!  Copyright W. Ryssens & M. Bender
!
!===============================================================================
!
! Some technical notes:
!
! a) The densities are represented as vectors on the mesh, instead of 3D arrays. 
!    Advantages: 
!      *) Two less indices, which is needed to naively store densities
!         and their derivatives for N2/3LO. Fortran 90 was limited to 
!         arrays of rank 7, while Fortran 2008 lifted this limit to 15
!         many compilers (i.e. gfortran) do not yet support it.
!      *) Less involved coding in the functional.f90 file, since there
!         are less ':' to be put in the sums, and all of the spatial
!         indices are combined into one.
!      *) It MIGHT result in faster summing of densities, but this has never 
!         been a problem on the mean-field level.
!
!    Disadvantages: 
!      *) This needed a wrapper routine for the derivative functions Derive_grad
!         and derive_lap. Currently everything is accomplished by defining 
!         3D pointers, which is hopefully more efficient than a call to RESHAPE.
!         If this ever takes up a significant fraction of computation time, one
!         can think of explicitly writing the 1D derivative code.
!
! b) Note that the storage scheme for derivatives is not yet implemented on the 
!    level of densities, only on the level of derivatives of densities.
!    Thus
!           D_N_N is fully stored with indices (nx*ny*nz,3,3,2)
!    But 
!           Der_Der_D_I_I is stored as (nx*ny*nz,7,2)
!===============================================================================
use compilation
use geninfo
use wavefunctions, only: HFPsi, HFdPsi, HFddPsi, HFdddpsi
use wavefunctions, only: occupations, HFBlocks, blocks
use wavefunctions, only: upairing, vpairing
use derivatives 
use preconditioning 

implicit none

    !---------------------------------------------------------------------------
    ! Type declaration of the various densities
$DECLARATION   

    !---------------------------------------------------------------------------
    ! Density-mixing parameter default value
    real(KIND=dp) :: denmix = 0.75_dp
    
    !---------------------------------------------------------------------------
    ! Previous value(s) of the density rho.
    real(KIND=dp), allocatable :: D_I_I_hist(:,:,:)
    !---------------------------------------------------------------------------
    ! The amount of iterations to keep in memory for the density mixing
    integer           :: memory = 1
    !---------------------------------------------------------------------------
    ! Precondition the update of the density or not. 
    integer :: den_precon=0, den_mix=0
    !---------------------------------------------------------------------------
    ! Preconditioning matrices for the D_I_I
    real*8, allocatable :: preconX_den(:,:,:,:)
    real*8, allocatable :: preconY_den(:,:,:,:)
    real*8, allocatable :: preconZ_den(:,:,:,:) 
    !---------------------------------------------------------------------------
    ! Temporary
    real*8 :: Cden, epsden, denmom
contains

subroutine readdensit

    namelist /densit/ denmix, memory, den_precon, Cden, epsden, den_mix, denmom

    read(unit=*, nml = densit)

end subroutine readdensit

subroutine densit(iteration)
    !---------------------------------------------------------------------------
    ! Calculate all of the densities
    !---------------------------------------------------------------------------
    integer      :: i, it, wave
    integer, intent(in) :: iteration
    real(KIND=dp):: weight
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Allocation and initialization
$INITIALIZATION

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Save old density for next iteration and mixing.
    ! Note that this is only necessary at the moment for the ordinary rho
    ! density, it is the one that can make calculations unstable.
    if(.not. allocated(D_I_I_hist)) then
        allocate(D_I_I_hist(nx*ny*nz,2,memory)) ; D_I_I_hist = 0.0_dp
    endif   
    do i=1,memory-1
        D_I_I_hist(:,:,memory-i+1) = D_I_I_hist(:,:,memory-i)
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Saving the input density for mixing    
    D_I_I_hist(:,:,1) = D_I_I
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Zero the current density
    D_I_I   = 0.0
    D_Nm_Nm = 0.0

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculation by summing the densities
    do wave=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave.le.sum(HFBlocks(1:Blocks/2))) it = 1
        weight = occupations(wave)
        do i=1,mv
$EXPRESSION
        enddo
    enddo
   
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
!    ! Calculation by summing of the pairing densities
!    do wave=1,nwt
!        ! Isospin is neutron in the first half of blocks, proton in the rest
!        it = 2
!        if(wave.le.sum(HFBlocks(1:Blocks/2))) it = 1
!        weight = upairing(wave) * vpairing(wave)
!        if( iteration.ne.0)  then
!             weight = weight * (1-denmix)
!        endif
!        do i=1,mv
!!   PAIRINGEXPRESSION
!        enddo
!    enddo     

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! The asked for mixing+preconditioning scheme.
    call MassageDensity(iteration)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Calculation of the 'derived' densities, densities obtainable by 
    ! deriving other ones. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    do it=1,2
$DERIVATION  
    enddo  
end subroutine densit

subroutine MassageDensity(iteration)
    !---------------------------------------------------------------------------
    ! Operate on the density before feeding it into the rest of the program.
    !
    ! Parameters controlling the behaviour
    !
    !  memory : number of previous iterations to keep in memory
    !   
    !   
    !---------------------------------------------------------------------------
    
    integer, intent(in)   :: iteration
    real(KIND=dp),pointer :: res(:,:,:), Pres(:,:,:)
    real(KIND=dp)         :: C, eps, particles(2), mixparam
    real(KIND=dp), target :: resid(nx*ny*nz,2), Presid(nx*ny*nz,2)
    integer               :: it, sx, sy, sz, i,j,k, succes, N, iter
    
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable, save :: denupdates(:,:)
    
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable,save :: DIIS_matrix(:,:), RHS(:),residuals(:,:,:)
    integer, allocatable ,save :: PivotInfo(:)
    
    real(KIND=dp), allocatable :: TMP(:,:)
    real(KIND=dp)              :: Work(100)
    
    if(.not.allocated(preconX_den)) then
        allocate(preconX_den(nx,nx,2,2))        ; preconX_den= 0.0_dp
        allocate(preconY_den(ny,ny,2,2))        ; preconY_den= 0.0_dp
        allocate(preconZ_den(nz,nz,2,2))        ; preconZ_den= 0.0_dp
    endif
    
    if(iteration.eq.0) return
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Compute the residual
    resid = D_I_I - D_I_I_hist(:,:,1)
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Perform preconditioning if asked for.
    if(den_precon.eq.1) then 
        !-----------------------------------------------------------------------
        ! Decide on the preconditioning constants eps and c
        C   = Cden 
        eps = epsden ! Guess?

        particles(1) = neutrons
        particles(2) = protons
        !-----------------------------------------------------------------------
        ! Apply the preconditioner
        do it=1,2
            !-------------------------------------------------------------------
            ! Invert the derivatives.
            call InvertDerivatives(eps,C,preconX_den(:,:,:,it), &
            &                            preconY_den(:,:,:,it), &
            &                            preconZ_den(:,:,:,it))
        
            ! Remap to 3D
            res(1:nx,1:ny,1:nz)  =>  resid(1:nx*ny*nz,it)
            Pres(1:nx,1:ny,1:nz) => Presid(1:nx*ny*nz,it)
            
            do i=1,ny*nz
                Pres(:,i,1) =                                                  &
                &                       matmul(preconX_den(:,:,2,it),res(:,i,1))
            enddo
            do k=1,nz
                do i=1,nx
                    Pres(i,:,k) = Pres(i,:,k) +                                &
                    &                   matmul(preconY_den(:,:,2,it),res(i,:,k))
                enddo
            enddo
            do i=1,nx*ny
                Pres(i,1,:) = Pres(i,1,:) +                                    &
                &                       matmul(preconZ_den(:,:,2,it),res(i,1,:))
            enddo
        enddo
        resid = Presid 
    endif
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Perform mixing
    select case(den_mix) 
    case(0) 
        ! Simple linear mixing at the moment.
        D_I_I = D_I_I_hist(:,:,1) + (1-denmix) * resid

    case(1)
        !-----------------------------------------------------------------------
        ! DIIS mixing with (N) iterations
        if(.not.allocated(DIIS_matrix)) then
            allocate(DIIS_matrix(memory+1, memory+1)); DIIS_matrix=0
            allocate(RHS(memory+1))                  ; RHS=0
            allocate(PivotInfo(memory+1))            ; PivotInfo =0
            allocate(residuals(mv, 2, memory))       ; residuals = 0
        endif
        !-----------------------------------------------------------------------
        ! Save the output density as history.
        ! This is crucial to maintain enough information in the density 
        ! mixing algorithm.
        D_I_I_hist(:,:,1) = D_I_I
    
        !N is the number of mixed iterations
        N = min(Iteration, memory)
        
        ! DIIS matrix elements are constructed from past residuals
        do i=1,N-1
            residuals(:,:,N-i+1) = residuals(:,:,N-i) 
        enddo
        residuals(:,:,1)= resid 
        
        do i=1,N
            do j=i,N
                DIIS_matrix(i,j) = sum(residuals(:,:,i) * residuals(:,:,j))*dv
                DIIS_matrix(j,i) = DIIS_matrix(i,j)
            enddo
        enddo

        if(Iteration.le.1) then
          !---------------------------------------------------------------------
          !Return linear damping when not enough info is present
          D_I_I = D_I_I_hist(:,:,1) + (1-denmix) * resid
          return
        endif
        
        !-----------------------------------------------------------------------
        ! Small gradient solver to solve the optimization problem.
        !
        ! Needs to be cleaned up.
        !
        ! Initial guess 
        RHS(1:N)  = 1/sqrt(1.0*N)
        do iter=1,100
            RHS(1:N) = RHS(1:N)-200*RHS(1:N)*matmul(DIIS_matrix(1:N,1:N),RHS(1:N)**2)
            ! Normalize
            RHS(1:N) = RHS(1:N)/sqrt(sum(RHS(1:N)**2))
        enddo
        print *, 'RHS', RHS(1:N)
!-------------------------------------------------------------------------------
!        Quarantined code to do old-style DIIS.
!------------------------------------------------------------------------------- 
!
!        RHS                  =   0.0_dp
!        DIIS_Matrix(N+1,1:N) =  -1.0_dp
!        DIIS_Matrix(1:N,N+1) =  -1.0_dp
!        DIIS_Matrix(N+1,N+1) =   0.0_dp
!        RHS(N+1)             =  -1
!        
!        call DSYSV('L',N+1,1,DIIS_Matrix(1:N+1,1:N+1),N+1,pivotinfo,RHS(1:N+1),&
!        &           N+1,work,size(Work),Succes)
!        if(Succes.ne.0) then
!            print *, 'Error while solving the DIIS linear system.', Succes
!            stop
!        endif
!------------------------------------------------------------------------------- 

        !-----------------------------------------------------------------------
        ! Construct the new density.
        ! Note that we add an extra point into the last direction to keep
        ! on adding new information.
        D_I_I = 0.0
        do i=1,N
            D_I_I = D_I_I + RHS(i)**2*D_I_I_hist(:,:,i)
        enddo

    case(2)
         ! Density mixing with some momentum added in :)
         if(.not.allocated(DenUpdates)) then
            allocate(DenUpdates(mv,2)) ; DenUpdates = 0
         endif
         DenUpdates = (1-denmix) * resid + denmom*DenUpdates
         D_I_I = D_I_I_hist(:,:,1) + DenUpdates
!        !-----------------------------------------------------------------------
!        ! Broyden mixing algorithm
!        if(.not.allocated(residuals)) then
!            allocate(residuals(mv, 2, memory))       ; residuals = 0
!            allocate(B_matrix(memory, memory))       ; B_matrix = 0
!            allocate(beta(memory, memory))           ; beta = 0
!        endif
!        
!        !-----------------------------------------------------------------------
!        ! N is the number of mixed iterations
!        N = min(Iteration, memory)
!        do i=1,N-1
!            residuals(:,:,N-i+1) = residuals(:,:,N-i) 
!        enddo
!        residuals(:,:,1)= resid 
!        
!        !---------------------------------------------------------------------
!        D_I_I = D_I_I_hist(:,:,1) + (1-denmix) * residuals(:,:,1)
!        
!        !Return linear damping when not enough info is present
!        if(Iteration.le.1) return
!        !-----------------------------------------------------------------------
!        ! Build the matrix of residual products (not yet normalized)
!        do i=2,N
!            do j=i,N
!                B_matrix(i,j) = sum(                                           &
!                &               (residuals(:,:,i) - residuals(:,:,i-1))        &
!                &               (residuals(:,:,j) - residuals(:,:,j-1)))*dv
!                
!                B_matrix(i,j) = B_matrix(j,i)
!            enddo
!        enddo
!        ! Build the 
        
        
    end select
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Safeguard
    where(D_I_I.lt.1d-10) D_I_I = 0 
end subroutine MassageDensity
    
end module densities
