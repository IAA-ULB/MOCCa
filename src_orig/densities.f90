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
    real(KIND=dp) :: denmix = 0.65_dp
    
    !---------------------------------------------------------------------------
    ! Previous value(s) of the density rho.
    real(KIND=dp), allocatable :: D_I_I_hist(:,:,:)
    !---------------------------------------------------------------------------
    ! The amount of iterations to keep in memory for the density mixing
    integer           :: memory = 1
    !---------------------------------------------------------------------------
    ! The type of mixing to perform.
    character(len=20) :: mixing='linear'
    
    !---------------------------------------------------------------------------
    ! Preconditioning matrices for the D_I_I
    real*8, allocatable :: preconX_den(:,:,:,:)
    real*8, allocatable :: preconY_den(:,:,:,:)
    real*8, allocatable :: preconZ_den(:,:,:,:) 
contains

subroutine readdensit

    namelist /densit/ denmix, memory, mixing

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
    ! Mixing the ordinary density via the asked for mixing+preconditioning 
    ! scheme.
    if(iteration.ne.0) call MassageDensity
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Save old density for next iteration.
    if(.not. allocated(D_I_I_hist)) then
        allocate(D_I_I_hist(nx*ny*nz,2,memory)) ; D_I_I_hist = 0.0_dp
    endif   
    D_I_I_hist(:,:,1) = D_I_I
    if(any(D_I_I.lt.0.0_dp)) print *, 'negative density'
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Calculation of the 'derived' densities, densities obtainable by 
    ! deriving other ones. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    do it=1,2
$DERIVATION  
    enddo  
end subroutine densit

subroutine MassageDensity
    !---------------------------------------------------------------------------
    ! Operate on the density before feeding it into the rest of the program.
    !
    !
    !---------------------------------------------------------------------------
    
    real(KIND=dp),pointer :: res(:,:,:), Pres(:,:,:)
    real(KIND=dp)         :: C, eps, particles(2)
    real(KIND=dp), target :: resid(nx*ny*nz), Presid(nx*ny*nz)
    integer               :: it, sx, sy, sz, i,j,k
    
    if(.not.allocated(preconX_den)) then
        allocate(preconX_den(nx,nx,2,2))        ; preconX_den= 0.0_dp
        allocate(preconY_den(ny,ny,2,2))        ; preconY_den= 0.0_dp
        allocate(preconZ_den(nz,nz,2,2))        ; preconZ_den= 0.0_dp
    endif
    
    if(mixing .ne. 'precon') return
    
    !---------------------------------------------------------------------------
    ! Decide on the preconditioning constants eps and c
    ! For now guess c to be B3 of SLy4
    C   = - 10 
    eps = - 300 ! Guess?

    particles(1) = neutrons
    particles(2) = protons
    print *, 'Mixing is happening!'
    !---------------------------------------------------------------------------
    ! Apply the preconditioner
    do it=1,2
        !-----------------------------------------------------------------------
        ! Invert the derivatives.
        call InvertDerivatives(eps,C,preconX_den(:,:,:,it), &
        &                            preconY_den(:,:,:,it), &
        &                            preconZ_den(:,:,:,it))
    
        ! Compute the residual
        resid = D_I_I(:,it) - D_I_I_hist(:,it,1)
        ! Remap to 3D
        res(1:nx,1:ny,1:nz)  =>  resid(1:nx*ny*nz)
        Pres(1:nx,1:ny,1:nz) => Presid(1:nx*ny*nz)
        
        do i=1,ny*nz
            Pres(:,i,1) =                                                      &
            &                       matmul(preconX_den(:,:,2,it),res(:,i,1))
        enddo
        do k=1,nz
            do i=1,nx
                Pres(i,:,k) = Pres(i,:,k) +                                  &
                &                   matmul(preconY_den(:,:,2,it),res(i,:,k))
            enddo
        enddo
        do i=1,nx*ny
            Pres(i,1,:) = Pres(i,1,:) +                                      &
            &                       matmul(preconZ_den(:,:,2,it),res(i,1,:))
        enddo
        

        print *, sum(res)*dv, sum(Pres)*dv, sum(D_I_I(:,it))*dv      
!        D_I_I(:,it) = D_I_I_hist(:,it,1) + Presid(:)
!        print *, 'norm', sum(Presid**2)*dv
!        ! Renormalizing
!        D_I_I(:,it) = D_I_I(:,it)* (particles(it)/(sum(D_I_I(:,it)))/dv)
    enddo
    !---------------------------------------------------------------------------
    
end subroutine MassageDensity
    
end module densities
