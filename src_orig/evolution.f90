module evolution
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
! Module that governs the evolution of the single-particle wavefunctions from 
! one iteration to the next. 
!
! Currently possible:
!   a) Gradient descent aka imaginary time step.
!
!===============================================================================

    use wavefunctions
    use functional
    use preconditioning

    implicit none
    
    !---------------------------------------------------------------------------
    ! Parameters of the iteration scheme
    real(KIND=dp):: dt    =  0.01
    real(KIND=dp):: hbar  =  6.58211928_dp
    !---------------------------------------------------------------------------
    ! Norm of the gradient
    real(KIND=dp) :: gradientnorm
    !---------------------------------------------------------------------------
    !Maximum number of iterations and number of iterations to skip printing of
    ! the code in the evolve subroutine
    integer :: MaxIter=100, PrintIter=10
    !---------------------------------------------------------------------------
    ! Precondition, whether to use the PG preconditioner
    character(len=20) :: Precondition = 'None'
    !---------------------------------------------------------------------------
    ! Strategy for evolution of the spwfs
    character(len=20) :: Strategy = 'IMTIME'
    !---------------------------------------------------------------------------
    !Procedure that determines the evolution of a Spwf under imaginary time.
    abstract interface
      subroutine Evolve_interface(Iteration)
        integer, intent(in)       :: iteration
      end subroutine
    end interface
    procedure(Evolve_Interface),pointer :: Evolve    
    !---------------------------------------------------------------------------
    ! Procedure pointer for the preconditioning
    procedure(Precondition_PG),pointer :: Precon 

    !---------------------------------------------------------------------------
    ! Default value of the momentum factor.
    real(KIND=dp) :: momfactor=0.1
    !---------------------------------------------------------------------------
    ! Index of the spwf which is maximised instead of minimised
    integer :: maxind = 0
    
    !---------------------------------------------------------------------------
    ! Inverse of the second order derivative matrices with appropriate constants
    real*8, allocatable :: preconX(:,:,:,:)
    real*8, allocatable :: preconY(:,:,:,:)
    real*8, allocatable :: preconZ(:,:,:,:) 
contains
    
    subroutine ReadEvolution
    !---------------------------------------------------------------------------
    ! Read the information on the evolution of the spwfs. 
    !
    !
    !---------------------------------------------------------------------------
        use geninfo

        namelist /evolution/ dt, maxiter, printiter, precondition, strategy,   &
        &                    momfactor

        read(unit=*, nml=evolution)
        !-----------------------------------------------------------------------
        !  Assign the correct preconditioner
        Precondition = to_upper(Precondition )
        if(adjustl(Precondition) .eq. 'PG' ) then
            Precon => Precondition_PG
        else 
            Precon => Precondition_none
        endif
        !-----------------------------------------------------------------------
        ! Assign the correct evolution routine
        Strategy = to_upper(Strategy)
        if(adjustl(Strategy) .eq. 'IMTIME' ) then
            Evolve => Evolve_graddesc
        elseif(adjustl(Strategy) .eq. 'MOMENTUM') then
            Evolve => Evolve_momentum
        else
            stop ('STRATEGY NOT RECOGNIZED.')
        endif

    end subroutine ReadEvolution

    subroutine PrintEvolution
        !-----------------------------------------------------------------------
        ! Print the information on the evolution strategy.
        !
        !-----------------------------------------------------------------------

        1 format(90('-'))
        2 format(' Evolution strategy: ', a20 )
        3 format('    Parameters:  dt= ', f7.4)        
        4 format(' Preconditioning   : ', a20 )
    
        print 1
        print 2, adjustl(Strategy)
        print 3, dt
        print 4, adjustl(Precondition)
        print 1
    end subroutine PrintEvolution

    subroutine Evolve_graddesc(iteration)
        !-----------------------------------------------------------------------
        ! 
        ! a) For every wave-function do a gradient step
        !    
        !    psi => ( 1 - dt/hbar h ) psi
        ! b) Calculate values:
        !    < psi | h   | psi >
        !    < psi | h^2 | psi >
        ! b) Orthonormalize within symmetry blocks
        !-----------------------------------------------------------------------
        
        use wavefunctions
        
        integer, intent(in) :: iteration
        integer             :: wave, iso,k,i
        real(KIND = dp)     :: hpsi(nx*ny*nz,4)
        
        gradientnorm = 0.0_dp

        ! Calculate the preconditioning matrices
        if(Precondition .ne. 'NONE' ) call CalculatePreconditioners()
        
        do wave=1,nwt
            if(wave .lt. nwn) then
                iso = -1
            else
                iso = +1
            endif

            hpsi = sphamil( hfpsi(:,:,wave)     ,                              &
            &              hfdpsi(:,:,:,wave)   ,                              &
            &              hfddpsi(:,:,:,wave),                                &
            &              hfdddpsi(:,:,:,wave),                               &
            &              sx(:,wave), sy(:,wave), sz(:,wave),iso)

            spenergies(wave)  = sum(hfpsi(:,:,wave) * hpsi(:,:)) * dv
            dispersions(wave) = sum( hpsi(:,:)**2)  * dv   -spenergies(wave)**2          

            gradientnorm = gradientnorm + occupations(wave) * &
            & sum((spenergies(wave) * hfpsi(:,:,wave) - hpsi(:,:))**2)*dv

            hpsi =   hpsi - spenergies(wave) * hfpsi(:,:,wave)
            hpsi =   Precon(hpsi, sx(:,wave), sy(:,wave), sz(:,wave), iso)    

            hfpsi(:,:,wave) = hfpsi(:,:,wave) -  dt/hbar * hpsi
        enddo
    
        gradientnorm = sqrt(gradientnorm)/(neutrons + protons)  
        !if(abs(gradientnorm) .lt.5d-5) stop
        call GramSchmidt
    
    end subroutine Evolve_graddesc

    subroutine Evolve_momentum(iteration)
        !-----------------------------------------------------------------------
        ! 
        ! a) For every wave-function do a gradient step, but with and added 
        !    momentum term 
        !    psi^(i+1) => ( 1 - dt/hbar h ) psi^(i) + gamma * deltapsi
        !    
        !    where deltapsi is the difference
        !       psi^(i) - psi^(i-1)
        !
        !    gamma is currently a simple constant = 0.9
        ! 
        ! b) Calculate values:
        !    < psi | h   | psi >
        !    < psi | h^2 | psi >
        ! b) Orthonormalize within symmetry blocks
        !-----------------------------------------------------------------------
        
        use wavefunctions
        
        integer, intent(in)   :: iteration
        integer               :: wave, iso,k,i
        real(KIND = dp)       :: hpsi(nx*ny*nz,4), olde
        real(KIND = dp), allocatable, save :: Updates(:,:,:)      

        if(.not.allocated(Updates)) then
            allocate(Updates(nx*ny*nz,4,nwt))
            Updates = 0.0_dp
        endif

        gradientnorm = 0.0_dp
        if(Precondition .ne. 'NONE' ) call CalculatePreconditioners()

        do wave=1,nwt
            if(wave .lt. nwn) then
                iso = -1
            else
                iso = +1
            endif
            !-------------------------------------------------------------------
            ! Calculate the single-particle hamiltonian.
            hpsi = sphamil( hfpsi(:,:,wave)     ,                              &
            &              hfdpsi(:,:,:,wave)   ,                              &
            &              hfddpsi(:,:,:,wave)  ,                              &
            &              hfdddpsi(:,:,:,wave) ,                              &
            &              sx(:,wave), sy(:,wave), sz(:,wave),iso)

            spenergies(wave)  = sum(hfpsi(:,:,wave) * hpsi(:,:)) * dv
            dispersions(wave) = sum(hpsi(:,:)**2)*dv - spenergies(wave)**2          
            
            gradientnorm = gradientnorm + occupations(wave) * &
            & sum((spenergies(wave) * hfpsi(:,:,wave) - hpsi(:,:))**2)*dv
            !-------------------------------------------------------------------
            ! Take of the part that is propagation in its own direction.
            hpsi =   hpsi - spenergies(wave) * hfpsi(:,:,wave)
            !-------------------------------------------------------------------
            ! Precondition the update if necessary.
            hpsi =   Precon(hpsi, sx(:,wave), sy(:,wave), sz(:,wave), iso)  
            !-------------------------------------------------------------------
            ! Add some history and 'momentum' to the update. 
            updates(:,:,wave) = momfactor*updates(:,:,wave) - dt/hbar * hpsi
            
            !-------------------------------------------------------------------
            ! Update the wavefunctions.
            hfpsi(:,:,wave) = hfpsi(:,:,wave) + updates(:,:,wave)
        enddo
    
        gradientnorm = sqrt(gradientnorm)/(neutrons + protons)  
        ! Orthonormalize
        call GramSchmidt
    
    end subroutine Evolve_momentum

    function Precondition_PG(psi, px, py, pz, iso) result(Ppsi)
        !-----------------------------------------------------------------------
        ! Apply a suitable preconditioner to the spwf.
        !-----------------------------------------------------------------------

        use functional

        real(KIND=dp), intent(in), target  :: psi(nx*ny*nz,4)
        real(KIND=dp), target              :: Ppsi(nx*ny*nz,4)
    
        integer, intent(in)   :: px(4),py(4),pz(4), iso
        integer               :: i,j,k, l, sx, sy, sz,  it
        real(KIND=dp),pointer :: p3(:,:,:,:), Pp3(:,:,:,:)
        
        it = (iso+3)/2
        
        p3(1:nx,1:ny,1:nz,1:4) => psi
        Pp3(1:nx,1:ny,1:nz,1:4) => Ppsi

        do l=1,4
            sx = (px(l) + 3)/2 ! These are equal to 
            sy = (py(l) + 3)/2 !    1    if pi =   -1  or 0
            sz = (pz(l) + 3)/2 !    2    if pi =   +1 
            
            do i=1,ny*nz
                Pp3(:,i,1,l) =                                                 &
                &                       matmul(preconX(:,:,sx,it),p3(:,i,1,l))
            enddo   
            do k=1,nz
                do i=1,nx
                    Pp3(i,:,k,l) = Pp3(i,:,k,l) +                              &
                    &                   matmul(preconY(:,:,sy,it),p3(i,:,k,l))
                enddo
            enddo
            do i=1,nx*ny
                Pp3(i,1,:,l) = Pp3(i,1,:,l) +                                  &
                &                       matmul(preconZ(:,:,sz,it),p3(i,1,:,l))
            enddo
        enddo
    end function Precondition_PG

    function Precondition_None(psi, px, py, pz, iso) result(Ppsi)
        !-----------------------------------------------
        ! Apply a suitable preconditioner to the spwf.
        !----------------------------------------------

        real(KIND=dp), intent(in), target :: psi(nx*ny*nz,4)
        real(KIND=dp)                     :: Ppsi(nx*ny*nz,4)
        integer, intent(in)               :: px(4),py(4),pz(4), iso
        
        Ppsi = psi
    end function Precondition_None
    
    subroutine CalculatePreconditioners
        !-----------------------------------------------------------------------
        ! Find suitable constants for use in the preconditioners and employ
        ! to calculate the preconditioning matrices.
        !-----------------------------------------------------------------------
    
        integer       ::  pm,i, loca,k, it, startind, endind
        real(KIND=dp) :: epsilon0, inproduct
        
        if(.not.allocated(preconx)) then
            allocate(preconx(nx,nx,2,2))
            allocate(precony(ny,ny,2,2))
            allocate(preconz(nz,nz,2,2))
        endif
    
        do it=1,2
            !-------------------------------------------------------------------
            ! Find a proper value for epsilon0
            epsilon0 = 0.0_dp
            
            if (it .eq. 1) then
                startind = 1
                endind   = nwn
            else
                startind = nwn+1
                endind   = nwt
            endif
            !-------------------------------------------------------------------
            ! Find the minimum sp. energy for this nucleon species.
            do i=startind, endind
                if(spenergies(i) .lt.  epsilon0) then
                    epsilon0 = spenergies(i)
                    loca = i
                endif
            enddo
            !-------------------------------------------------------------------
            ! Calculate the kinetic energy of this particular level.
            Inproduct = 0.0_dp
            do k=1,4          
                    do i=1,mv
                           Inproduct = Inproduct + HFPsi(i,k,loca) *  & 
                           &  ( HFddPsi(i,1,k,loca) + &
                           &    HFddPsi(i,4,k,loca) + &
                           &    HFddPsi(i,6,k,loca))
                    enddo
            enddo
            ! Epsilon is the potential energy, i.e. E_spwf - E_kin
            epsilon0 =   epsilon0 + hbm(it) * Inproduct * dv
            !-------------------------------------------------------------------
            ! Precalculate the inverse of the matrices
            !
            !  ( epsilon - hbar/2m * Delta)^{-1}
            ! 
            call InvertDerivatives(epsilon0, -hbm(it),preconX(:,:,:,it),       &
            &                                         preconY(:,:,:,it),       &
            &                                         preconZ(:,:,:,it))
                                           
        enddo
        
    end subroutine CalculatePreconditioners
    
!    subroutine InvertDerivatives
!    !---------------------------------------------------------------------------
!    ! Construct the inverse matrices of the second Lagrange derivative 
!    ! matrices.
!    !---------------------------------------------------------------------------
!    
!    integer :: pivotx(nx)
!    integer :: pivoty(ny)
!    integer :: pivotz(nz)    
!    integer :: ierror, pm,i, loca,k, it, startind, endind

!    real(KIND=dp), allocatable:: toinvert(:,:)
!    real(KIND=dp) :: work(nz)
!    real(KIND=dp) :: epsilon0, inproduct, epsilon0old(2) = 0.0_dp
!    
!    if(.not.allocated(invlaplaX)) then
!        allocate(invlaplaX(nx,nx,2,2))
!        allocate(invlaplaY(ny,ny,2,2))
!        allocate(invlaplaZ(nz,nz,2,2))
!    endif

!    do it=1,2
!        !-----------------------------------------------------------------------
!        ! Find a proper value for epsilon0
!        epsilon0 = 0.0_dp
!        
!        if (it .eq. 1) then
!            startind = 1
!            endind   = nwn
!        else
!            startind = nwn+1
!            endind   = nwt
!        endif

!        do i=startind, endind
!            if(spenergies(i) .lt.  epsilon0) then
!                epsilon0 = spenergies(i)
!                loca = i
!            endif
!        enddo
!        Inproduct = 0.0_dp
!        do k=1,4          
!                do i=1,mv
!                       Inproduct = Inproduct + HFPsi(i,k,loca) *  & 
!                       &  ( HFddPsi(i,1,k,loca) + &
!                       &    HFddPsi(i,4,k,loca) + &
!                       &    HFddPsi(i,6,k,loca))
!                enddo
!        enddo
!        epsilon0 =   epsilon0 - hbm(it) * Inproduct * dv
!        epsilon0 = -  epsilon0 /hbm(it)


!        !-----------------------------------------------------------------------
!        ! Invert the shifted Laplacians
!        do pm=1,2
!            invLaplaX(:,:,pm,it) =  laplaX(:,:,pm)
!            do i=1,nx
!                invLaplaX(i,i,pm,it) = invLaplaX(i,i,pm,it) - epsilon0
!            enddo
!            call dgetrf (nx, nx, invLaplaX(:,:,pm,it), nx,pivotx, ierror) 		
!            call dgetri (nx, invLaplaX(:,:,pm,it), nx, pivotx, work, nx, ierror) 
!        enddo

!        do pm=1,2
!            invLaplaY(:,:,pm,it) = laplaY(:,:,pm)
!            do i=1,ny
!                invLaplaY(i,i,pm,it) =  invLaplaY(i,i,pm,it) - epsilon0
!            enddo
!            call dgetrf (ny, ny, invLaplaY(:,:,pm,it), ny,pivoty, ierror) 		
!            call dgetri (ny, invLaplaY(:,:,pm,it), nx, pivoty, work, ny, ierror) 
!        enddo

!        do pm=1,2
!            invLaplaZ(:,:,pm,it) = laplaZ(:,:,pm)
!            do i=1,nz
!                invLaplaZ(i,i,pm,it) =  invLaplaZ(i,i,pm,it) - epsilon0
!            enddo
!            call dgetrf (nz, nz, invLaplaZ(:,:,pm,it), nz,pivotz, ierror) 		
!            call dgetri (nz, invLaplaZ(:,:,pm,it), nz, pivotz, work, nz, ierror) 
!        enddo

!        if(ierror.ne.0) then
!            print *, 'Error in inverting the shifted laplacians.'
!        endif
!    enddo
! end subroutine InvertDerivatives
end module evolution
