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
! Currently possible are 
!
! IMTIME => Gradient Descent/Imaginary Time
! HEAVYB => Heavy-ball dynamics
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
    ! Norm of the gradient and weighted sum of the dispersionss
    real(KIND=dp) :: gradientnorm, d2h
    !---------------------------------------------------------------------------
    !Maximum number of iterations and number of iterations to skip printing of
    ! the code in the evolve subroutine
    integer :: MaxIter=100, PrintIter=10
!    !---------------------------------------------------------------------------
!    ! Precondition, whether to use the PG preconditioner
!    character(len=20) :: Precondition = 'None'
    !---------------------------------------------------------------------------
    ! Strategy for evolution of the spwfs
    ! Valid choices: 
    !   IMTIME    => Gradient Descent/Imaginary Time
    !   HEAVYBALL => Heavy-ball dynamics
    character(len=20) :: Strategy = 'HEAVYBALL'
    !---------------------------------------------------------------------------
    ! Allow Tantalus to estimate the runtime parameters of the algorithm 
    ! or stay faithful to those specified by the user.
    logical :: EstimateParams = .true.
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
    !procedure(Precondition_PG),pointer :: Precon 
    !---------------------------------------------------------------------------
    ! Default value of the momentum factor.
    real(KIND=dp) :: momentum=0.0
    !---------------------------------------------------------------------------
    ! Inverse of the second order derivative matrices with appropriate constants
    real*8, allocatable :: preconX(:,:,:,:)
    real*8, allocatable :: preconY(:,:,:,:)
    real*8, allocatable :: preconZ(:,:,:,:) 

contains
    
    subroutine ReadEvolution(file_number)
        !-----------------------------------------------------------------------
        ! Read the information on the evolution of the spwfs. 
        !
        !
        !-----------------------------------------------------------------------
        use geninfo

        integer(dp), intent(in), optional   :: file_number   

        namelist /evolution/ dt, maxiter, printiter, strategy, momentum,       &
        &                    estimateparams


        if(present(file_number)) then
          read(unit=file_number, nml=evolution)
        else
          read(unit=*, nml=evolution)
        endif
        !-----------------------------------------------------------------------
        !  Assign the correct preconditioner
!        Precondition = to_upper(Precondition )
!        if(adjustl(Precondition) .eq. 'PG' ) then
!            Precon => Precondition_PG
!        else 
!            Precon => Precondition_none
!        endif
        !-----------------------------------------------------------------------
        ! Assign the correct evolution routine
        Strategy = to_upper(Strategy)
        if(adjustl(Strategy) .eq. 'IMTIME' ) then
            Evolve => Evolve_graddesc
        elseif(adjustl(Strategy) .eq. 'HEAVYBALL') then
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

        1 format(80('-'))
        2 format(' Evolution strategy: ', a20 )
        3 format('   dt= ', f7.4, ' mu= ', f7.4 )        
        4 format('   Estimate (dt,mu)  : ', a3)
!        5 format(' Preconditioning   : ', a20 )
    
        print 1
        print 2, adjustl(Strategy)
        
        if( EstimateParams) then
          print 4, 'YES'
        else 
          print 4, ' NO'
          print 3, dt, momentum
        endif
        
!        print 5, adjustl(Precondition)
    end subroutine PrintEvolution

    subroutine Evolve_graddesc(iteration)
        !-----------------------------------------------------------------------
        ! 
        ! a) For every wave-function do a gradient step
        !    
        !    psi => ( 1 - dt/hbar h ) psi
        !
        ! b) Calculate values:
        !
        !    < psi | h   | psi >
        !    < psi | h^2 | psi >
        !
        ! c) Orthonormalize within symmetry blocks
        !
        !-----------------------------------------------------------------------
        
        use wavefunctions
        
        integer, intent(in) :: iteration
        integer             :: wave, iso, iter
        real(KIND = dp)     :: hpsi(nx*ny*nz,4)
        
        gradientnorm = 0.0_dp
        d2h          = 0.0_dp

        iter = iteration      ! To get around the unused variable warnings
                              ! from compilers. Note that the variable needs to
                              ! be declared for the procedure pointers to work.

        ! Calculate the preconditioning matrices
!        if(Precondition .ne. 'NONE' ) call CalculatePreconditioners()
        
        do wave=1,nwt
            if(wave .le. nwn) then
                iso = -1
            else
                iso = +1
            endif

            hpsi = sphamil( hfpsi(:,:,wave)     ,                              &
            &              hfdpsi(:,:,:,wave)   ,                              &
            &              hfddpsi(:,:,:,wave),                                &
            &              hfdddpsi(:,:,:,wave),                               &
            &              sx(:,wave), sy(:,wave), sz(:,wave),iso,.false.)

            spenergies(wave)  = sum(hfpsi(:,:,wave) * hpsi(:,:)) * dv
            dispersions(wave) = sum( hpsi(:,:)**2)  * dv   -spenergies(wave)**2          

            select case(pairingtype)
            case(0,1)
              d2h          = d2h + rho_can(wave)*dispersions(wave)
              gradientnorm = gradientnorm + rho_can(wave) *                    &
              & sum((spenergies(wave) * hfpsi(:,:,wave) - hpsi(:,:))**2)*dv
            case(2) 
              d2h          = d2h + rho_pairing(wave,wave)*dispersions(wave)
              gradientnorm = gradientnorm + rho_pairing(wave,wave) *           &
              & sum((spenergies(wave) * hfpsi(:,:,wave) - hpsi(:,:))**2)*dv
            end select
              
            hpsi =   hpsi - spenergies(wave) * hfpsi(:,:,wave)
            !hpsi =   Precon(hpsi, sx(:,wave), sy(:,wave), sz(:,wave), iso)    

            hfpsi(:,:,wave) = hfpsi(:,:,wave) -  dt/hbar * hpsi
        enddo
    
        gradientnorm = sqrt(gradientnorm)/(neutrons + protons)  
        call GramSchmidt
    end subroutine Evolve_graddesc

    subroutine Evolve_momentum(iteration)
        !-----------------------------------------------------------------------
        ! 
        ! a) For every wave-function do a gradient step, but with and added 
        !    momentum term 
        !
        !    psi^(i+1) => ( 1 - dt/hbar h ) psi^(i) + mu * deltapsi
        !    
        !    where deltapsi is the difference
        !       psi^(i) - psi^(i-1)
        ! 
        ! b) Calculate values:
        !    < psi | h   | psi >
        !    < psi | h^2 | psi >
        !
        ! c) Orthonormalize within symmetry blocks
        !-----------------------------------------------------------------------
        
        use wavefunctions
        
        integer, intent(in)   :: iteration
        integer               :: wave, iso
        real(KIND = dp)       :: hpsi(nx*ny*nz,4)
        ! Store the change in the spwfs from last iteration
        real(KIND = dp), allocatable, save :: Updates(:,:,:)      

        if(.not.allocated(Updates)) then
            allocate(Updates(nx*ny*nz,4,nwt))
            Updates = 0.0_dp
        endif

!        if(Precondition .ne. 'NONE') then
!          print *, 'Preconditioning not supported with heavy-ball.'
!        endif  
        if(EstimateParams) call IterativeEstimation(iteration)

        gradientnorm = 0.0_dp
        d2h = 0.0_dp
        
        do wave=1,nwt
            if(wave .le. nwn) then
                iso = -1
            else
                iso = +1
            endif
            !-------------------------------------------------------------------
            ! Calculate the action of the single-particle hamiltonian.
            hpsi = sphamil( hfpsi(:,:,wave)     ,                              &
            &              hfdpsi(:,:,:,wave)   ,                              &
            &              hfddpsi(:,:,:,wave)  ,                              &
            &              hfdddpsi(:,:,:,wave) ,                              &
            &              sx(:,wave), sy(:,wave), sz(:,wave),iso,.false.)
          
            !-------------------------------------------------------------------
            spenergies(wave)  = sum(hfpsi(:,:,wave) * hpsi(:,:)) * dv
            dispersions(wave) = sum(hpsi(:,:)**2)*dv - spenergies(wave)**2          
            
            select case(pairingtype)
            case(0,1)
              d2h          = d2h + rho_can(wave)*dispersions(wave)
              gradientnorm = gradientnorm + rho_can(wave) *                    &
              & sum((spenergies(wave) * hfpsi(:,:,wave) - hpsi(:,:))**2)*dv
            case(2) 
              d2h          = d2h + rho_pairing(wave,wave)*dispersions(wave)
              gradientnorm = gradientnorm + rho_pairing(wave,wave) *           &
              & sum((spenergies(wave) * hfpsi(:,:,wave) - hpsi(:,:))**2)*dv
            end select
            
            !-------------------------------------------------------------------
            ! Remove the part that is propagation in its own direction.
            hpsi =   hpsi - spenergies(wave) * hfpsi(:,:,wave)
            !-------------------------------------------------------------------
            ! Add some history and 'momentum' to the update. 
            updates(:,:,wave) = momentum*updates(:,:,wave) - dt/hbar * hpsi
            !-------------------------------------------------------------------
            ! Update the wavefunctions.
            hfpsi(:,:,wave) = hfpsi(:,:,wave) + updates(:,:,wave)
        enddo
    
        gradientnorm = sqrt(gradientnorm)/(neutrons + protons) 
        d2h          = d2h/(neutrons+protons)
        ! Orthonormalize
        call GramSchmidt
    
    end subroutine Evolve_momentum

    subroutine IterativeEstimation(Iteration)
      !-------------------------------------------------------------------------
      ! Estimate optimum parameters (dt,mu) of the iterative process to try and
      ! achieve optimal convergence. 
      !-------------------------------------------------------------------------
            
      use wavefunctions
      
      1 format (a20, 99f10.3)
      2 format ('-----------------------------------------------------------')
      3 format (' Warning: maximum value on the mesh could not be estimated.')
      4 format (' maxE = ', f10.3,  ' convergence =', es10.3)
     
      integer, intent(in)              :: iteration
      
      real(KIND=dp), allocatable, save :: maxspwf(:,:)
      real(KIND=dp), allocatable, save :: update(:,:), actionofh(:,:)
      real(KIND=dp), allocatable, save ::   dmax(:,:,:)
      real(KIND=dp), allocatable, save ::  ddmax(:,:,:)
      real(KIND=dp), allocatable, save :: dddmax(:,:,:)
      
      integer       :: estiter, iter, sxm(4), sym(4), szm(4), ii, i
      real(KIND=dp) :: con, maxE, compare, relE, kappa
      !-------------------------------------------------------------------------
      ! Step 1: Solve the auxiliary problem for the largest single-particle 
      !         ennergy on the mesh
      if(Iteration .eq.1) then
          ! Initialize with a random spwf at the start.
          allocate(maxspwf(nx*ny*nz,4)) 
          allocate(update(nx*ny*nz,4)) ; allocate(actionofh(nx*ny*nz,4))
          allocate(dmax(nx*ny*nz,3,4))
          allocate(ddmax(nx*ny*nz,6,4))
          allocate(dddmax(nx*ny*nz,10,4))
            
          call random_number(maxspwf)                        ! randomize
          maxspwf = 1.0/sqrt(sum(maxspwf**2)*dv) * maxspwf   ! normalize
      endif

      estiter = 500
      update  = 0.0

      ! For now, assume positive parity, +i signature neutron.
      sxm(1) =  1 ; sym(1) = +1 ; szm(1) = +1
      sxm(2) = -1 ; sym(2) = -1 ; szm(2) = +1 
      sxm(3) = -1 ; sym(3) = +1 ; szm(3) = -1
      sxm(4) =  1 ; sym(4) = -1 ; szm(4) = -1

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Iterative estimation of the maximal energy
      do iter=1,estiter
          !---------------------------------------------------------------------
          ! Note that onthefly = .true., making sphamil take care of the 
          ! derivatives itself
          actionofh = sphamil(maxspwf, dmax, ddmax, dddmax,sxm,sym,szm,1,.true.)
          con       = maxE
          maxE      = sum(actionofh * maxspwf) * dv
          con       = con - maxE
          !---------------------------------------------------------------------
          ! Simple power iteration seems to better than gradient descent          
          maxspwf  = actionofh 
          !---------------------------------------------------------------------
          ! Normalize
          maxspwf = 1.0/sqrt(sum(maxspwf**2)*dv) * maxspwf   ! normalize
          !---------------------------------------------------------------------
          ! Don't be to picky about convergence, within the order of an MeV is
          ! good enough.
          if(abs(con).lt. 1d-2) exit
      enddo

      if(iter.eq.estiter+1) then
       print 1
       print 2
       print 3
       print 4, maxE, con
      endif
      !-------------------------------------------------------------------------
      ! Step 2: estimate the minimal relevant energy
      relE = 100000000
      select case(PairingType)
      case(0)   
          ! HF calculation: look at the difference in spwf energies
          do i=1,nwt
            if(abs(rho_can(i)).lt.0.5) cycle
            do ii=1,nwt
                if(abs(rho_can(ii)).gt.0.5) cycle
                compare = spenergies(ii)  - spenergies(i) 
                if(compare .gt. 0.0) then
                  relE = min(relE, compare)
                endif
            enddo
          enddo
      case(1,2)
          ! BCS or HFB calculation: look at quasiparticle energies
          relE = minval(QPenergies)
      end select
      ! Safeguard the difference
      if(relE .lt. 0.05) relE = 0.05
      !-------------------------------------------------------------------------
      ! Step 3: use these estimations to determine a value for dt and mu.
      maxE  = maxE - minval(spenergies)
      kappa = relE/maxE
      
      if(Iteration.eq.1) then
        ! Don't mess up with a too large step at the start of the iterations
        dt = 2.0/maxE*hbar* 0.90
        momentum = 0.0
      else        
        momentum = ((sqrt(kappa) - 1)/(sqrt(kappa)+1))**2
        dt    = 4.0/(maxE+relE+2*sqrt(maxE*relE))*hbar*0.90
      endif  
  end subroutine IterativeEstimation
!===============================================================================
! Projection on the feasible subspace routine
!===============================================================================  
  subroutine FeasibleProject()
  !-----------------------------------------------------------------------------
  ! Subroutine performing one (or more) alternate step for the alternating
  ! constraints. The idea is a a simple gradient step in the direction of a
  ! satisfied constraint, meaning that the objective function being minimized is
  !
  !  ( < O > - O_target )^2
  ! 
  ! and we update the single-particle wavefunctions according to 
  !
  ! psi = ( 1 - epsilon \hat{O} ) psi
  !  
  ! with
  !
  ! epsilon = 1/2 * ( <C> - C )/( < C^2 >)
  !
  ! where < C >^2 is the one-body part of the two-body operator C.
  !-----------------------------------------------------------------------------

   use wavefunctions
   use moments

   type(Moment),pointer  :: Current
   real(KIND=dp)         :: multipole(nx*ny*nz,2), update(nx*ny*nz,2)
   real(KIND=dp)         :: O2, value, des
   integer               :: it, wave, k

   Current    => Root
   multipole = 0.0_dp
   call compcutoff()
   
   do while(associated(Current%Next))
    Current => Current%next
   
    if(Current%ConstraintType.lt.2) cycle
    O2    = sum(Current%Squared)                    ! < C^2 >
    Value = sum(Current%Value)                      ! Current value of <C>
    Des   = Current%Constraint                      ! Desired final value
    
    update = 0.0
    !-----------------------------------------------------------------------
    !Calculate the update
    do it=1,2
        Update(:,it) = 0.5*(Value-Des)/O2*Cutoff(:,it)*Current%SpherHarm
    enddo
    multipole = multipole + Update
   enddo
   !---------------------------------------------------------------------------
   ! With the update in hand, we update the spwfs
   do wave=1,nwt

      it = 1
      if(wave .gt. nwn) it = 2

      !Substituting the correction
      do k=1,4
        HFPsi(:,k,wave) = (1 - multipole(:,it))*HFPsi(:,k,wave)
      enddo
    enddo
   !---------------------------------------------------------------------------
   ! Finally, orthonormalisation
   call Gramschmidt
  end subroutine feasibleproject
  
!===============================================================================
! Preconditioning routines
!===============================================================================

!    function Precondition_PG(psi, px, py, pz, iso) result(Ppsi)
!        !-----------------------------------------------------------------------
!        ! Apply a suitable preconditioner to the spwf.
!        !-----------------------------------------------------------------------

!        use functional

!        real(KIND=dp), intent(in), target  :: psi(nx*ny*nz,4)
!        real(KIND=dp), target              :: Ppsi(nx*ny*nz,4)
!    
!        integer, intent(in)   :: px(4),py(4),pz(4), iso
!        integer               :: i,j,k, l, sx, sy, sz,  it
!        real(KIND=dp),pointer :: p3(:,:,:,:), Pp3(:,:,:,:)
!        
!        it = (iso+3)/2
!        
!        p3(1:nx,1:ny,1:nz,1:4) => psi
!        Pp3(1:nx,1:ny,1:nz,1:4) => Ppsi

!        do l=1,4
!            sx = (px(l) + 3)/2 ! These are equal to 
!            sy = (py(l) + 3)/2 !    1    if pi =   -1  or 0
!            sz = (pz(l) + 3)/2 !    2    if pi =   +1 
!            
!            do i=1,ny*nz
!                Pp3(:,i,1,l) =                                                 &
!                &                       matmul(preconX(:,:,sx,it),p3(:,i,1,l))
!            enddo   
!            do k=1,nz
!                do i=1,nx
!                    Pp3(i,:,k,l) = Pp3(i,:,k,l) +                              &
!                    &                   matmul(preconY(:,:,sy,it),p3(i,:,k,l))
!                enddo
!            enddo
!            do i=1,nx*ny
!                Pp3(i,1,:,l) = Pp3(i,1,:,l) +                                  &
!                &                       matmul(preconZ(:,:,sz,it),p3(i,1,:,l))
!            enddo
!        enddo
!    end function Precondition_PG

!    function Precondition_None(psi, px, py, pz, iso) result(Ppsi)
!        !-----------------------------------------------
!        ! Apply a suitable preconditioner to the spwf.
!        !----------------------------------------------

!        real(KIND=dp), intent(in), target :: psi(nx*ny*nz,4)
!        real(KIND=dp)                     :: Ppsi(nx*ny*nz,4)
!        integer, intent(in)               :: px(4),py(4),pz(4), iso
!        
!        Ppsi = psi
!    end function Precondition_None
    
!    subroutine CalculatePreconditioners
!        !-----------------------------------------------------------------------
!        ! Find suitable constants for use in the preconditioners and employ
!        ! to calculate the preconditioning matrices.
!        !-----------------------------------------------------------------------
!    
!        integer       :: i, loca,k, it, startind, endind
!        real(KIND=dp) :: epsilon0, inproduct
!        
!        if(.not.allocated(preconx)) then
!            allocate(preconx(nx,nx,2,2))
!            allocate(precony(ny,ny,2,2))
!            allocate(preconz(nz,nz,2,2))
!        endif
!    
!        do it=1,2
!            !-------------------------------------------------------------------
!            ! Find a proper value for epsilon0
!            epsilon0 = 0.0_dp
!            
!            if (it .eq. 1) then
!                startind = 1
!                endind   = nwn
!            else
!                startind = nwn+1
!                endind   = nwt
!            endif
!            !-------------------------------------------------------------------
!            ! Find the minimum sp. energy for this nucleon species.
!            do i=startind, endind
!                if(spenergies(i) .lt.  epsilon0) then
!                    epsilon0 = spenergies(i)
!                    loca = i
!                endif
!            enddo
!            !-------------------------------------------------------------------
!            ! Calculate the kinetic energy of this particular level.
!            Inproduct = 0.0_dp
!            do k=1,4          
!                    do i=1,mv
!                           Inproduct = Inproduct + HFPsi(i,k,loca) *  & 
!                           &  ( HFddPsi(i,1,k,loca) + &
!                           &    HFddPsi(i,4,k,loca) + &
!                           &    HFddPsi(i,6,k,loca))
!                    enddo
!            enddo
!            ! Epsilon is the potential energy, i.e. E_spwf - E_kin
!            epsilon0 =   epsilon0 + hbm(it) * Inproduct * dv
!            !-------------------------------------------------------------------
!            ! Precalculate the inverse of the matrices
!            !
!            !  ( epsilon - hbar/2m * Delta)^{-1}
!            ! 
!            call InvertDerivatives(epsilon0, -hbm(it),preconX(:,:,:,it),       &
!            &                                         preconY(:,:,:,it),       &
!            &                                         preconZ(:,:,:,it))
!                                           
!        enddo
!        
!    end subroutine CalculatePreconditioners
end module evolution
