module evolution
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
!
! Module that governs the evolution of the single-particle wavefunctions from 
! one SCF iteration to the next. This evolution is determined principally by
! one single keyword in the Evolution namelist.
! 
! This can currently take the following values:
!     IMTIME => Gradient Descent/Imaginary Time
!     HEAVYB => Heavy-ball dynamics
!
! Note that setting these keywords has many consequences, from convergence speed
! to memory requirements, the way a calculation scales to multiple MPI ranks 
! up to the very definition of what constitutes an "iteration". Be sure to 
! spend some time thinking about what option to choose!
!
! These settings are "catch-all" and by default define many things about the 
! evolution, althought several other flags can determine specific aspects of 
! evolution. 
!
! 1. ortho_strategy      : 
! 2. subspace_rotation   :
! 3. spwf_preconditioning:
! 
! Not all combinations of options are valid inputs!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Hephaestos keywords
!
! N2               : $N2
! N3               : $N3
!===============================================================================

    use wavefunctions
    use functional
    use preconditioning
    use timing

    implicit none

    !---------------------------------------------------------------------------
    ! Parameters of the iteration scheme
    real(KIND=dp):: dt    =  0.01
    real(KIND=dp):: hbar  =  6.58211928_dp
    real(KIND=dp):: stepsize_safety = 0.9
    !---------------------------------------------------------------------------
    ! Default value of the momentum factor.
    real(KIND=dp) :: momentum=0.0
    !---------------------------------------------------------------------------
    ! Norm of the gradient and weighted sum of the dispersions
    real(KIND=dp) :: gradientnorm, d2h
!    !---------------------------------------------------------------------------
!    ! Precondition, whether to use the PG preconditioner
!    character(len=20) :: Precondition = 'None'
    !---------------------------------------------------------------------------
    ! Strategy determination for the evolution of the spwfs
    ! Strategy = the "global" selection of algorithm 
    !   IMTIME    => Gradient Descent/Imaginary Time
    !   HEAVYBALL => Heavy-ball dynamics
    !   HBSANE    => "Sane" heavy-ball dynamics
    character(len=20) :: Strategy = 'HEAVYBALL'
    ! Orthonormalisation strategy
    !   GRAMSCHMIDT => Gram-Schmidt "sequential" orthonormalisation
    !   CHOLESKY    => Cholesky decomposition
    character(len=20)               :: ortho_strategy = 'GramSchmidt'
    ! Subspace rotation 
    !    whether or not to throw in an explicit diagonalisation of the 
    !    single-particle hamiltonian in the subspace spanned by the spwfs
    !    in memory.
    logical :: subspace_rotation = .false.
    !---------------------------------------------------------------------------
    !Procedure that determines the evolution of a Spwf under imaginary time.
    abstract interface
      subroutine Evolve_interface(Iteration)
        integer, intent(in)       :: iteration
      end subroutine
    end interface
    procedure(Evolve_Interface),pointer :: Evolve    
    !---------------------------------------------------------------------------
    ! Allow Tantalus to estimate the runtime parameters of the heavy-ball 
    ! algorithm for the linear subproblem or stay faithful to those specified 
    ! by the user. 
    logical :: EstimateParams     = .true.
    !---------------------------------------------------------------------------
    ! Procedure pointer for the preconditioning
    !procedure(Precondition_PG),pointer :: Precon 
    !---------------------------------------------------------------------------
    ! Inverse of the second order derivative matrices with appropriate constants
    real*8, allocatable :: preconX(:,:,:,:)
    real*8, allocatable :: preconY(:,:,:,:)
    real*8, allocatable :: preconZ(:,:,:,:) 
contains
    
    subroutine ReadEvolution(file_number)
        !-----------------------------------------------------------------------
        ! Read the information on the evolution of the spwfs. 
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Input:
        !   file_number : optional integer. If present, read from (open) channel
        !                 with this number. If absent, read from STDIN.
        !-----------------------------------------------------------------------
        use geninfo

        integer(dp), intent(in), optional   :: file_number
#if(USE_MPI>0)
        integer                             :: mpi_err
#endif

        namelist /evolution/ dt, momentum,                                     &
        &                    gradient_stepsize, gradient_mu,                   &
        &                    maxiter, printiter, strategy,                     &
        &                    estimateparams, estimategradparams,               &
        &                    gradient_safety, efficientHFB, ortho_strategy,    &
        &                    stepsize_safety, freezeiter
        !-----------------------------------------------------------------------
        ! Only the very first MPI rank reads the input
        if(MPI_RANK.eq.0) then
          if(present(file_number)) then
            read(unit=file_number, nml=evolution)
          else
            read(unit=*, nml=evolution)
          endif
        endif
#if(USE_MPI > 0)
        !-----------------------------------------------------------------------
        ! Broadcasting from rank 0 to the rest
        call MPI_BCAST(dt       , 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(momentum , 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(gradient_stepsize, 1, MPI_REAL8, 0, &
        &                                               MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(gradient_safety  , 1, MPI_REAL8, 0, &
        &                                               MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(stepsize_safety  , 1, MPI_REAL8, 0, &
        &                                               MPI_COMM_WORLD, mpi_err)

        call MPI_BCAST(estimateparams,     1, MPI_LOGICAL, 0, &
        &                                               MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(estimategradparams, 1, MPI_LOGICAL, 0, &
        &                                               MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(efficientHFB      , 1, MPI_LOGICAL, 0, &
        &                                               MPI_COMM_WORLD, mpi_err)
        
        call MPI_BCAST(strategy  ,len(strategy), MPI_CHARACTER, 0, & 
        &                                               MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(ortho_strategy, len(ortho_strategy), MPI_CHARACTER, 0, &
        &                                               MPI_COMM_WORLD, mpi_err)

        call MPI_BCAST(maxiter   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_BCAST(printiter , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
        !-----------------------------------------------------------------------
#endif
        !-----------------------------------------------------------------------
        ! Bookkeeping that each rank should do
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! a) Assign the correct preconditioner
!        Precondition = to_upper(Precondition )
!        if(adjustl(Precondition) .eq. 'PG' ) then
!            Precon => Precondition_PG
!        else 
!            Precon => Precondition_none
!        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! b) assign the correct evolution routine
        Strategy = to_upper(Strategy)
        if(adjustl(Strategy) .eq. 'IMTIME' ) then
            Evolve => Evolve_graddesc
        elseif(adjustl(Strategy) .eq. 'HEAVYBALL') then
            Evolve => Evolve_momentum
        elseif(adjustl(Strategy) .eq. 'HBSANE') then
            Evolve => Evolve_momentum_sane
        else
            call stp('STRATEGY NOT RECOGNIZED.')
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! c) set diagsphamil 
        ! - If we use the heavy-ball algorithm for the pairing subproblem, we 
        !   limit the heavy-ball algorithm in the linear subproblem to 
        !   optimising the relevant subspace and not in diagonalising the 
        !   individual spwfs.
        if(pairingscheme .eq. 1) then
          diagsphamil = .false.
        else
          diagsphamil = .true.
        endif
        ! - if efficientHFB is true, we also do not diagonalise h
        if(efficientHFB) diagsphamil = .false.
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! d) assign the correct orthonormalisation routine
        ortho_strategy = to_upper(ortho_strategy)
        if(adjustl(ortho_strategy) .eq. 'GRAMSCHMIDT' ) then
            Orthonormalize => GramSchmidt
        elseif(adjustl(ortho_strategy) .eq. 'LOEWDIN') then
            call stp('Implementation error')
!            Orthonormalize => Loewdin
        elseif(adjustl(ortho_strategy) .eq. 'CHOLESKY') then
            Orthonormalize => Cholesky_orthonormalisation
        else
            call stp('Orthonormalisation strategy not recognized.')
        endif
        ! - - - 
    end subroutine ReadEvolution

    subroutine PrintEvolution
        !-----------------------------------------------------------------------
        ! Print the information on the way the spwfs are evolved in this 
        ! calculation.
        !-----------------------------------------------------------------------

        1 format(80('-'))
        2 format(' Evolution strategy: ', a20 )
        3 format('   dt= ', f7.4, ' mu= ', f7.4 )
       31 format('   maxiter =', i5, ' printiter = ', i5)        
       32 format('   of which freezeiter= ', i5, ' do change the potentials.')
        4 format('   Estimate (dt,mu) linear subproblem  : ', a3)
       41 format('   Safety factor for linear subproblem : ', f7.4)
       42 format('   Estimate (dt,mu) pairing subproblem : ', a3)
       43 format('   Safety HFB-gradient                 : ', f7.4)

!        5 format(' Preconditioning   : ', a20 )
        6 format(' Diagonalise the s.p. hamiltonian: ', a3)
        7 format(' EfficientHFB : ACTIVE! ')
        8 format(' Orthonormalisation strategy: ', a20)

        print 1
        print 2, adjustl(Strategy)
        print 31, maxiter, printiter
        print 32, freezeiter
        if( EstimateParams) then
          print 4, 'YES'
          print 41, stepsize_safety
        else 
          print 4, ' NO'
          print 3, dt, momentum
        endif

        if( EstimateGRADParams) then
          print 42, 'YES'
          print 43, gradient_safety
        else 
          print 42, ' NO'
          print 3, gradient_stepsize, gradient_mu
        endif
!        print 5, adjustl(Precondition)

        if(efficientHFB) print 7
        if(diagsphamil) then
            print 6, 'YES'
        else
            print 6, 'NO'
        endif
        print 8, ortho_strategy

    end subroutine PrintEvolution

!===============================================================================
! Evolution routines 
!===============================================================================

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

#if(USE_MPI > 0)
        call stp('Subroutine evolve_graddesc is not ready for use in MPI calculations.')
#endif
        do wave=1,nwt
            if(wave .le. nwn) then
                iso = -1
            else
                iso = +1
            endif

            hpsi = sphamil( hfpsi(:,:,wave)     ,                              &
            &              hfdpsi(:,:,:,wave)   ,                              &
            &              hfddpsi(:,:,:,wave),                                &
$N3         &              hfdddpsi(:,:,:,wave),                               &
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

        call orthonormalize

    end subroutine Evolve_graddesc

    subroutine Evolve_momentum(iteration)
        !-----------------------------------------------------------------------
        ! 
        ! Evolution of the single-particle wavefunctions in memory through
        ! heavy-ball evolution. 
        !
        ! (0) If asked for, estimate the evolution parameters dt and mu.
        ! (1) For every wave-function do a gradient step, but with and added 
        !     momentum term 
        !
        !    psi^(i+1) => psi^(i) - dt/hbar * g^(i)  + mu * deltapsi
        !    
        !    where 
        !      *) deltapsi is the difference  psi^(i) - psi^(i-1)
        !      *) g^(i) is the relevant gradient of the s.p. wavefunction.
        !
        !
        !    g^(i) can have two different forms:
        !      (i) If the code is attempting to diagonalize the s.p. hamiltonian
        !          completely, then the gradient is the "full" gradient:
        !
        !          g^(i) =  h^(i) | psi^(i) > - eps^i  | psi^(i) > 
        !          eps^i =  < psi^(i) | h^(i) | psi^(i) >
        !          
        !          where the subtraction just takes away the gradient in the
        !          direction of the spwf itself.
        !
        !     (ii) If the code is trying to construct the space spanned by 
        !          the lowest eigenstates of h (without resolving the 
        !          individual eigenstates by iteration), then we take 
        !     
        !           g^(i) =  h^(i) | psi^(i) > - \sum_j \lambda_j | psi_j > 
        !           lambda_j =  < psi^(j) | h^(i) | psi^(i) >
        !
        !           where the sum runs over ALL single-particle states in 
        !           memory. Hence, we are only exploring "new" directions that
        !           can actually lower the Routhian, without shuffling around
        !           single-particle states in memory.
        !
        ! (2) To end the evolution, orthonormalize all the s.p.w.fs.
        !
        ! As side-effects, the code calculates 
        !
        !   (a) The Hartree-Fock transformation linking the spwfs in memory
        !       to the Hartree-Fock basis. If we diagsphamil =.true., then this
        !       is the trivial transformation. Otherwise, this array is 
        !       constructed by diagonalising the sp. hamiltonian in the 
        !       subspace of the spwfs in memory, BEFORE the heavy-ball 
        !       evolution. Hence, it diagonalizes
        !
        !         h_ij = <psi_i | h | \psi_j >
        !
        !   (b) spenenergies array
        !       If diagonalizing the sp hamiltonian:
        !            spenergies(i )=  < psi^(i) | h^(i) | psi^(i) >     
        !       Else, the eigenvalues of the sphamiltonian obtained from the 
        !       diagonalisation mentioned above.
        !
        !   (c) gradientnorm
        !       Convergence measure that measures the size of the gradients.
        !   (d) d2h:
        !       Weighted dispersion of the spwfs, only calculated when
        !       diagsphamil = .true.
        !-----------------------------------------------------------------------

        use wavefunctions

        integer, intent(in)   :: iteration
        integer               :: wave, iso, B, si, N, wave2, lwork, ifail
        integer               :: wg, wg2
#if(USE_MPI>0)
        integer               :: mpi_err
#endif
        real(KIND=dp), allocatable :: work(:)
        real(KIND=dp)              :: hpsi(nx*ny*nz,4) 

        call start_timer(T_evolution)

        if(.not.allocated(Momentum_Updates)) then
            ! we only store the history for the LOCALLY stored wavefunctions
            allocate(Momentum_Updates(nx*ny*nz,4,nwt_local))
            Momentum_Updates = 0.0_dp
        endif

        if(.not.allocated(current_sph)) then 
            allocate(current_sph(nwt,nwt)) ; current_sph = 0.0d0
        endif

        if(EstimateParams) call IterativeEstimation(iteration)

        si           = 0
        gradientnorm = 0.0_dp
        d2h          = 0.0_dp
        hftransfo    = 0.0d0
        spenergies   = 0.0d0
        dispersions  = 0.0d0
        current_sph  = 0.0d0
        do B=1,8                       !<---- this loops over local spwf indices
          N = HFblocks(B) ; if(N.eq.0) cycle
          iso = -1
          if(B.gt.4) iso = +1
          do wave=si+1,si+N     ! = local index of the spwf
            wg = spwf_map(wave) ! = global index of the spwf
            !-------------------------------------------------------------------
            ! Calculate the action of the single-particle hamiltonian.
            hpsi = sphamil( hfpsi(:,:,wave)     ,                              &
            &              hfdpsi(:,:,:,wave)   ,                              &
            &              hfddpsi(:,:,:,wave)  ,                              &
$N3         &              hfdddpsi(:,:,:,wave) ,                              &
            &              sx(:,wave), sy(:,wave), sz(:,wave),iso,.false.)

            if(diagsphamil) then
              ! If we are diagonalising the s.p. hamiltonian, we use hpsi to
              ! calculate
              !       spenergies  : <psi|h  | psi> 
              !       dispersions : <psi|h h| psi> - spenergies**2       
              spenergies(wg)  = sum(hfpsi(:,:,wave) * hpsi(:,:)) * dv
              dispersions(wg) = sum(hpsi(:,:)**2)*dv - spenergies(wg)**2              
              ! d2h can be calculated as a convergence measure in this case
              select case(pairingtype)
              case(0,1)
                d2h          = d2h + rho_can(wg)*dispersions(wg)
              case(2) 
                d2h          = d2h + rho_pairing(wg,wg)*dispersions(wg)
              end select
            endif

#if(PASTA == 0)
            !-------------------------------------------------------------------
            ! We construct the matrix elements of the single-particle 
            ! hamiltonian in the basis of s.p. wavefunctions in memory.
#if(USE_MPI>0)
            do wave2=wave,si+N           ! local index
                wg2 = spwf_map(wave2)    ! global index

                select case (balancing_strategy)
                case(1)
                  ! In the case of symmetry-block-wise load balancing, we can
                  ! safely assume all relevant wavefunctions are represented
                  ! on the current MPI rank and we do not need more complicated
                  !  things.
                  current_sph(wg2,wg)  = sum(hfpsi(:,:,wave2) * hpsi(:,:))* dv
                  current_sph(wg ,wg2) = current_sph(wg2,wg)
                case DEFAULT
                  call stp('Subroutine evolve_momentum is not yet ready for &
                  &         MPI calculations with balancing_strategy different &
                  &         from 1.')
                end select
            enddo
#else
            do wave2=wave,si+N           ! local index
                wg2 = spwf_map(wave2)    ! global index
                current_sph(wg2,wg)  = sum(hfpsi(:,:,wave2) * hpsi(:,:))* dv
                current_sph(wg ,wg2) = current_sph(wg2,wg)
            enddo
#endif
#else 
            !-------------------------------------------------------------------
            ! ... but this is very costly when nwt is large, i.e. when doing 
            ! pasta  calculations, so we do something simple instead.
            !-------------------------------------------------------------------
            current_sph(wg,wg) = spenergies(wg)
#endif
            !-------------------------------------------------------------------
            if(diagsphamil) then
              ! We are diagonalising the s.p. hamiltonian completely, i.e.
              ! also in the space spanned by the s.p. wavefunctions in memory
              ! We simply move in the direction of h|psi\rangle
              hpsi =   hpsi - spenergies(wg) * hfpsi(:,:,wave)
            else
              ! We do not diagonalise the s.p. hamiltonian, we only try to 
              ! construct the subspace spanned by its lowest eigenstates. 
              ! Hence, our update should only take us into "new" directions.
              ! So, we orthogonalise the update direction to all spwfs in
              ! storage.
              call start_timer(T_Hortho)
              do wave2=si+1,si+N          ! local index
                wg2 = spwf_map(wave2)     ! global index
                hpsi = hpsi - current_sph(wg,wg2)*hfpsi(:,:,wave2)
              enddo
              call stop_timer(T_Hortho)
              ! The norm of the gradient can always be calculated as a 
              ! convergence measure.
              gradientnorm = gradientnorm + sum(hpsi**2)*dv
            endif

            !-------------------------------------------------------------------
            ! Add some history and 'momentum' to the update. 
            momentum_updates(:,:,wave) = &
            &               momentum*momentum_updates(:,:,wave) - dt/hbar * hpsi
          enddo

          do wave=si+1, si+N
            !-------------------------------------------------------------------
            ! Update the wavefunctions.
            hfpsi(:,:,wave) = hfpsi(:,:,wave) + momentum_updates(:,:,wave)
          enddo
          si = si + N
        enddo                     !<---- end of the loop over local spwf indices

        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Bookkeeping for convergence criteria
        gradientnorm = sqrt(gradientnorm) 
        d2h          = d2h/(neutrons+protons)

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Orthonormalize the new spwf basis.
        call orthonormalize

#if(USE_MPI > 0)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Collecting all arrays on all MPI ranks. The ALLREDUCE calls are valid, 
        ! since we zeroed the initial arrays at the top of this routine.
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! nwt-scalars
        call MPI_ALLREDUCE(MPI_IN_PLACE,d2h         , 1, MPI_REAL8, MPI_SUM, &
        &                                                MPI_COMM_WORLD,mpi_err)
        call MPI_ALLREDUCE(MPI_IN_PLACE,gradientnorm, 1, MPI_REAL8, MPI_SUM, &
        &                                                MPI_COMM_WORLD,mpi_err)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! nwt-vectors
        call MPI_ALLREDUCE(MPI_IN_PLACE,dispersions , nwt, MPI_REAL8,          &
        &                                       MPI_SUM, MPI_COMM_WORLD,mpi_err)
        if(diagsphamil) then
          call MPI_ALLREDUCE(MPI_IN_PLACE,spenergies,                          &
          &                  nwt, MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,mpi_err)
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! nwt-Matrices
        ! Note, this can be done block-wise in order to save on communication
        !       it CANNOT be included in the previous loop since it relies on 
        !       global indices
        call MPI_ALLREDUCE(MPI_IN_PLACE,current_sph , nwt**2, MPI_REAL8,       &
        &                                       MPI_SUM, MPI_COMM_WORLD,mpi_err)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#endif

        !-----------------------------------------------------------------------
        ! Treatment of the HF transformation
        ! Note: this has been separated from the previous loop since it really
        !       needs to consider the entire sphamiltonian in a given subblock.
        !       Hence a new loop over symmetry blocks.
        !-----------------------------------------------------------------------
        si = 0
        do B=1,8                      !<---- this loops over global spwf indices
          N = HFblocks_global(B) ; if(N.eq.0) cycle
          iso = -1
          if(B.gt.4) iso = +1
          if(diagsphamil) then
              ! The HF-transfo we use is trivial, i.e. the s.p. wavefunctions
              ! in memory constitute the HF basis.
              do wave=1,N
                  hftransfo(si+wave,si+wave) = 1.0d0
              enddo
          else
              ! The s.p. wavefunctions in memory do not constitute the HF-basis.
              ! We diagonalise the single-particle hamiltonian (*) to 
              ! obtain the unitary transformation into that basis.
              !
              ! (*) Actually, the single-particle hamiltonian BEFORE the 
              !     heavy-ball evolution.
              call start_timer(T_HFdiag)

              HFtransfo(si+1:si+N,si+1:si+N) = current_sph(si+1:si+N, si+1:si+N)

              lwork = -1; allocate(work(1))
              call DSYEV( 'V', 'U', N, HFtransfo(si+1:si+N,si+1:si+N), N, &
              &                       spenergies(si+1:si+N),work,lwork,ifail)
              lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
              call DSYEV( 'V', 'U', N, HFtransfo(si+1:si+N,si+1:si+N), N, &
              &                       spenergies(si+1:si+N),work,lwork,ifail)
              deallocate(work)

              call stop_timer(T_HFdiag)
          endif

          si = si + N
        enddo

        call stop_timer(T_evolution)

    end subroutine Evolve_momentum
    
    subroutine Evolve_momentum_sane(iteration)
    !
    !
    !
    !
    !
      use wavefunctions

      integer, intent(in)        :: iteration

      integer                    :: si, m, B, N, iso
      real(KIND=dp), allocatable :: hpsi(:,:,:)

      call start_timer(T_evolution)

      if(.not.allocated(Momentum_Updates)) then
          ! we only store the history for the LOCALLY stored wavefunctions
          allocate(Momentum_Updates(nx*ny*nz,4,nwt_local))
          Momentum_Updates = 0.0_dp
      endif

      if(.not.allocated(current_sph)) then 
          allocate(current_sph(nwt,nwt)) ; current_sph = 0.0d0
      endif

      current_sph = 0.0d0
      hftransfo   = 0.0d0
      d2h         = 0.0d0
      if(EstimateParams) call IterativeEstimation(iteration)

      !-------------------------------------------------------------------------
      ! Step 1: construct all updates
      si = 0
      do B=1,8
        N = HFBlocks(B) ; if(N.eq.0) cycle
        iso = -1        ; if(B.gt.4) iso = +1

        allocate(hpsi(mv,4,N))

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Obtain the action of the s.p.h. on the spwfs
        call apply_sphamil_block(N,HFpsi(:,:,si+1:si+N),hpsi,&
        &                          sx(:,si+1),sy(:,si+1),sz(:,si+1),iso)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Construct the residual
        do m=1,N
          hpsi(:,:,m)=hpsi(:,:,m) - sum(hpsi(:,:,m)*HFpsi(:,:,si+m))*dv*HFPsi(:,:,si+m)
        enddo
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Add some history and 'momentum' to the update. 
        momentum_updates(:,:,si+1:si+N) = & 
        &     - dt/hbar * hpsi      +  momentum*momentum_updates(:,:,si+1:si+N) 
        ! actual update
        deallocate(hpsi)
        si = si + N
      enddo
      !-------------------------------------------------------------------------
      ! Step 2: perform the update
      HFPSI = HFPSI + momentum_updates
      !-------------------------------------------------------------------------
      ! Step 3: orthonormalize
      call orthonormalize
      !-------------------------------------------------------------------------
      ! Step 4: diagonalize the s.p. hamiltonian block-by-block
      si = 0
      do B=1,8
        N = HFBlocks(B) ; if(N.eq.0) cycle
        iso = -1        ; if(B.gt.4) iso = +1

        allocate(hpsi(mv,4,N))
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Obtain the action of the s.p.h. on the spwfs
        call apply_sphamil_block(N,HFpsi(:,:,si+1:si+N),hpsi,&
        &                          sx(:,si+1),sy(:,si+1),sz(:,si+1),iso)
        
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Calculate the norm of the residuals
        do m=1,N
          dispersions(si+m) = sum((hpsi(:,:,m) - sum(hpsi(:,:,m)*HFpsi(:,:,si+m))*dv*HFPsi(:,:,si+m))**2)*dv
          d2h          = d2h + rho_can(si+m)*dispersions(si+m)
       enddo
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Actually diagonalize
        call diag_sph(N,N,HFPsi(:,:,si+1:si+N),hpsi, &
        &             momentum_updates(:,:,si+1:si+N), &
        &             sx(:,si+1),sy(:,si+1),sz(:,si+1),iso,&
        &             spenergies(si+1:si+N))
        deallocate(hpsi)
        si = si + N
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Populate current_sph and hftransfo
        do m=1,N
          current_sph(si+m, si+m) = spenergies(si+m)
          hftransfo(si+m,si+m)    = 1.0d0
        enddo

        si = si + N
      enddo
      d2h          = d2h/(neutrons+protons)

      call stop_timer(T_evolution)
    end subroutine Evolve_momentum_sane 
    
!===============================================================================
! Utility routines 
!===============================================================================

    subroutine apply_sphamil_block(m, x,hx,sx,sy,sz,iso)
      !-------------------------------------------------------------------------
      ! Apply the single-particle hamiltonian as specified by the potentials  
      ! currently in memory to a set of vectors in single-particle space with 
      ! common symmetry properties.
      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !
      ! Input:
      !       m :  number of vectors to compute h|psi> for
      !       x :  a set of vectors in s.p. space, i.e. a matrix of dimension
      !         (nx*ny*nz,4,m)
      ! sx/sy/sz:  the symmetries under plane reflections of the vectors, i.e.
      !            4 integers in each case.
      !      iso:  the isospin of the vectors.
      !
      ! Output: 
      !   hx :  the application of the s.p. hamiltonian to the vectors, a new
      !         matrix of dimension (nx*ny*nz,4,m)
      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Important to note: this routine does NOT reuse derivatives, i.e. it 
      ! tells the sphamiltonian to apply all derivative matrix multiplications 
      ! on the fly. 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Optional TODO: write an additional routine that takes as input the
      ! fully vectorized spwfs, i.e. psi(4*nx*ny*nz,m).
      !-------------------------------------------------------------------------
      integer, intent(in)                  ::  m, sx(4), sy(4), sz(4), iso
      real(KIND=dp), intent(in), target    ::  x(mv,4,m)
      real(KIND=dp), intent(out), target   :: hx(mv,4,m)
      integer                              :: wave
      real(KIND=dp)                        :: dummy_d (mv,3,4)
      real(KIND=dp)                        :: dummy_dd(mv,6,4)
 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do wave=1,m
        hx(:,:,wave)=sphamil(x(:,:,wave),dummy_d, dummy_dd, sx,sy,sz,iso,.true.)
      enddo
    end subroutine apply_sphamil_block

    subroutine diag_sph(m,n,x,hx,upd,sxb,syb,szb,iso,eigenvalues)
      !------------------------------------------------------------------------
      ! 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !
      ! Input:
      !       m : number of vectors passed in that span the reduced space
      !       n : number of eigenstates to construct
      !       x : a set of vectors in s.p. space, 
      !           i.e. a matrix of dimension (nx*ny*nz,4,m)
      !      hx : the application of h on the vectors x, 
      !           i.e. another matrix of dimension (nx*ny*nz,4,m)
      !      upd: a set of vectors that needs to undergo the same unitary 
      !           transformation as x
      ! sx/y/zb :  the symmetries under plane reflections of the vectors,
      !            i.e. 4 integers in each case.
      !      iso:  the isospin of the vectors.
      !
      ! Output: 
      !           x: the first n columns are the lowest n eigenstates of h 
      !              in the reduced space
      !          hx: the first n columns are the application of h on the 
      !              lowest n eigenvectors.
      ! eigenvalues: eigenvalues of the s.p. hamiltonian in the reduced space
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Technical notes:
      !  This routine accepts spinors on the mesh in the format
      !      x(4*mv, m)
      !  while the rest of the code operates
      !      x(mv, 4,m)
      !  The reason is to (i) make the Lapack calls more transparent and (ii)
      !  to aid compiler vectorisation.
      !------------------------------------------------------------------------

      integer, intent(in)          ::  m, n, sxb(4), syb(4), szb(4), iso
      real(KIND=dp), intent(inout) ::  x(mv*4,m)
      real(KIND=dp), intent(inout) :: hx(mv*4,m), upd(mv*4,m)
      real(KIND=dp), intent(out)   :: eigenvalues(n)
      real(KIND=dp)                :: sph(m,m), temp(4*mv,m), tempe(m)

      integer                      :: lwork, info
      real(KIND=dp), allocatable   :: work(:)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Populate the full matrix sph 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! This BLAS call does not exploit the symmetry of the problem, it 
      ! explicitly computes all elements of sph even though we know that 
      !    sph(i,j) = sph(j,i)^*
      ! TODO: figure out whether this BLAS call outperforms a symmetric 
      !       implementation with many matrix-vector calls. 
      !
      call DGEMM('t', 'n', m, m, 4*mv, dv, x(:,:), 4*mv, &
      &                                   hx(:,:), 4*mv,0.0d0,sph,m)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Diagonalise the sphamiltonian
      ! 
      ! First inquire about working memory
      allocate(work(1))
      call DSYEV('V','L', m ,sph,m,tempe,work,-1,info)
      lwork=work(1)
      deallocate(work)
      allocate(work(lwork))
      ! ..... and now do the actual work
      call DSYEV('V','L', m ,sph,m,tempe,work,lwork,info)
      if(info.ne.0) then
        call stp('Issue with diagonalising in eval_sph.')
      endif
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Now we construct the n lowest eigenvectors
      temp = x! temporary copy
      call DGEMM('n','n',4*mv,n,m, 1.0d0,temp, 4*mv,sph(:,1:n), m, 0.0d0, & 
      &                                        x(:,1:n), 4*mv)
      ! ... and aply the same transformation to upd
      temp = upd ! temporary copy
      call DGEMM('n','n',4*mv,n,m, 1.0d0,temp, 4*mv,sph(:,1:n), m, 0.0d0, & 
      &                                        upd(:,1:n), 4*mv)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Bookkeeping
      eigenvalues = tempe(1:n)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine diag_sph

    subroutine IterativeEstimation(Iteration)
      !-------------------------------------------------------------------------
      ! Estimate optimum parameters (dt,mu) of the heavy-ball iterative process
      ! to try and achieve optimal convergence rate.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Input:
      !   Iteration : integer, outer iteration counter. Used to decide whether 
      !               to allocate things, yes or no. 
      !-------------------------------------------------------------------------

      use wavefunctions

      1 format (a20, 99f10.3)
      2 format ('-------------------------------------------------------------')
      3 format (' Warning: maximum eigenvalue of h could not be estimated.  ')
      4 format (' Isospin = ', i2, 'maxE = ', f10.3,  ' convergence =', es10.3)

      integer, intent(in)              :: iteration

      real(KIND=dp), allocatable, save :: maxspwf(:,:,:)
      real(KIND=dp), allocatable, save :: update(:,:), actionofh(:,:)
      real(KIND=dp), allocatable, save ::   dmax(:,:,:)
      real(KIND=dp), allocatable, save ::  ddmax(:,:,:)
      real(KIND=dp), allocatable, save :: dddmax(:,:,:)

      integer       :: estiter, iter, ii, i, it, iso
      real(KIND=dp) :: con(2), maxE, compare, relE, kappa, Es(2)
      !-------------------------------------------------------------------------
      ! Step 1: Solve the auxiliary problem for the largest single-particle 
      !         energy on the mesh
      if(Iteration .eq.1) then

          if(allocated(maxspwf))   deallocate(maxspwf)
          if(allocated(update))    deallocate(update)
          if(allocated(actionofh)) deallocate(actionofh)
          if(allocated(dmax))      deallocate(dmax)
          if(allocated(ddmax))     deallocate(ddmax)
          if(allocated(dddmax))    deallocate(dddmax)

          ! Initialize with a random spwf at the start.
          allocate(maxspwf(nx*ny*nz,4,2)) 
          allocate(update(nx*ny*nz,4)) ; allocate(actionofh(nx*ny*nz,4))
          allocate(dmax(nx*ny*nz,3,4))
          allocate(ddmax(nx*ny*nz,6,4))
          allocate(dddmax(nx*ny*nz,10,4))

          call random_number(maxspwf)                        ! randomize
          do it=1,2
            maxspwf(:,:,it) = &                                      ! normalize
                        & 1.0/sqrt(sum(maxspwf(:,:,it)**2)*dv) * maxspwf(:,:,it)
          enddo
      endif

      estiter = 500
      update  = 0.0
      maxE    = 100.0 ! Initialize some value to avoid compiler complaints
      Es      = 100.0
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Iterative estimation of the maximal energy: evolve two single-particle
      ! wavefunctions (one for each isospin) to guess at the maximal eigenvalue
      ! of the single-particle hamiltonian.
      Es(2) = 0.0d0
      do it = 1,2
        con(it) = 1
        iso     = 2*it-3
        do iter=1,estiter
          !---------------------------------------------------------------------
          ! Two notes on this call to sphamil
          ! - onthefly = .true., such that derivatives are calculated
          ! - sx/y/z_max are set in the set_spwf_symmetries routine and are
          !   assumed to be the reflection quantum numbers of the very first
          !   symmetry block.
          actionofh = sphamil(maxspwf(:,:,it), dmax, ddmax,                    &
$N3       &                                        dddmax,                     &
          &                                     sx_max,sy_max,sz_max,iso,.true.)
          con(it)   = Es(it)
          Es(it)    = sum(actionofh * maxspwf(:,:,it)) * dv
          con(it)   = con(it) - Es(it)
          !---------------------------------------------------------------------
          ! Simple power iteration seems to better than gradient descent          
          maxspwf(:,:,it)  = actionofh 
          !---------------------------------------------------------------------
          ! Normalize
          maxspwf(:,:,it) = 1.0/sqrt(sum(maxspwf(:,:,it)**2)*dv)*maxspwf(:,:,it)
          !---------------------------------------------------------------------
          ! Don't be to picky about convergence, within the order of an MeV is
          ! good enough.
          if(abs(con(it)).lt. 1d-2) exit
        enddo
      enddo
      
      if(any(abs(con) .gt. 1d-2)) then
          print 1
          print 2
          print 3
          print 4, -1, Es(1), con(1)
          print 4, +1, Es(2), con(2)
      endif
      ! Take the maximum value of both isospins
      maxE = maxval(Es)
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
          relE = minval(abs(QPenergies))
      end select
      ! Safeguard the difference
      if(relE .lt. 0.20) relE = 0.50
      !-------------------------------------------------------------------------
      ! Step 3: use these estimations to determine a value for dt and mu.
      maxE  = maxE - minval(spenergies)
      kappa = relE/maxE
      
      if(Iteration.eq.1) then
        ! Don't mess up with a too large step at the start of the iterations
        dt = 2.0/maxE*hbar * stepsize_safety
        momentum = 0.0
      else        
        momentum = ((sqrt(kappa) - 1)/(sqrt(kappa)+1))**2
        dt    = 4.0/(maxE+relE+2*sqrt(maxE*relE))*hbar *  stepsize_safety
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
   real(KIND=dp)         :: mpsi(nx*ny*nz,4), jpsi(nx*ny*nz,4)
   real(KIND=dp)         :: O2, value, des, scale, crankfactor(3)
   integer               :: wave, k, B, si, N, it, i

   call start_timer(T_feasible)

   !----------------------------------------------------------------------------
   ! (i) The contribution of the multipole moments to the update
   Current    => Root
   multipole = 0.0_dp
   call compcutoff()
   
   do while(associated(Current%Next))
    Current => Current%next
   
    if(Current%ConstraintType.lt.2) cycle

    select case(Current%isoswitch)
    case(0)
      Value = sum(Current%Value)                    ! Total value
      O2    = sum(Current%Squared)                  ! < C^2 >
    case(1,2)
      it    = Current%isoswitch
      Value = Current%Value(it)                     
      O2    = Current%Squared(it)                  ! < C^2 >
    end select
    Des   = Current%Constraint                      ! Desired final value
    scale = Current%Scalefactor                     ! Scale factor
    
    update = 0.0
    select case(Current%isoswitch)
    case(0)
      do it=1,2
        Update(:,it) = 0.5*(Value-Des)/O2*Cutoff(:,it)*Current%SpherHarm*scale
      enddo
    case(1,2)
      ! only one nucleon species feels the constraint
      it = Current%isoswitch
      Update(:,it) = 0.5*(Value-Des)/O2*Cutoff(:,it)*Current%SpherHarm*scale
    end select
    multipole = multipole + Update
   enddo
   
   !----------------------------------------------------------------------------
   ! (ii) The contribution of the cranking constraints to the update
   do i=1,3
     if(CrankType(i).ne.1) cycle ! Only include cranking for cranktype=1
     CrankFactor(i)= 0.5*(TotalAngMom(i)-CrankValues(i))/J2_sp(i)
     ! Rescale with a factor
     Crankfactor(i) = Crankfactor(i)*CrankScaleFactor(i)
   enddo
   
   !----------------------------------------------------------------------------
   ! We evolve the spwfs
   si = 0   
   do B=1,8
    N = HFBlocks(B) ; if(N.eq.0) cycle
    it = 1 ;  if(B.gt.4) it = 2

    do wave=1,N
      do k=1,4
        mpsi(:,k) = multipole(:,it) * HFPsi(:,k,si+wave)
      enddo
      
      jpsi = 0.0d0
      do i=1,3
        if(cranktype(i) .ne. 1) cycle
        ! Add J_i | psi >
        jpsi = jpsi + crankfactor(i) &
        &          * AngMomOperator(HFpsi(:,:,si+wave), HFdpsi(:,:,:,si+wave),i)
      enddo
      do k=1,4
        jpsi(:,k) = cutoff(:,it) * jpsi(:,k)      
      enddo
      
      ! Substituting the correction
      HFPsi(:,:,si+wave) = HFPsi(:,:,si+wave) - mpsi - jpsi
 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
      !HFPsi(:,:,si+wave) = HFPsi(:,:,si+wave) - 2 * mpsi
      ! The factor two is a historical accident, and could be of course 
      ! accomodated by a redefinition of the Update above, but I prefer to 
      ! include it here and leave a trace of this happy (?) mistake.
      
      ! 22/07/21: turns out the factor two was not a happy mistake. For 
      ! quadrupole constraints in EV8/CR8-mode, the code worked fine. For 
      ! EV4-like calculations, this turned out to be too aggressive.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
    enddo
    si = si + N
   enddo
   !---------------------------------------------------------------------------
   ! Finally, orthonormalisation
   call orthonormalize

   call stop_timer(T_feasible)

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

    subroutine clean_evolution()
      if(allocated(preconx)) deallocate(preconx)
      if(allocated(precony)) deallocate(precony)
      if(allocated(preconz)) deallocate(preconz)
      if(allocated(momentum_updates)) deallocate(momentum_updates)
    end subroutine clean_evolution
end module evolution
