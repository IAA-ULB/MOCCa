module HFB_gradient
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
 ! This module implements all routines needed to gradient-step our way to the
 ! solution of the HFB problem. The main "driving" routine is located HFB.f90.
 !
 ! Subroutines
 ! -------------------
 ! - subroutine gradient_step
 ! - subroutine find_fermi_brent
 ! - subroutine Brent_bisection
 ! - subroutine ortho_bogo
 ! - subroutine diagonalise_H11H20
 ! - occupation_blocks
 ! - subroutine diag_by_occ_block
 ! - symmetrize_bogo
 ! - diagonalise_H_free
 ! - occupation_blocks
 ! - permute_columns
 !
 ! Functions
 ! -------------------
 ! - buildgrad
 ! - precongrad
 ! - GradUpdate
 ! - calcH20
 ! - calcH11
 ! - calcN20
 ! - calcN11
 ! - disp_bogo
 ! - particle_number_bogo
 !
 !==============================================================================
 ! TODO: 
 ! 
 ! a. the overall combined workings of HFB.f90 and HFB_gradient.f90 are
 !    more complicated than they need to be for a historical reason. 
 ! 
 !    Historically, I implemented the gradient solver only for the trivial
 !    occupation factors: f_mu=0. The current implementation thus juggles the 
 !    columns in the Bogoliubov matrix, configuration matrix and other things, 
 !    resulting in great confusion.
 ! 
 !    However, this is not necessary: the new and improved gradient solver
 !    accepts now f_mu different from zero. If it works correctly, this module
 !    should then no longer care about the ordering of quasiparticle columns
 !    or anything. This will eliminate quite some index juggling and the 
 !    whole notion of 'grad_blocks'.
 !
 !
 !==============================================================================
 
  use geninfo
  use wavefunctions

  implicit none
  
  !-----------------------------------------------------------------------------
  ! Parameters for the heavy-ball evolution in the pairing subproblem.
  real(KIND=dp) :: gradient_stepsize = 0.03, gradient_mu = 0.6
  real(KIND=dp) :: gradient_safety   = 0.1
  !-----------------------------------------------------------------------------
  ! The norm of the gradient for every isospin
  !
  !      N  = ||  H^{20} - lambda N^{20} || + || H^{11} - lambda N^{1} ||
  !
  ! where || is the Frobenius norm of the matrix and in the second part only
  ! the non-diagonal matrix elements are counted. 
  real(KIND=dp) :: HFBgradnorm(2) = 0
  !-----------------------------------------------------------------------------
  ! History of the gradients for the heavy-ball evolution. 
  ! Saved here, as the routines below can be used for either isospin
  real(KIND=dp), allocatable :: Z_updates(:,:,:)
  !-----------------------------------------------------------------------------
  ! Effective block size of the matrices in the gradient solver. Note that this
  ! does not NEED to be equal to the HFBlocks, as with blocking we can change
  ! the relative balance of signature +/- states in the U & V columns.  
  integer, save :: grad_blocks(8) = 0
  !-----------------------------------------------------------------------------
  ! To use (or not) preconditioning of the gradient in this module.
  logical :: gradient_precon = .false.
  !-----------------------------------------------------------------------------
  ! Allow the code to estimate evolution parameters by itself
  logical :: Estimategradparams = .true.

contains 

  pure function buildgrad(H20,N20,H11,N11,occ,lambda,Eqp,precon) result(grad)
    !---------------------------------------------------------------------------
    ! Construct the gradient of the energy. 
    !
    ! The "bare" gradient is
    !     g  = ( g^11  g^20 )
    !
    ! with 
    !
    !     g^11  = [ H^{11} - lambda N^{11}, f]
    !     g^20  = ( H^{20} - lambda N^{20} ) - { f, H^{20} - lambda N^{20}}
    !
    ! where
    !
    !    {,} is an anticommutator and f is the matrix of qp occupations.
    !
    ! The 20-component of the gradient can be preconditioned by employing an 
    ! approximation to the second derivative of the energy:
    !
    !     [Pg]_mn = g_mn/max(E^qp_m + E^qp_n , 2.0)
    ! 
    ! where this equation should only be used in the quasi-particle basis, i.e.
    ! the basis that diagonalizes H^{11}. 
    !
    ! Input: 
    !   H20   : matrix elements of the two-qp part of the HFB Hamiltonian
    !   N20   : matrix elements of the two-qp part of the number operator
    !   lambda: Fermi energy of the species under consideration
    !   occ   : occupation factors f
    !   precon: logical, indicating whether to precondition or not
    !
    ! Output: 
    !   grad  : matrix elements of the gradient (Z-matrix), built out of the 
    !           ingredients passed in.
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable :: grad(:,:,:)
    real(KIND=dp), intent(in)  :: H20(:,:), N20(:,:), H11(:,:), N11(:,:)
    real(KIND=dp), intent(in)  :: occ(:), Eqp(:), lambda
    logical, intent(in)        :: precon 
    integer                    :: i,j, N

    N = size(H20,1)
    allocate(grad(N,N,2)) 

    grad(:,:,1) = H11 - lambda * N11
    grad(:,:,2) = H20 - lambda * N20
    
    do i=1, N
      do j=1, N
        grad(i,j,1) = grad(i,j,1) * (         occ(j) - occ(i))
        grad(i,j,2) = grad(i,j,2) * ( 1.0d0 - occ(j) - occ(i))
      enddo
    enddo
    
    if(precon) grad = precon_grad(grad, Eqp)
  end function buildgrad

  pure function precon_grad(grad, Eqp)  result(Pgrad)
    !---------------------------------------------------------------------------
    ! Precondition the gradient with estimates of the quasiparticle energies.
    !
    !    [Pg]^11_mn = g^11_mn/max(E^qp_m - E^qp_n , 2.0)
    !    [Pg]^20_mn = g^20_mn/max(E^qp_m + E^qp_n , 2.0)
    !
    ! Input: 
    !   grad : gradient to be preconditioned
    !   Eqp  : quasiparticles to use in the preconditioner
    !
    ! Output:
    !   Pgrad: preconditioned gradient
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! See:
    ! * L. M. Robledo and G. F. Bertsch, 
    !   Application of the gradient method to Hartree-Fock-Bogoliubov theory. 
    !   Physical Review C, 84(1), 014312 (2011).  
    !   https://doi.org/10.1103/PhysRevC.84.014312
    !
    ! * S. Perez-Martin, 
    !   Microscopic description of atomic nuclei with an odd number of nucleons. 
    !   PhD Thesis, Universidat Autonoma de Madrid (2006).
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: grad(:,:,:), Eqp(:)
    real(KIND=dp)              :: fac
    real(KIND=dp), allocatable :: Pgrad(:,:,:)
    integer   :: i,j

    Pgrad = grad
    do i=1, size(grad,1)
        do j=1, size(grad,1)
          ! 11-part
          fac = max(abs(Eqp(i) - Eqp(j)), 2.0d0)
          Pgrad(i,j,1) = Pgrad(i,j,1)/fac
          ! 20-part
          fac = max(Eqp(i)+ Eqp(j), 2.0d0)
          Pgrad(i,j,2) = Pgrad(i,j,2)/fac
        enddo
    enddo

  end function precon_grad

  subroutine gradient_step(h,gaps,targetN, Bogo, occ, Eqp,lambda, alpha,mu, &
  &                        prev, precon, gradnorm, blocks,   maxiter, ifail)
    !---------------------------------------------------------------------------
    ! Perform one (or more) heavy-ball evolution steps, starting from an initial
    ! Bogoliubov transformation with given single-particle hamiltonian h and 
    ! pairing gaps Delta. 
    !
    !  -> 1. Find an appropriate Fermi energy to fix the particle number
    !  |  2. Update U and V matrices with a heavy-ball Thouless step
    !  |  3. Update the momentum array for the next iteration
    !  <- 4. Check for convergence or maximum number of iterations
    !     5. Symmetrize the Bogoliubov transformation. 
    !     6. Diagonalise the HFB Hamiltonian in every occupation block, to 
    !        obtain quasiparticle energies.
    !    
    ! Assuming a specific formatting of the Bogoliubov transformation:
    !
    !              (  V^*   U   )
    !       B  =   (            )
    !              (  U^*   V   )
    !
    ! this module ONLY evolves the right-most half of the Bogoliuv matrix in 
    ! steps 1-4, such that we construct the left-hand part of the Bogoliubov
    ! matrix explicitly in step 5.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    ! Input:
    !---------
    ! h         : single-particle hamiltonian 
    ! gaps      : pairing gaps
    ! targetN   : targetted number of particles
    ! Bogo      : Bogoliubov transformation on input
    ! occ       : occupation factors f
    ! Eqp       : an estimate for the quasi-particle energies for the 
    !             preconditioning of the evolution (if requested)
    ! lambda    : Fermi energy
    ! alpha     : step-size for the heavy-ball evolution
    ! mu        : momentum for the heavy-ball evolution
    ! precon    : whether or not to precondition the update
    ! blocks    : size of the quantum number blocks for this evolution 
    ! maxiter   : number of steps to do (recommended: 1)
    !
    ! Output:
    !---------
    !    Bogo : new Bogoliubov transformation
    !    Eqp  : estimated qp energies based on the diagonalisation of  H^{11}
    !   lambda: new value for the Fermi energy
    ! gradnorm: Frobenius norm of the gradient that was used as step
    !  ifail  : if 0, succes. If 1, something went wrong.
    !---------------------------------------------------------------------------  
    integer, intent(in)          :: blocks(4), maxiter
    integer, intent(out)         :: ifail
    real(KIND=dp), intent(in)    :: targetN, h(:,:), gaps(:,:), alpha, mu
    real(KIND=dp), intent(in)    :: occ(:)
    logical, intent(in)          :: precon
    real(KIND=dp), intent(inout) :: Bogo(:,:), lambda
    real(KIND=dp), intent(inout) :: Eqp(:),  prev(:,:,:), gradnorm
    real(KIND=dp), allocatable   :: H20(:,:), N20(:,:), grad(:,:,:)
    real(KIND=dp), allocatable   :: N11(:,:), H11(:,:)

    real(KIND=dp)                :: particles,  normN, lambda_corr
    integer                      :: iter
    logical                      :: converged
  
    converged = .false. 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! If we want to precondition the evolution, we should attempt to diagonalise
    ! as much as possible the HFB Hamiltonian without affecting the many-body 
    ! state, such that we get estimates for the quasiparticle energies.
    if(precon) then
      H11 = calcH11(Bogo, h, gaps, blocks) ; N11 = calcN11(Bogo, blocks)
      H20 = calcH20(Bogo, h, gaps, blocks) ; N20 = calcN20(Bogo, blocks)
      call diagonalise_H_free(Bogo,prev,H11,N11,H20,N20,lambda,occ,blocks,Eqp)
    endif
    
    do iter=1,maxiter
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! We calculate the relevant matrices to build the gradient 
        H20 = calcH20(Bogo, h,gaps,blocks) ; H11 = calcH11(Bogo, h, gaps, blocks)
        N20 = calcN20(Bogo, blocks)        ; N11 = calcN11(Bogo, blocks) 
        
        ! First, we check if the pairing has collapsed
        normN     = calc_gradN(N11,N20,occ)   
        if(normN.gt.1d-4) then
          ! 1. If the pairing has not collapsed, we try different values of the 
          ! Fermi energy to arrive at a correct particle number AFTER the 
          ! heavy-vall step. 
          call find_fermi_brent(Bogo, H20, N20, H11, N11, prev, occ, Eqp,   &
          &        lambda, particles, targetN, alpha, mu, precon, blocks, ifail)
        endif
        
        ! 2. Building the gradient update 
        ! (with the old Fermi energy if the pairing collapsed)
        grad = buildgrad(H20, N20, H11,N11, occ, lambda, Eqp, precon)
        ! ...  and executing it
        bogo = GradUpdate(grad, prev, bogo, alpha, mu, blocks)
        ! 3. ...  and saving it for the next iteration 
        prev = -alpha * grad + mu * prev
        ! ... we recalculate the deviation of the particle number
        particles = particle_number_bogo(bogo, occ, blocks) - targetN
        ! .... as well as the norm of the gradient 
        gradnorm  = sqrt(sum(grad**2))
        N20 = calcN20(Bogo, blocks) ; N11 = calcN11(Bogo, blocks)
        normN     = calc_gradN(N11,N20,occ)
        if(abs(particles).gt.pairing_prec .and. normN .gt. 1d-4) then
          N20 = calcN20(Bogo, blocks) ; N11 = calcN11(Bogo, blocks)
          H20 = 0.0                   ; H11 = 0.0
          lambda_corr = 0.1
          call find_fermi_brent(Bogo, H20, N20, H11, N11, prev, occ, Eqp,   &
          &                     lambda_corr,particles, targetN, 1.0d0, 0.0d0,  &
          &                     precon, blocks, ifail)
          grad = buildgrad(H20, N20, H11, N11, occ, lambda_corr, Eqp, precon)
          bogo = GradUpdate(grad, prev, bogo, 1.0d0, 0.0d0, blocks)
        endif

        ! 4. Check for convergence if this is process is repeated multiple times
        if(gradnorm .lt. 1d-6) converged = .true.

        if(converged) then
          exit
        endif
    enddo
    !---------------------------------------------------------------------------
    ! 5. Construct the left-most half of the Bogoliubov matrix by symmetry
    call symmetrize_bogo(bogo, blocks)
    !---------------------------------------------------------------------------
    ! 6. Perform an additional transformation of the Bogoliubov transformation 
    ! and the prev array: bring all matrix elements of the HFB Hamiltonian on 
    ! which the energy does not depend to zero. This also gives us an estimate 
    ! of the quasiparticle energies.
    H11 = calcH11(Bogo, h, gaps, blocks) ; N11 = calcN11(Bogo, blocks)
    H20 = calcH20(Bogo, h, gaps, blocks) ; N20 = calcN20(Bogo, blocks)
    call diagonalise_H_free(Bogo,prev,H11,N11,H20,N20,lambda,occ,blocks,Eqp)

  end subroutine gradient_step
  
  pure function calc_gradN(N11,N20,occ) result(norm)
      !-------------------------------------------------------------------------
      ! Calculate the total norm of the gradient in the N direction.
      !
      ! Input:
      !    N11 : 11-component of the particle number operator
      !    N20 : 20-component of the particle number operator
      !    occ : quasiparticle occupation factors
      !
      ! Output:
      !    norm: norm of the gradient in the N direction
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in) :: N11(:,:), N20(:,:), occ(:) 
      integer                   :: i, j, N
      real(KIND=dp)             :: norm
      
      norm = 0.0d0
      
      N = size(occ)
      do i=1,N
        do j=1,N
          norm = norm + (N20(i,j) * (1 - occ(i) - occ(j)))**2 &
          &           + (N11(i,j) * (occ(j) - occ(i)))**2
        enddo
      enddo  
      norm = sqrt(norm)  
      
  end function calc_gradN

  pure subroutine find_fermi_brent(Bogo, H20, N20, H11, N11, prev, occ, Eqp,&
  &                               lambda,particles, targetN, alpha, mu, precon,&
  &                               blocks, ifail)
    !---------------------------------------------------------------------------
    ! Find a Fermi energy such that the particle number is (on average) correct
    ! AFTER the heavy-ball evolution, using Brents method:
    ! 
    ! https://en.wikipedia.org/wiki/Brent%27s_method
    ! 
    ! which combines bisection, secant method and inverse quadratic 
    ! interpolation.The original source is probably
    ! R. P. Brent (1973), "Chapter 4: An Algorithm with Guaranteed Convergence
    ! for Finding a Zero of a Function", Algorithms for Minimization without
    ! Derivatives, Englewood Cliffs, NJ: Prentice-Hall.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input :  
    ! ----------
    !   Bogo     : Bogoliubov transformation. Only the right-most half is 
    !              referenced.
    !   H20,H11  : Quasiparticle components of the HFB Hamiltonian.
    !   N20,N11  : Quasiparticle components of the particle number operator.
    !   prev     : Previous update of the U,V matrices
    !   occ      : occupation factors of the quasiparticles
    !   Eqp      : quasiparticle energies, only used when preconditioning
    !   lambda   : Initial guess for the Fermi energy
    !   targetN  : Targetted number of particles
    !   alpha    : stepsize for the heavy-ball step
    !   mu       : momentum parameter for the heavy-ball step
    !   precon   : whether or not to precondition
    !   blocks   : structure of the symmetry blocks
    !
    ! Output:
    ! ----------
    !   lambda   : Fermi energy that leads to the targetted particle number
    !              after completion of the heavy-ball step
    !   particles: calculated number of particles
    !   ifail    : error code for the initial bracketing phase
    !              0 : succesfull
    !              10: bracketing failed 
    !  
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! The original routine for the direct HFB solver is very heavily 
    ! inspired/copy-pasted by the routines implemented in MOCCa by M. Bender.
    ! This routine, in turn, is almost a copy of that one. 
    !-------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: H20(:,:), N20(:,:), H11(:,:), N11(:,:)
    real(KIND=dp), intent(in)    :: targetN, Eqp(:)
    integer, intent(in)          :: blocks(4)
    real(KIND=dp), intent(in)    :: alpha, mu, prev(:,:,:), occ(:)
    logical, intent(in)          :: precon
    integer, intent(out)         :: ifail

    real(KIND=dp), intent(out)   :: particles
    real(KIND=dp), intent(inout) :: lambda, bogo(:,:)
    real(KIND=dp), allocatable   :: gradA(:,:,:), gradB(:,:,:)
    real(KIND=dp), allocatable   :: nbA(:,:), nbB(:,:)

    real(KIND=dp)                :: InitialBracket(2), FA, FB, N
    integer                      :: idir, idirsig, FailCount
    logical                      :: Success

    idir = 0 ; idirsig = 1

    !---------------------------------------------------------------------------
    ! STEP 1: set up an initial bracket
    !---------------------------------------------------------------------------
    gradA = buildgrad(H20, N20, H11, N11, occ, lambda, Eqp, precon)
    nbA   = GradUpdate(gradA, prev, bogo, alpha, mu, blocks)
    N     = particle_number_bogo(nbA, occ, blocks) - targetN
    
    ! Check if this guess for lambda is good enough
    if(abs(N).lt.pairing_prec .or. alpha .eq. 0.0d0) then
      return
    endif
    ! Use present Fermi energy as starting point and check the direction
    ! where the zero of <N>-N0 can be expected.
    ! If <N>-N0 <  0, search at higher values.
    ! If <N>-N0 >= 0, search at lower  values.
    ! Initialize InitialBracket(it,1) = A, InitialBracket(it,2) = B with 
    ! present  Fermi energy.
    ! "dir" is the label of the InitialBracket(it,idir) that has to be moved,
    ! "idirsig" is the sign of steps needed to go into that direction.
    InitialBracket(:) = lambda
    if ( N .lt. 0.0_dp ) then 
      idir   =  2 ;  idirsig =  1
    else 
      idir   =  1 ;  idirsig = -1
    endif
  
    ! Try to find a boundary that brackets the Fermi energy in the direction 
    ! into which the Fermi energy has to be changed.
    FailCount = -1 ;  Success = .false.

    do while(.not. Success)
        FailCount = FailCount + 1
        
        ! update moving boundary and recalculate particle numbers at both.
        InitialBracket(idir) = &
        &                 InitialBracket(idir) + idirsig * 0.1_dp*(FailCount+1)

        gradA = buildgrad(H20,N20,H11,N11,occ,InitialBracket(1),Eqp,precon)
        gradB = buildgrad(H20,N20,H11,N11,occ,InitialBracket(2),Eqp,precon)

        nbA = GradUpdate(gradA, prev, bogo, alpha, mu, blocks)
        nbB = GradUpdate(gradB, prev, bogo, alpha, mu, blocks)

        FA = particle_number_bogo(nbA, occ, blocks) - targetN
        FB = particle_number_bogo(nbB, occ, blocks) - targetN
        
        ! diagnostic printing for convergence analysis (usually commented out)
        !print '(" Bracketing ",i4,1l2,(2(f13.8,es16.7)))',        &
        !      & FailCount,Success,InitialBracket(1),FA, InitialBracket(2),FB          
        ! code failure (Fermi energy has changed by 30 MeV)
        if (Failcount .gt. 76) then
!          print '(/," A = ", f13.8, " FA = ",1es12.4,              &
!               &    " B = ", f13.8, " FB = ",1es12.4)',            &
!               &     InitialBracket(1),FA,InitialBracket(2),FB 
          ifail = 10
          return
          !stop 'FindFermiBrent: Search for InitialBracket failed.'
        endif
        ! check if root is bracketed for isospin it after the update
        if( FA*FB .lt. 0.0_dp ) then 
            ! Correct Bracket found!
            Success = .true.
        endif
      enddo

      call Brent_bisection(initialbracket(1),initialbracket(2),FA,FB,bogo,     & 
      &                    H20, N20, H11, N11, occ, prev, Eqp, lambda,      &
      &                    particles, targetN, alpha, mu, precon, blocks,200)

  end subroutine find_fermi_brent

  pure subroutine Brent_bisection(X1,X2,FX1, FX2, Bogo, H20, N20, H11, N11,   &
  &                              occ, prev, Eqp, lambda,particles, targetN, &
  &                              alpha,  mu, precon, blocks, depth)
    !---------------------------------------------------------------------------
    ! This routine searches for the Fermi energy
    ! by Brent's methods https://en.wikipedia.org/wiki/Brent%27s_method
    ! which combines bisection, secant method and inverse quadratic 
    ! interpolation. The original source is probably
    ! R. P. Brent (1973), "Chapter 4: An Algorithm with Guaranteed Convergence
    ! for Finding a Zero of a Function", Algorithms for Minimization without
    ! Derivatives, Englewood Cliffs, NJ: Prentice-Hall, 
    !---------------------------------------------------------------------------
    ! see pages 1188 - 1189 of http://apps.nrbook.com/fortran/index.html
    ! W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery,
    ! Numerical Recipes in Fortran in Fortran 90, Second Edition (1996).
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: H20(:,:), N20(:,:), H11(:,:), N11(:,:)
    real(KIND=dp), intent(in)    :: alpha, mu, targetN, Eqp(:)
    real(KIND=dp), intent(in)    :: X1 , X2, FX1 , FX2, prev(:,:,:), occ(:)
    integer, intent(in)          :: blocks(4), depth
    logical, intent(in)          :: precon 

    real(KIND=dp), intent(out)   :: particles
    real(KIND=dp), intent(inout) :: lambda, bogo(:,:)
    
    real(KIND=dp), parameter     :: eps = 1.d-9
    
    real(KIND=dp), allocatable   :: grad(:,:,:), newbogo(:,:)
    real(KIND=dp)                :: A , B, C , FA, FB , FC
    real(KIND=dp)                :: D , E, S , P  , Q , R 
    real(KIND=dp)                :: Num , Tol , XM 
    integer                      :: FailCount
    logical                      :: Found

    A  = X1 ; B  = X2 
    FA = FX1; FB = FX2
    Found = .false.    
    if (A .eq. B) then 
      !-------------------------------------------------------------------------
      ! This signals that FA = FB is zero within the tolerance.
      ! Either near-converged HFB or HF case of completely broken-down pairing
      ! which also satisfies FA = FB = 0 within an interval. The 
      ! latter case cannot be handled by the algorithm below.
      !-------------------------------------------------------------------------
      Found = .true.
    endif

    C = B ; FC = FB 
    E = -1000000 ; D = -1000000
  
    FailCount = -1

    do while(.not.Found) 
      FailCount = FailCount + 1
      if ( ( FB .gt. 0.0_dp .and. FC .gt. 0.0_dp ) .or. & 
         & ( FB .lt. 0.0_dp .and. FC .lt. 0.0_dp ) )  then
        C  = A     ;  FC = FA
        D  = B - A ;  E  = D
      endif
      if ( abs(FC) .lt. abs(FB) ) then
        A  = B ;  FA = FB
        B  = C ;  FB = FC
        C  = A ;  FC = FA
      endif
      !-------------------------------------------------------------------------
      ! Convergence check
      ! Note (W.R.): I have tightened convergence a bit compared to the values
      !              in MOCCa by M.B. 
      !-------------------------------------------------------------------------
      Tol  = 2.0_dp * eps * abs(B) + 0.05_dp * Pairing_prec
      XM   = 0.5_dp * (C-B)
      !----------------------------------------------------------------
      ! Note: the tolerance is on the precision of the Fermi energy,
      ! NOT the nearness of the particle number to the targeted value.
      !----------------------------------------------------------------
      if ( abs(XM) .le. Tol .or. FB .eq. 0.0_dp ) then
        Lambda =  B
        Found = .true. 
        cycle
      endif
      if ( abs(E) .ge. Tol .and. abs(FA) .gt. abs(FB) ) then
        S = FB/FA
        if ( A .eq. C ) then
          P = 2.0_dp * XM * S
          Q = 1.0_dp - S
        else
          Q = FA/FC
          R = FB/FC
          P = S * (2.0_dp * XM * Q * (Q-R) & 
                  &    - (B-A)*(R-1.0_dp))
          Q = (Q-1.0_dp)*(R-1.0_dp)*(S-1.0_dp)
        endif
        if ( P .gt. 0.0_dp ) Q = -Q
        P = abs(P)
        if (2.0_dp * P .lt. min(3.0_dp*XM*Q - abs(Tol*Q),abs(E*Q))) then
          E = D
          D = P / Q
        else
          D = XM
          E = D 
        endif
      else
        D = XM
        E = D 
      endif
      A  = B 
      FA = FB
      B  = B + merge(D,sign(Tol,XM),abs(D) .gt. Tol)    
  
      !-------------------------------------------------------------------------
      ! B is present best guess for the fermi energy, FB the corresponding 
      ! particle number.
      !-------------------------------------------------------------------------
      grad = buildgrad(H20, N20, H11, N11, occ, B, Eqp, precon)
      newbogo = GradUpdate(grad, prev, bogo, alpha, mu, blocks)
      Num  = particle_number_bogo(newbogo, occ, blocks)
      FB   = Num - targetN
      particles = FB

      !-------------------------------------------------------------------------
      ! diagnostic printing for convergence analysis (usually commented out)
      !-------------------------------------------------------------------------
      ! NOTE: B is the the best guess for the zero of F. A has been the previous
      ! "closest" interval boundary that is not updated after Found
      ! is set to .true. The actual zero might therefore be outside the 
      ! interval [A,B]. If so, the true zero is typically closer to B than 
      ! A is to B.
      !-------------------------------------------------------------------------
      ! Note further: as A and B are swapped from time to time, B might be 
      ! smaller than A when printed here
      !-------------------------------------------------------------------------
!      print '(" BrentBisection ",i4,(1l2,2(f13.8,es16.7),f14.8))',    &
!           & FailCount, Found,A,FA,B,FB,Num
!     
      if ( FailCount .gt. Depth ) then
        Lambda    = B ; particles = FB 
        return 
      endif
    enddo
    ! Output
    Lambda    = B ; particles = FB 
  end subroutine Brent_bisection

  pure function GradUpdate(grad, prev, bogo, alpha, mu, blocks) result(NB)
    !---------------------------------------------------------------------------
    ! Update the bogoliubov transformation with a heavy-ball step:
    !
    !     U'  =  U - alpha * U * Z^11 + mu * dU^11  [11-update]
    !              - alpha * V * Z^20 + mu * dU^20  [20-update]
    !  
    !     V'  =  V - alpha * V * Z^11 + mu * dV^11  [11-update]
    !              - alpha * U * Z^20 + mu * dV^20  [20-update]
    ! 
    ! followed by an orthonormalisation via a Gramm-Schmidt algorithm.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    ! Input:
    !------------
    !     alpha   : step size
    !     mu      : momentum of the evolution 
    !     grad    : direction of the update, Z
    !               divided in 11 and 20 part
    !     prev    : previous update, (dU,dV)^T
    !               divided in 11 and 20 part
    !     bogo    : current value of the Bogoliubov transformation
    !     blocks  : structure of the symmetry-blocks
    !
    ! Output:
    ! ---------
    !     NB      : updated Bogoliubov transformation
    ! 
    !  Note that this routine only deals with the right-most half part of the 
    !  Bogoliubov transformation: 
    ! 
    !    B_complete = ( V^*  U )    B_here = ( U )
    !                 ( U^*  V )             ( V )
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: alpha, mu 
    real(KIND=dp), intent(in)    :: grad(:,:,:), prev(:,:,:)
    integer, intent(in)          :: blocks(4)
    real(KIND=dp), intent(in   ) :: bogo(:,:)
    real(KIND=dp), allocatable   :: V(:,:), U(:,:), NB(:,:)
    integer                      :: B, sb, N, N2, T, si

    NB = bogo
    sb = 0  ; si = 0
    do B=1,4,2
      N  = blocks(B)   ; if(N.eq.0) cycle
      N2 = blocks(B+1)
      T = N + N2
          
      ! Extract U and V matrices to make things more legible
      U = bogo(sb+  1:sb+  T,sb+T+1:sb+2*T)
      V = bogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! 11-component updates
      ! U'  = U - alpha U Z^11 + mu dU 
      NB(sb  +1:sb+T ,sb+T+1:sb+2*T) = NB(sb  +1:sb+T ,sb+T+1:sb+2*T)          &
      &                       - alpha * matmul(U,grad(si+1:si+T,si+1:si+T,1))  &
      &                       + mu    * matmul(U,prev(si+1:si+T,si+1:si+T,1))
      ! V'  = V - alpha V Z^11 + mu dV 
      NB(sb+T+1:sb+2*T,sb+T+1:sb+2*T) = NB(sb+T+1:sb+2*T,sb+T+1:sb+2*T)        &
      &                       - alpha * matmul(V,grad(si+1:si+T,si+1:si+T,1))  &
      &                       + mu    * matmul(V,prev(si+1:si+T,si+1:si+T,1))
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! 20-component updates
      !
      ! U'  = U - alpha V Z^20 + mu dU 
      NB(sb  +1:sb+T ,sb+T+1:sb+2*T) = NB(sb  +1:sb+T ,sb+T+1:sb+2*T)          &
      &                       - alpha *matmul(V,grad(si+1:si+T,si+1:si+T,2))   &
      &                       + mu    *matmul(V,prev(si+1:si+T,si+1:si+T,2))

      ! V'  = V - alpha U Z^20  + mu dV
      ! Note the additional minus sign for time-reversal conserved calculations
      NB(sb+T+1:sb+2*T,sb+T+1:sb+2*T) = NB(sb+T+1:sb+2*T,sb+T+1:sb+2*T)        &
$NTR      &               -alpha * matmul(U,grad(si+1:si+T,si+1:si+T,2))       &
$TR       &               +alpha * matmul(U,grad(si+1:si+T,si+1:si+T,2))       &
$NTR      &               + mu   * matmul(U,prev(si+1:si+T,si+1:si+T,2))
$TR       &               - mu   * matmul(U,prev(si+1:si+T,si+1:si+T,2))

      si = si +   T
      sb = sb + 2*T
    enddo
    ! - - - - - - - - - - - - - - - - -
    ! Don't forget to orthonormalise 
    call ortho_bogo(NB, blocks)
    
  end function GradUpdate

  pure subroutine ortho_bogo(bogo, blocks)
    !---------------------------------------------------------------------------
    ! Gramm-Schmidt routine to orthogonalise the Bogoliubov transformation. 
    ! This routine only MODIFIES the right-most half of the Bogoliubov 
    ! transformation. However, this Gramm-Schmidt procedure also ORTHOGONALIZES 
    ! against the (hypothetical) left-most half of the Bogoliubov transformation.
    !
    ! Input:
    !    Bogo  : Bogoliubov transformation
    !    blocks: Sizes of the symmetry blocks
    !
    ! Output: 
    !    Bogo : orthonormalized Bogoliubov transformation
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(inout) :: bogo(:,:)
    integer, intent(in)          :: blocks(:)
    integer                      :: i,j,k, N, N2, sb, T, B
    real(KIND=dp)                :: overlap
    real(KIND=dp), allocatable   :: partner(:)

    sb = 0
    do B=1,size(blocks),2
      N = blocks(B) ; if(N.eq.0) cycle
      N2= blocks(B+1)
      T = N+N2

      allocate(partner(2*T))

      do j=1,T

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Construct  the partner state ( V^*_j ) for readability
        !                              ( U^*_j )
$TR     partner(  1:  T) =-bogo(sb+T+1:sb+2*T,sb+T+j)
$NTR    partner(  1:  T) =+bogo(sb+T+1:sb+2*T,sb+T+j)
        partner(T+1:2*T) = bogo(sb  +1:sb+  T,sb+T+j)
      
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Normalise vector ( U_j )
        !                  ( V_j )
        overlap = 0.0_dp
        do i=1,2*T
            Overlap = Overlap + Bogo(sb+i,sb+T+j)**2
        enddo
        !print *, j,overlap
        Overlap    = 1.0_dp/sqrt(overlap)
        bogo(sb+1:sb+2*T,sb+T+j) = bogo(sb+1:sb+2*T,sb+T+j)*overlap
        
        overlap = 0.0_dp
        do i=1,2*T
            Overlap = Overlap + Bogo(sb+i,sb+T+j)*partner(i)
        enddo
        Bogo(sb+1:sb+2*T,sb+T+j) = Bogo(sb+1:sb+2*T,sb+T+j)    &
        &                        - Overlap * partner
        
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Orthogonalise all the rest against vector j
        do i=j+1,T
          Overlap = 0.0_dp
          do k=1,2*T
              Overlap = Overlap + bogo(sb+k,sb+T+j) * bogo(sb+k,sb+T+i)
          enddo
          !print *, i, overlap
          Bogo(sb+1:sb+2*T,sb+T+i) = Bogo(sb+1:sb+2*T,sb+T+i) &
          &                        - Overlap * Bogo(sb+1:sb+2*T,sb+T+j)
        enddo
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Orthogonalise all the rest against  ( V^*_j )
        !                                     ( U^*_j )
        do i=j+1,T
          Overlap = 0.0_dp
          do k=1,2*T
              Overlap = Overlap + partner(k) * bogo(sb+k,sb+T+i)
          enddo
          Bogo(sb+1:sb+2*T,sb+T+i) = Bogo(sb+1:sb+2*T,sb+T+i) &
          &                        - Overlap * partner
        enddo
      enddo
      sb = sb+ 2*T
      deallocate(partner)
    enddo
  end subroutine ortho_bogo

  pure function calcH20(Bogo,h,gaps,blocks) result(H20)
    !---------------------------------------------------------------------------
    ! Calculate the 2-quasi-particle-excitation component of H:
    ! 
    ! H20 =   U^{dagger} h   V^* - V^{\dagger} \Delta^* V^*
    !       - V^{dagger} h^t U^* + U^{\dagger} \Delta   U^*
    !
    ! Input:
    ! -------
    !   Bogo : Bogoliubov transformation (only right-half is used)
    !   h    : single-particle hamiltonian
    !   gaps : pairing gaps
    !   rho  : normal and anomalous density matrix 
    !          (only necessary when constraining the particle number dispersion)
    ! Output:
    ! -------
    !   H20   :  20-component of the HFB Hamiltonian
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)               :: H(:,:), bogo(:,:), gaps(:,:)
    integer, intent(in)                     :: blocks(4)
    
    real(KIND=dp), allocatable :: H20(:,:), U(:,:), V(:,:)
    real(KIND=dp), allocatable :: hV(:,:), hU(:,:), dV(:,:), dU(:,:)
!    real(KIND=dp), allocatable :: chi(:,:)
    integer                    :: B, N, N2, si, sb, T

    allocate(H20(sum(blocks), sum(blocks))) ; H20 = 0
    
    si = 0 ; sb = 0
    do B=1,4,2
      N  = blocks(B)    ; if(N.eq.0) cycle 
      N2 = blocks(B+1)
      T = N + N2

      ! Getting the U and V out to make the formulas explicit
      ! and the matrix multiplications memory-local
      U = Bogo(sb  +1:sb+  T,sb+T+1:sb+2*T)
      V = Bogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)

!      chi = -2*rho(si+1:si+T,si+1:si+T)
!      do i=1,T
!        chi(i,i) = chi(i,i) + 1
!      enddo

      !  h V^* and h^t U^*
      hV = matmul(   h(si+1:si+T,si+1:si+T), V) ! - 2*lambda2*matmul(chi, V)
      hU = matmul(   h(si+1:si+T,si+1:si+T), U) !- 2*lambda2*matmul(chi, U)
      
      ! d^* V^* and dU^*
      dV = matmul(gaps(si+1:si+T,si+1:si+T), V)
      dU = matmul(gaps(si+1:si+T,si+1:si+T), U)

      ! We reuse the defined symbols to save a matrix multiplication here
      U = transpose(U) ; V = transpose(V)

      ! We can save some effort here in the future, H20 is antisymmetric     
$NTR     H20(si+1:si+T, si+1:si+T)  = matmul(U, hV) + matmul(U, dU) &
$NTR                             &  - matmul(V, hU) - matmul(V, dV)
      ! Note the extra minus sign for time-reversal 
$TR      H20(si+1:si+T, si+1:si+T)  =-matmul(U, hV) + matmul(U, dU) &
$TR                              &  - matmul(V, hU) - matmul(V, dV)
      si = si +  T
      sb = sb +2*T
    enddo
  end function calcH20

  pure function calcH11(Bogo,h,gaps,blocks) result(H11)
    !---------------------------------------------------------------------------
    ! Calculate the 11 component of H 
    !  
    !   H^11 = U^T h U + U^T Delta V - V^T Delta U - V^t h V
    ! 
    ! Attention: the contribution of the Fermi energy is not included in this
    ! routine!
    !
    ! Input:
    ! -------
    !   Bogo  :  Bogoliubov transformation
    !   h     :  matrix elements of the single-particle hamiltonian
    !   gaps  :  matrix elements of the pairing gaps
    !   blocks:  structure of the symmetry blocks
    !
    ! Output:
    ! -------
    !   H11   :  11-component of the HFB Hamiltonian
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: h(:,:), bogo(:,:), gaps(:,:)
    real(KIND=dp), allocatable :: H11(:,:), U(:,:), V(:,:)
    real(KIND=dp), allocatable :: hV(:,:), hU(:,:), dV(:,:), dU(:,:)
    integer, intent(in)        :: blocks(4)
    integer                    :: B, N, N2, si, sb, T

    allocate(H11(sum(blocks), sum(blocks))) ; H11 = 0

    si = 0 ; sb = 0
    do B=1,4,2
      N  = blocks(B)    ; if(N.eq.0) cycle 
      N2 = blocks(B+1)
      T = N + N2

      ! Getting the U and V out to make the formulas explicit
      ! and the matrix multiplications memory-local
      U = Bogo(sb  +1:sb+  T,sb+T+1:sb+2*T)
      V = Bogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)

      !  h V^* and h U^*
      hV = matmul(   h(si+1:si+T,si+1:si+T), V) !- lambda * V
      hU = matmul(   h(si+1:si+T,si+1:si+T), U) !- lambda * U
      ! W.R. 19/07/22: N11 is now calculated separately from this routine

      ! d^* V^* and dU^*
      dV = matmul(gaps(si+1:si+T,si+1:si+T), V)
      dU = matmul(gaps(si+1:si+T,si+1:si+T), U)

      ! We reuse the defined symbols to save a matrix multiplication here
      U = transpose(U) ; V = transpose(V)

      ! We can save some effort here in the future, H20 is antisymmetric     
      H11(si+1:si+T, si+1:si+T)  = matmul(U, hU) &
                              &  + matmul(U, dV) &
$NTR                          &  - matmul(V, dU) & 
$TR                           &  + matmul(V, dU) & 
                              &  - matmul(V, hV)

      si = si +  T
      sb = sb +2*T
    enddo

  end function calcH11

  pure function calcN20(Bogo, blocks) result(N20)
    !---------------------------------------------------------------------------
    ! Calculate the 20-component of the particle number operator
    !
    !   N^{20} = U^dagger V^* - V^{\dagger} U^*
    !
    ! Input
    ! -------
    !  Bogo   : Bogoliubov transformation, only the right-most half is used
    !  blocks : structure of the symmetry blocks
    !
    ! Output:
    ! -------
    !  N20    : 20-component of the particle number operator
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: Bogo(:,:)
    real(KIND=dp), allocatable :: N20(:,:)
    real(KIND=dp), allocatable :: U(:,:), V(:,:)
    integer, intent(in)        :: blocks(4)

    integer :: B, si, N, N2, T, sb

    allocate(N20(sum(blocks), sum(blocks))) ; N20 = 0

    sb = 0 ; si = 0
    do B=1,4,2
      N = blocks(B) ; if(N.eq.0) cycle
      N2= blocks(B+1)
      T = N + N2

      U = Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T)
      V = Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)

      ! Note the extra minus sign when time-reversal is conserved
$TR   N20(si+1:si+T,si+1:si+T) =-matmul(transpose(U),V)-matmul(transpose(V),U)
$NTR  N20(si+1:si+T,si+1:si+T) = matmul(transpose(U),V)-matmul(transpose(V),U)
        
      si = si +     T
      sb = sb + 2 * T
    enddo
  end function calcN20

  pure function calcN11(Bogo, blocks) result(N11)
    !---------------------------------------------------------------------------
    ! Calculate the 11-component of the particle number operator
    !
    !   N^{20} = U^dagger U - V^{\dagger} V^*
    !
    ! Input:
    ! -------
    !   Bogo  : Bogoliubov transformation, only the right-most half is used
    !   blocks: symmetry block structure 
    ! 
    ! Output:
    ! -------
    !   N11   : 11-component of N
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: Bogo(:,:)
    real(KIND=dp), allocatable :: N11(:,:)
    real(KIND=dp), allocatable :: U(:,:), V(:,:)
    integer, intent(in)        :: blocks(4)

    integer :: B, si, N, N2, T, sb

    allocate(N11(sum(blocks), sum(blocks))) ; N11 = 0

    sb = 0 ; si = 0
    do B=1,4,2
      N = blocks(B) ; if(N.eq.0) cycle
      N2= blocks(B+1)
      T = N + N2

      U = Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T)
      V = Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)

      N11(si+1:si+T,si+1:si+T) = matmul(transpose(U),U)-matmul(transpose(V),V)

      si = si +     T
      sb = sb + 2 * T
    enddo
  end function calcN11

  pure function particle_number_bogo(bogo, occ, blocks) result(part)
    !-------------------------------------------------------------------------
    ! Calculate the particle number associated with a given Bogoliubov 
    ! transformation: sum the diagonal elements of rho.
    !
    ! Input:
    ! -------
    !   bogo  : bogoliubov transformation, only the right-most half is used
    !   occ   : quasiparticle occupation factors
    !   blocks: symmetry-block structure
    !
    ! Output:
    ! -------
    !   part  : total number of particles
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: bogo(:,:), occ(:)
    real(KIND=dp)             :: part
    real(KIND=dp), allocatable:: V(:,:), U(:,:)
    integer, intent(in)       :: blocks(4)

    integer :: B,  N, N2, T, si, sb, j
    
    part = 0.0d0
    si = 0 ; sb   = 0
    do B=1,4,2
      N = blocks(B) ; if(N.eq.0) cycle
      N2= blocks(B+1)
      T = N + N2

      U = Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T)
      V = Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)

      !    rho   = U   f U^\dagger + V^* (1 - f) V^T
      do j=1,T
        part = part + sum(        occ(si+j)  * U(:,j)**2) &
        &           + sum( (1.0d0-occ(si+j)) * V(:,j)**2)
      enddo
      si = si +     T
      sb = sb + 2 * T
   enddo

   ! Time*reversal factor 2
$TR  part = 2.0d0 * part
  end function particle_number_bogo
  
  function disp_bogo(bogo, occ, blocks) result(disp)
    !---------------------------------------------------------------------------
    ! Calculate the particle number dispersion
    !
    !   Delta N =  2 Tr( rho (1 - rho))
    !
    ! Input:
    ! -------
    !   bogo  : bogoliubov transformation, only the right-most half is used
    !   occ   : quasiparticle occupation factors
    !   blocks: symmetry-block structure
    !
    ! Output:
    ! --------
    !   disp  : particle number dispersion
    !
    !
    ! TODO: expression not correct yet when f!= 0
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: bogo(:,:), occ(:,:)
    real(KIND=dp)             :: disp
    real(KIND=dp), allocatable:: V(:,:), r(:,:), c(:,:)
    integer, intent(in)       :: blocks(4)

    integer :: B,  N, N2, T, sb, i
    
    disp = 0.0d0
    sb = 0
    do B=1,4,2
      N = blocks(B) ; if(N.eq.0) cycle
      N2= blocks(B+1)
      T = N + N2

      V = Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)
      
      r = matmul(transpose(V), V)
      c = matmul(r, r)
      
      do i=1, T
        disp = disp + r(i,i) - c(i,i) 
      enddo
      sb   = sb + 2 * T
   enddo
   disp  = 2 * disp
   ! Time*reversal factor 2
$TR  disp = 2 * disp 
  end function disp_bogo

  subroutine diagonalise_H_free(Bogo,prev,H11,N11,H20,N20,lambda,occ,blocks,Eqp)
    !---------------------------------------------------------------------------
    ! Diagonalise the HFB Hamiltonian in subspaces of quasiparticles with 
    ! equal occupation factors.
    !
    !  -> If all occupation factors are 0, like for even-even nuclei, then 
    !     we diagonalise H^11.
    !  -> If some occupations factors are 1, we diagonalise parts of H^11 and 
    !     parts of H^20 __separately__.
    !  -> If some occupations factors are 1/2 (Equal filling), we diagonalise
    !     parts of H^11 and H20 __together__.
    !
    ! If we are not performing equal filling, then this means diagonalising 
    ! the 11-component of H, i.e. H^{11}.
    !
    ! On output, this produces estimates for the quasiparticles. 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   Bogo   : Bogoliubov transformation
    !            ALL matrix elements are referenced.
    !   prev   : previous updates of the gradient step
    !   H11    : 11-component of the HFB Hamiltonian
    !   N11    : 11-component of the particle number operator
    !   H20    : 20-component of the HFB Hamiltonian
    !   N20    : 20-component of the particle number operator
    !   lambda : Fermi energy 
    !   occ    : quasiparticle occupation factors f
    !   block  : symmetry-block structure of the matrices
    !  
    ! Output:
    !   Bogo   : Bogoliubov transformation with the quasiparticles transformed
    !            unitarily within occupation blocks.
    !   Eqp    : Diagonal matrix elements of H^11 in the new basis
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: H11(:,:), N11(:,:), H20(:,:), N20(:,:)
    real(KIND=dp), intent(in)    :: occ(:), lambda
    real(KIND=dp), intent(inout) :: Bogo(:,:), prev(:,:,:)
    real(KIND=dp), intent(out)   :: Eqp(:)
    integer, intent(in)          :: blocks(4)
    
    real(KIND=dp), allocatable   :: block_transfo(:,:), A11(:,:), A20(:,:)
    real(KIND=dp), allocatable   :: temp(:,:), tempE(:), work(:), transfo(:,:)
    real(KIND=dp), allocatable   :: prev_full(:,:)
    integer :: k, sb, B, T, N, N2, si, s, NB, lwork, occ_blocks(3), ifail
    
    si = 0 ; sb = 0
    do B=1,4,2
      N = blocks(B) ; if(N.eq.0) cycle
      N2= blocks(B+1) 
      T = N + N2

      allocate(block_transfo(T,T), transfo(2*T,2*T))
      block_transfo = 0 ; transfo = 0.0d0
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Initial bookkeeping 
            
      ! Copy the relevant symmetry subblocks
      A11 = H11(si+1:si+T, si+1:si+T) - lambda * N11(si+1:si+T,si+1:si+T)
      A20 = H20(si+1:si+T, si+1:si+T) - lambda * N20(si+1:si+T,si+1:si+T)

      ! Separate things into blocks of equal occupation
      call occupation_blocks(T, occ(si+1:si+T), occ_blocks, block_transfo)

      ! Reorder the matrices according to occupation
      A11 = matmul(                     A11, block_transfo)
      A11 = matmul(transpose(block_transfo), A11)
      A20 = matmul(                     A20, block_transfo)
      A20 = matmul(transpose(block_transfo), A20)

      !-------------------------------------------------------------------------
      ! It is easier (for my insight) to work with the full prev matrix in this
      ! routine, instead of with the reduced one.
      allocate(prev_full(2*T,2*T))
      prev_full(  1:  T,  1:  T) = prev(si+1:si+T, si+1:si+T, 1)
      prev_full(T+1:2*T,T+1:2*T) = prev(si+1:si+T, si+1:si+T, 1)
      prev_full(  1:  T,T+1:2*T) = prev(si+1:si+T, si+1:si+T, 2)
      prev_full(T+1:2*T,  1:  T) =+prev(si+1:si+T, si+1:si+T, 2)

      !-------------------------------------------------------------------------
      ! Transform Bogoliubov matrix into this basis (as well as the prev array)
      ! through a simple column permutation. 
      call permute_columns(Bogo, prev_full, block_transfo, sb, T, .false.)
      !-------------------------------------------------------------------------
      ! Diagonalise the subblocks with f = 0 and f = 1
      s = 0   ! starting index of the subblock, to be incremented
      do k=1,2
        NB = occ_blocks(k)
        if( NB .eq. 0) cycle   ! Don't work if nothing in this occupation-block

        ! diagonalize A11 in this particular occupation block
        lwork = -1; allocate(work(1))
        call DSYEV( 'V', 'U', NB, A11(s+1:s+NB,s+1:s+NB), NB, &
        &                        Eqp(si+s+1:si+NB),work,lwork,ifail)
        lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
        call DSYEV( 'V', 'U', NB, A11(s+1:s+NB,s+1:s+NB), NB, &
        &                        Eqp(si+s+1:si+NB),work,lwork,ifail)
        deallocate(work)

        if(ifail.ne.0) then
          print *, 'Issue with DSYEV call in diagonalise_H_free.'
          print *, 'IFAIL = ', ifail
          stop
        endif
        ! Populate the transformation in this subblock        
        transfo(s  +1:s  +NB,s  +1:s  +NB) =  A11(s+1:s+NB,s+1:s+NB)
        transfo(s+T+1:s+T+NB,s+T+1:s+T+NB) =  A11(s+1:s+NB,s+1:s+NB)
        transfo(s  +1:s  +NB,s+T+1:s+T+NB) =  0.0d0
        transfo(s+T+1:s+T+NB,s  +1:s  +NB) =  0.0d0
        ! Note that the 20-part of this transformation is zero

        s = s + occ_blocks(k)
      enddo

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Diagonalise the f=0.5 subblock
      NB = occ_blocks(3)
      if(NB .ne. 0) then
          allocate(temp(2*NB, 2*NB)) ; temp  = 0.0d0
          allocate(tempE(2*NB))      ; tempE = 0.0d0

          ! Putting in the H^{11} part
          temp(   1:  NB,   1:  NB) =-A11(s+1:s+NB,s+1:s+NB)
          temp(NB+1:2*NB,NB+1:2*NB) =+A11(s+1:s+NB,s+1:s+NB)
          ! Putting in the H^{20} part
          temp(   1:  NB,NB+1:2*NB) =+A20(s+1:s+NB,s+1:s+NB)
          temp(NB+1:2*NB,   1:  NB) =+A20(s+1:s+NB,s+1:s+NB)
          ! TODO: There is likely some signs to be considered for T-reversal

          lwork = -1; allocate(work(1))
          call DSYEV( 'V', 'U', 2*NB, temp,2*NB, tempE,work,lwork,ifail)
          lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
          call DSYEV( 'V', 'U', 2*NB, temp,2*NB, tempE,work,lwork,ifail)
          Eqp(s+1:s+NB)     = tempE(NB+1:2*NB)          

          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! Reorder the negative energy qps: the diagonalisation routine 
          ! produces the ordering (-E_n, ..., -E_1, E_1, ...., E_n) but we
          ! need the textbook ordering (-E_1, ...., E_n, E_1, ..., E_n) here
          temp(:,1:NB) = temp(:,NB:1:-1) 

          !---------------------------------------------------------------------          
          ! Populate the transformation in this occupation block        
          transfo(s  +1:s+  NB,s  +1:s+  NB) = temp(   1:  NB,   1:  NB)
          transfo(s+T+1:s+T+NB,s  +1:s+  NB) = temp(NB+1:2*NB,   1:  NB)
          
          transfo(s  +1:s+  NB,s+T+1:s+T+NB) = temp(   1:  NB,NB+1:2*NB) 
          transfo(s+T+1:s+T+NB,s+T+1:s+T+NB) = temp(NB+1:2*NB,NB+1:2*NB)
          !---------------------------------------------------------------------          
          deallocate(temp, work, tempE)
      endif
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Transform the Bogoliubov transformation with the unitary transformations
      ! in all three occupation-blocks.
      Bogo(sb +1:sb+2*T, sb+1:sb+2*T) = matmul( &
      &           Bogo(sb+1:sb+2*T, sb+1:sb+2*T), transfo)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! If the Bogoliubov transformation transforms as 
      !
      !       Bogo' = Bogo U
      !
      ! Then the previous updates = the gradient, transforms as
      !
      !       U^dagger prev  U
      !
      prev_full = matmul(         prev_full,   transfo)
      prev_full = matmul(transpose(transfo), prev_full)
  
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Finally, we transform back to the original ordering of quasiparticles
      ! in the conventions of the rest of the code
      call permute_columns(Bogo, prev_full, block_transfo, sb, T, .true.)      

      ! and we fill the reduced prev matrix with the complete values we put.
      prev(si+1:si+T,si+1:si+T,1) = prev_full(  1:  T,  1:  T)
      prev(si+1:si+T,si+1:si+T,2) = prev_full(  1:  T,T+1:2*T)

      deallocate(block_transfo, transfo, prev_full)

      si = si +   T
      sb = sb + 2*T
    enddo

  end subroutine diagonalise_H_free
  
  subroutine permute_columns(Bogo, prev, transfo, sb, T, invert)
    !---------------------------------------------------------------------------
    ! Transform ONE SYMMETRY subblock of the Bogoliubov transformation as well 
    ! as the previous update (prev) into an "occupation-block" diagonal form.
    !
    !  ( U )  =>  (U^1_1, ...,U^1_n,U^0_1, ...,U^0_n,U^1/2_1, ...,U^1/2_n)
    !  ( V )      (V^1_1, ...,V^1_n,V^0_1, ...,V^0_n,V^1/2_1, ...,V^1/2_n)
    ! 
    ! where U/V^{1}   are the components of quasiparticles with occupation 1
    !       U/V^{0}   are the components of quasiparticles with occupation 0
    !       U/V^{1/2} are the components of quasiparticles with occupation 1/2
    !
    ! This subroutine merely applies the transformation generated by the 
    ! subroutine occupation_blocks. 
    !
    ! Additionally, this routine transforms FROM the convention of the rest
    ! of the program regarding quasiparticle energy ordering: 
    !
    !  (-E_N, -E_N-1, ..., -E_1, E_1, E_2, ...., E_N)
    !
    ! TO the text book ordering:
    !
    !  (-E_1, -E_2,   ..., -E_N, E_1, E_2, ...., E_N)
    !
    ! which is easier to use in the diagonalise_H_free subroutine.
    !
    ! The entirety of the transformation (occupation + energy ordering) can also
    ! be inverted by using the flag invert.
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input : 
    ! ---------
    !     Bogo   : Bogoliubov transformation. ALL elements are referenced.
    !     prev   : Previous Bogoliubov update, in complete form. 
    !              All matrix elements are referenced.
    !     transfo: column permutation, in the form of a matrix
    !              output of subroutine occupation_blocks
    !     sb     : starting index of the symmetry block in the Bogo matrix
    !     T      : size of this particular symmetry block
    !     invert : whether to do the inverse of transfo or not
    !
    ! Output:
    ! ---------
    !     Bogo   : reordered according to the description
    !     prev   : reordered according to the description
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(inout) :: Bogo(:,:)
    real(KIND=dp), intent(inout) :: prev(:,:)
    real(KIND=dp), intent(inout) :: transfo(:,:)
    logical, intent(in)          :: invert
    integer, intent(in)          :: sb, T

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This is the "normal" column permutation working on the right-hand side
    ! of the Bogoliubov and prev matrices
    if(invert) then
      transfo = transpose(transfo)
    endif

    Bogo(sb+1:sb+2*T, sb+T+1:sb+2*T) = matmul( &
    &                                 Bogo(sb+1:sb+2*T, sb+T+1:sb+2*T), transfo) 

    prev(  1:  T,T+1:2*T) = matmul(          prev(  1:  T,T+1:2*T), transfo) 
    prev(T+1:2*T,T+1:2*T) = matmul(          prev(T+1:2*T,T+1:2*T), transfo) 
    prev(  1:  T,T+1:2*T) = matmul(transpose(transfo),prev(  1:  T,T+1:2*T)) 
    prev(T+1:2*T,T+1:2*T) = matmul(transpose(transfo),prev(T+1:2*T,T+1:2*T)) 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! We can reorder prev with simple permutations, as we construct it as 
    !
    !    grad^11  grad^20 
    !    grad^20  grad^11
    !
    ! which means it already assumes the 'textbook' ordering of qps.
    prev(  1:  T,  1:  T) = matmul(              prev(  1:  T,   1: T), transfo) 
    prev(T+1:2*T,  1:  T) = matmul(              prev(T+1:2*T,  1:  T), transfo) 
    prev(  1:  T,  1:  T) = matmul(transpose(transfo),    prev(  1:  T,  1:  T)) 
    prev(T+1:2*T,  1:  T) = matmul(transpose(transfo),    prev(T+1:2*T,  1:  T)) 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! This is the same column permutation, but with the extra caveat that 
    ! we are switching conventions of the storage of the Bogoliubov transfo:
    !
    ! In most of the program, we have an ordering of qps that is 
    !   (-E_n, -E_n-1, E_n-2, ..., -E_1, E_1, .....,E_n)
    ! in order of increasing qp energy. 
    !
    ! To make the transformations below easier to read, we switch here the 
    ! ordering in the first half to: 
    !   (-E_1,-E_2, .........., E_1, E_2, ...., E_n)
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This is the reordering of transfo
    if(invert) then
      ! Go back to the original transfo before switching the rows
      transfo = transpose(transfo)
    endif
    transfo = transfo(T:1:-1,1:T)
    if(invert) then
      ! Inverting again
      transfo = transpose(transfo)
    endif
    Bogo(sb+1:sb+2*T,sb+1:sb+T) = matmul(Bogo(sb+1:sb+2*T,sb+1:sb+T), transfo) 

    ! Repair the reordering
    if(invert) transfo = transpose(transfo)
    transfo = transfo(T:1:-1,1:T)
  
  end subroutine permute_columns

  subroutine occupation_blocks(N, occ, sizes, transfo)
    !---------------------------------------------------------------------------
    ! Identify and reorder blocks in the Bogoliubov transformation with 
    ! identical occupation factors. We limit this at the moment to occupations 
    ! that are either 0, 1 or 0.5.
    ! 
    ! This routine counts the number of quasiparticles in each such 
    ! 'occupation-block' and constructs a simple transformation that can reorder
    ! the quasiparticles in the following ordering:
    !
    !          (f = 1 , f = 0, f = 0.5)
    !
    ! The current implementation uses a permutation matrix, i.e. a matrix 
    ! filled with zero's and a single 1 per column. This is a numerically
    ! inefficient way of implementing permutations, but it is visually easy as
    ! basis transformations are just calls to matmul.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    ! -------
    !         N : total size of occupations
    !       occ : occupation factors f
    ! 
    ! Output:
    ! -------
    !   sizes   : sizes of the subblocks with identical occupations
    !             Ordering: sizes(1) => F_ii = 0
    !                       sizes(2) => F_ii = 1
    !                       sizes(3) => F_ii = 0.5
    !   transfo : permutation matrix from 'current' to occupation-ordered
    !
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in)  :: occ(N)
    integer, intent(out)       :: sizes(3)
    real(KIND=dp), intent(out) :: transfo(N,N)
    integer                    :: indices(N,3),i, k, N, ind
    
    sizes = 0
    
    transfo = 0.0d0 ; indices = 0

    ! Count the number of qps in all three possible subblocks
    do i=1,N
      if(occ(i) .eq. 1.0d0) then
          ! Occupation factors 1.0
          sizes(1)             = sizes(1) + 1
          indices(sizes(1), 1) = i 
      elseif(occ(i) .eq. 0.0d0) then
          ! Occupation factors 0.0
          sizes(2)             = sizes(2) + 1
          indices(sizes(2), 2) = i 
      elseif(occ(i) .eq. 0.5d0) then
          ! Occupation factors 0.5
          sizes(3)             = sizes(3) + 1
          indices(sizes(3), 3) = i 
      else
          print *, 'Gradient solver only knows how to handle f=0, 0.5 or 1.'
          print ('(a20, f10.3)'), 'Offending entry = ', occ(i)
          print ('(99f10.3)'), occ(:)
          stop
      endif
    enddo    
    
    ! Construct the switching transformation
    ind = 1
    do k=1,3
      do i=1, sizes(k)
        transfo(indices(i,k),ind) = 1.0d0
        ind = ind + 1
      enddo  
    enddo
  
  end subroutine occupation_blocks

  subroutine symmetrize_bogo(Bogo, blocks)
    !---------------------------------------------------------------------------
    ! Symmetrize the Bogoliubov transformation: assuming the right-most half
    ! of (each symmetry block of) the array Bogo characterizes the correct
    ! Bogoliubov transformation, we construct the left-most half by symmetry.
    !
    ! Naively, we write
    !
    !   ( A  U )   =>  ( V^*  U )
    !   ( B  V )       ( U^*  V )
    !
    ! which is true up to a permutation of the left-most half. 
    ! The naive textbook ordering (shown above) follows
    !
    !  (-E_1, -E_2,   ..., -E_N, E_1, E_2, ...., E_N)
    !
    ! For historical reasons, we follow the ordering that is naturally obtained
    ! from diagonalisation routines:
    !
    !  (-E_N, -E_N-1, ..., -E_1, E_1, E_2, ...., E_N)
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    ! Input:
    ! ------
    !       Bogo  : Bogoliubov transformation. Only the right-most half of 
    !               every symmetry block is referenced. 
    !       blocks: symmetry structure of the matrix
    ! Output:
    ! ------
    !       Bogo  : Bogoliubov transformation, completely symmetrized. 
    !---------------------------------------------------------------------------
  
    real(KIND=dp), intent(inout) :: Bogo(:,:)
    integer, intent(in)          :: blocks(4)
    integer                      :: N, N2, T, sb, B, i
    
    sb = 0
    do B=1,4,2
      N = blocks(B)   ; if(N.eq.0) cycle
      N2= blocks(B+1)
      T = N + N2
  
      do i=1,T
        ! Note the extra minus sign when time-reversal is conserved.
$TR     Bogo(sb  +1:sb  +T, sb+T+1-i) =-Bogo(sb+T+1:sb+2*T, sb+T+i)   
$NTR    Bogo(sb  +1:sb  +T, sb+T+1-i) = Bogo(sb+T+1:sb+2*T, sb+T+i)
        Bogo(sb+T+1:sb+2*T, sb+T+1-i) = Bogo(sb  +1:sb+  T, sb+T+i)   
      enddo 
      sb = sb + 2*T
    enddo
    
  end subroutine symmetrize_bogo
  
end module HFB_gradient
