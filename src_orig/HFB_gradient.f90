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
 !
 !
 !------------------------------------------------------------------------------
 
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
  !      N  = ||  H^{20} - lambda N^{20} ||
  !
  ! where || is the Frobenius norm of the matrix. This is saved here, as the 
  ! routines below treat isospin independently.
  real(KIND=dp) :: HFBgradnorm(2) = 0
  !-----------------------------------------------------------------------------
  ! History of the gradients for the heavy-ball evolution. 
  ! Saved here, as the routines below can be used for either isospin
  real(KIND=dp), allocatable :: Z_updates(:,:)
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

  function buildgrad(H20, N20, lambda, Eqp, precon) result(grad)
    !---------------------------------------------------------------------------
    ! Construct the gradient of the energy. 
    !
    ! The "barre" gradient
    !     g  = H^{20} - lambda N^{20} 
    !
    ! which can be preconditioned by employing an approximation to the second
    ! derivative of the energy:
    !
    !     [Pg]_mn = g_mn/max(E^qp_m + E^qp_n , 2.0)
    ! 
    ! where this equation should only be used in the quasi-particle basis, i.e.
    ! the basis that diagonalizes H^{11}. 
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable :: grad(:,:)
    real(KIND=dp), intent(in)  :: H20(:,:), N20(:,:), Eqp(:), lambda
    logical                    :: precon 

    grad = H20 - lambda * N20
    if(precon) then
      grad = precon_grad(grad, Eqp)
    endif
    
  end function buildgrad

  subroutine gradient_step(h,gaps,targetN, Bogo,Eqp,lambda, alpha,mu,prev,     &
  &                        precon, gradnorm, blocks, maxiter, ifail)
    !---------------------------------------------------------------------------
    !
    ! Perform one (or more) heavy-ball evolution steps, starting from an initial
    ! Bogoliubov transformation with given single-particle hamiltonian h and 
    ! pairing gaps Delta. 
    !
    ! Input:
    !       h : single-particle hamiltonian 
    !     gaps: pairing gaps
    !  targetN: targetted number of particles
    !    Bogo : Bogoliubov transformation on input
    !    Eqp  : an estimate for the quasi-particle energies for the 
    !           preconditioning of the evolution (if requested)
    !   lambda: Fermi energy
    !   alpha : step-size for the heavy-ball evolution
    !      mu : momentum for the heavy-ball evolution
    !   precon: whether or not to precondition the update
    !   blocks: size of the quantum number blocks for this evolution 
    !  maxiter: number of steps to do (recommended: 1)
    !
    ! Output:
    !    Bogo : new Bogoliubov transformation
    !    Eqp  : estimated qp energies based on the diagonalisation of  H^{11}
    !   lambda: new value for the Fermi energy
    ! gradnorm: Frobenius norm of the gradient that was used as step
    !  ifail  : if 0, succes. If 1, something went wrong.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Note that we assume a specific formatting of the Bogoliubov transformation
    ! in this routine: 
    !
    !              (  V^*   U   )
    !       B  =   (            )
    !              (  U^*   V   )
    !
    ! where the right-hand half of the columns are the actual U and V 
    ! columns that have been 'picked' to construct the Bogoliubov reference 
    ! state. In the context of the HFB module, this means that this routine 
    ! only works correctly if the configuration matrix is trivial, that is to
    ! to say
    !
    !      C   = ( 0  0 )
    !            ( 0  1 )
    !
    ! This module ONLY evolves the right-most half of the Bogoliuv matrix!
    !---------------------------------------------------------------------------  
    integer, intent(in)          :: blocks(4), maxiter
    integer, intent(out)         :: ifail
    real(KIND=dp), intent(in)    :: targetN
    real(KIND=dp), intent(in)    :: h(:,:), gaps(:,:), alpha, mu
    logical, intent(in)          :: precon
    real(KIND=dp), intent(inout) :: Bogo(:,:), lambda
    real(KIND=dp), intent(inout) :: Eqp(:),  prev(:,:), gradnorm
    real(KIND=dp), allocatable   :: H20(:,:), N20(:,:), grad(:,:), H11(:,:)

    real(KIND=dp)                :: particles,  normN
    integer                      :: iter
    logical                      :: converged = .false.  
  
    converged = .false. 
    do iter=1,maxiter
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! We calculate the relevant matrices to build the gradient 
        H20 = calcH20(Bogo,h,gaps,blocks)
        N20 = calcN20(Bogo, blocks)

        ! First, we check if the pairing has collapsed: the Frobenius norm of 
        ! the N20 gives us an indication.
        normN = sqrt(sum(N20**2))        

        if(normN.gt.1d-4) then
          ! If the pairing has not collapsed, we try different values of the 
          ! Fermi energy to arrive at a correct particle number AFTER the 
          ! heavy-vall step. 
          call find_fermi_brent(Bogo, H20, N20, prev, Eqp, lambda,             &
          &                particles, targetN, alpha, mu, precon, blocks, ifail)
        endif
        ! Building the gradient update 
        ! (with the old Fermi energy if the pairing collapsed)
        grad = buildgrad(H20, N20, lambda, Eqp, precon)
        ! ...  and executing it
        bogo = GradUpdate(grad, prev, bogo, alpha, mu, blocks)
        ! ...  and saving it for the next iteration 
        prev = -alpha * grad + mu * prev
        ! and for good measure, we recalculate the (deviation of) the 
        ! particle number
        particles = particle_number_bogo(bogo, blocks) - targetN
        ! as well as the norm of the gradient 
        gradnorm  = sqrt(sum(grad**2))

        !expectedDE = expectedDE - alpha * sum((H20 - lambda * N20)**2)
    
        ! Check for convergence if this is process is repeated multiple times
        if(gradnorm .lt. 1d-6) converged = .true.

        if(converged) then
          exit
        endif
    enddo
    ! Perform an additional transformation of the Bogoliubov transformation to
    ! diagonalise H^11 and obtain another estimate for the QP energies
    H11 = calcH11(Bogo, h, gaps, lambda, blocks)
    call diagonalise_H11(bogo, prev, H11, blocks, Eqp)

  end subroutine gradient_step
  
  function precon_grad(grad, Eqp)  result(Pgrad)
    !---------------------------------------------------------------------------
    !
    !    [Pg]_mn = g_mn/max(E^qp_m + E^qp_n , 2.0)
    !
    !
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: grad(:,:), Eqp(:)
    real(KIND=dp)              :: fac
    real(KIND=dp), allocatable :: Pgrad(:,:)
    integer   :: i,j

    Pgrad = grad
    do i=1, size(grad,1)
        do j=1, size(grad,1)
          fac = max(Eqp(i)+ Eqp(j), 2.0d0)
          Pgrad(i,j) = Pgrad(i,j)/fac
        enddo
    enddo

  end function precon_grad

  subroutine find_fermi_brent(Bogo, H20, N20,  prev, Eqp, lambda,              &
  &                        particles, targetN, alpha, mu, precon, blocks, ifail)
    !---------------------------------------------------------------------------
    ! Find a Fermi energy such that the particle number is (on average) correct
    ! AFTER the heavy-ball evolution.
    ! 
    ! (NOTE: subroutine and documentation were closely copied from the one in 
    ! the HFB_direct module)
    !---------------------------------------------------------------------------
    ! This particular subroutine employs Brents method to fix the Fermi 
    ! energy. See
    ! 
    ! https://en.wikipedia.org/wiki/Brent%27s_method
    ! 
    ! which combines bisection, secant method and inverse quadratic 
    ! interpolation.The original source is probably
    ! R. P. Brent (1973), "Chapter 4: An Algorithm with Guaranteed Convergence
    ! for Finding a Zero of a Function", Algorithms for Minimization without
    ! Derivatives, Englewood Cliffs, NJ: Prentice-Hall,  
    !
    !-------------------------------------------------------------------------
    ! This routine is very heavily inspired/copy-pasted by the routines 
    ! implemented in MOCCa by M. Bender. 
    !-------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: H20(:,:), N20(:,:), Eqp(:)
    real(KIND=dp), intent(in)    :: targetN 
    integer, intent(in)          :: blocks(4)
    real(KIND=dp), intent(in)    :: alpha, mu, prev(:,:)
    logical, intent(in)          :: precon

    real(KIND=dp), intent(out)   :: particles
    real(KIND=dp), intent(inout) :: lambda, bogo(:,:)
    real(KIND=dp), allocatable   :: gradA(:,:), gradB(:,:), nbA(:,:), nbB(:,:)

    real(KIND=dp)                :: InitialBracket(2), FA, FB, N
    integer                      :: idir = 0 , idirsig = 1, FailCount, ifail
    logical                      :: Success

    !---------------------------------------------------------------------------
    ! STEP 1: set up an initial bracket
    !---------------------------------------------------------------------------
    gradA = buildgrad(H20, N20, lambda, Eqp, precon)
    nbA   = GradUpdate(gradA, prev, bogo, alpha, mu, blocks)
    N     = particle_number_bogo(nbA, blocks) - targetN
    
    ! Check if this guess for lambda is good enough
    if(abs(N).lt.pairing_prec .or. alpha .eq. 0.0d0) then
      !bogo = nbA
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

        gradA = buildgrad(H20, N20,InitialBracket(1), Eqp, precon)
        gradB = buildgrad(H20, N20,InitialBracket(2), Eqp, precon)

        nbA = GradUpdate(gradA, prev, bogo, alpha, mu, blocks)
        nbB = GradUpdate(gradB, prev, bogo, alpha, mu, blocks)

        FA = particle_number_bogo(nbA, blocks) - targetN
        FB = particle_number_bogo(nbB, blocks) - targetN

        ! diagnostic printing for convergence analysis (usually commented out)
        !print '(" Bracketing ",i4,1l2,(2(f13.8,es16.7)))',        &
        !      & FailCount,Success,InitialBracket(1),FA, InitialBracket(2),FB          
        ! code failure (Fermi energy has changed by 30 MeV)
        if (Failcount .gt. 76) then
          print '(/," A = ", f13.8, " FA = ",1es12.4,              &
               &    " B = ", f13.8, " FB = ",1es12.4)',            &
               &     InitialBracket(1),FA,InitialBracket(2),FB 
          ifail = 1
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
      &                    H20, N20, prev, Eqp, lambda,                        &
      &                    particles, targetN, alpha, mu, precon, blocks,200)

  end subroutine find_fermi_brent

  subroutine Brent_bisection(X1,X2,FX1, FX2, Bogo, H20, N20, prev, Eqp, lambda,& 
  &                        particles, targetN, alpha, mu, precon, blocks, depth)
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
    real(KIND=dp), intent(in)    :: H20(:,:), N20(:,:), Eqp(:)
    real(KIND=dp), intent(in)    :: alpha, mu, targetN
    real(KIND=dp), intent(in)    :: X1 , X2, FX1 , FX2, prev(:,:)
    integer, intent(in)          :: blocks(4), depth
    logical, intent(in)          :: precon 

    real(KIND=dp), intent(out)   :: particles
    real(KIND=dp), intent(inout) :: lambda, bogo(:,:)
    
    real(KIND=dp), allocatable   :: grad(:,:), newbogo(:,:)
    real(KIND=dp)                :: A , B, C , FA, FB , FC
    real(KIND=dp)                :: D , E, S , P  , Q , R 
    real(KIND=dp)                :: Num , Tol , XM 
    real(KIND=dp)                :: eps = 1.d-9
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
      grad = buildgrad(H20, N20,B, Eqp, precon)
      newbogo = GradUpdate(grad, prev, bogo, alpha, mu, blocks)
      Num  = particle_number_bogo(newbogo, blocks)
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
        print '(/," Warning: BrentBisection did not converge after ",i4," iterations")', & 
        &      FailCount
      endif
    enddo
    ! Output
    Lambda    = B ; particles = FB 
    !bogo      = newbogo
  end subroutine Brent_bisection

  function GradUpdate(grad, prev, bogo, alpha, mu, blocks) result(NB)
    !---------------------------------------------------------------------------
    ! Update the bogoliubov transformation with a heavy-ball step:
    !
    !     U'  =  U - alpha * V * Z + mu * dU 
    !     V'  =  V - alpha * U * Z + mu * dV
    ! 
    ! followed by an orthonormalisation via a Gramm-Schmidt algorithm.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    !     alpha   : step size
    !     mu      : momentum of the evolution 
    !     grad    : direction of the update, Z
    !     prev    : previous update, (dU,dV)^T
    !     bogo    : (on input) current value of the Bogoliubov transformation
    !      
    !     NB      : (on output) updated Bogoliubov transformation
    ! 
    !     Note that this routine only deals with part of the Bogoliubov 
    !     transformation, as does the rest of the module: 
    ! 
    !    B_complete = ( V^*  U )    B_here = ( U )
    !                 ( U^*  V )             ( V )
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: alpha, mu 
    real(KIND=dp), intent(in)    :: grad(:,:), prev(:,:)
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
      
      ! U'  = U - alpha V Z + mu dU 
      NB(sb  +1:sb+T ,sb+T+1:sb+2*T) = NB(sb  +1:sb+T ,sb+T+1:sb+2*T)          &
      &                       - alpha *matmul(V,grad(si+1:si+T,si+1:si+T))     &
      &                       + mu    *matmul(V,prev(si+1:si+T,si+1:si+T))

      ! V = V - alpha U Z + mu dV
      ! Note the additional minus sign for time-reversal conserved calculations
      NB(sb+T+1:sb+2*T,sb+T+1:sb+2*T) = NB(sb+T+1:sb+2*T,sb+T+1:sb+2*T)        &
$NTR      &               -alpha *matmul(U,grad(si+1:si+T,si+1:si+T))          &
$TR       &               +alpha *matmul(U,grad(si+1:si+T,si+1:si+T))          &
          &               + mu * matmul(U,prev(si+1:si+T,si+1:si+T))
      si = si +   T
      sb = sb + 2*T
    enddo
  
    !Don't forget to orthonormalise 
    call ortho_bogo(NB, blocks)
    
  end function GradUpdate

  subroutine ortho_bogo(bogo, blocks)
    !---------------------------------------------------------------------------
    ! Gramm-Schmidt routine to orthogonalise the Bogoliubov transformation. 
    ! As with the rest of this module, only operates on the right-most half
    ! of the Bogoliubov transformation. 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(inout) :: bogo(:,:)
    integer, intent(in)          :: blocks(4)
    integer                      :: i,j,k, N, N2, sb, T, B
    real(KIND=dp)                :: overlap

    sb = 0
    do B=1,4,2
      N = blocks(B) ; if(N.eq.0) cycle
      N2= blocks(B+1)

      T = N+N2
      do j=1,T
        ! Normalise vector ( U_j )
        !                  ( V_j )
        overlap = 0.0_dp
        do i=1,2*T
            Overlap = Overlap + Bogo(sb+i,sb+T+j)**2
        enddo
        !print *, j,overlap
        Overlap    = 1.0_dp/sqrt(overlap)
        bogo(sb+1:sb+2*T,sb+T+j) = bogo(sb+1:sb+2*T,sb+T+j)*overlap
        ! Orthogonalise all the rest against vector j
        do i=j+1,T
          Overlap = 0.0_dp
          do k=1,2*T
              Overlap = Overlap + bogo(sb+k,sb+T+j) * bogo(sb+k,sb+T+i)
          enddo
          !print *, i, overlap
          Bogo(:,sb+T+i) = Bogo(:,sb+T+i) - Overlap * Bogo(:,sb+T+j)
        enddo
        !print *
      enddo
      sb = sb+ 2*T
    enddo
  end subroutine ortho_bogo

  function calcH20(Bogo,h,gaps, blocks) result(H20)
    !---------------------------------------------------------------------------
    ! Calculate the 2-quasi-particle-excitation component of H:
    ! 
    ! H20 =   U^{dagger} h   V^* - V^{\dagger} \Delta^* V^*
    !       - V^{dagger} h^t U^* + U^{\dagger} \Delta   U^*
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)               :: H(:,:), bogo(:,:), gaps(:,:)
    integer, intent(in)                     :: blocks(4)


    real(KIND=dp), allocatable :: H20(:,:), U(:,:), V(:,:)
    real(KIND=dp), allocatable :: hV(:,:), hU(:,:), dV(:,:), dU(:,:)
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

      !  h V^* and h^t U^*
      hV = matmul(   h(si+1:si+T,si+1:si+T), V)
      hU = matmul(   h(si+1:si+T,si+1:si+T), U)

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

  function calcH11(Bogo,h,gaps,lambda,blocks) result(H11)
    !---------------------------------------------------------------------------
    ! Calculate the 11 component of H 
    !  
    !   H^11 = U^T h U + U^T Delta V - V^T Delta U - V^t h V
    ! 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: h(:,:), bogo(:,:), gaps(:,:)
    real(KIND=dp), intent(in)  :: lambda
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

      !  h V^* and h^t U^*
      hV = matmul(   h(si+1:si+T,si+1:si+T), V) - lambda * V
      hU = matmul(   h(si+1:si+T,si+1:si+T), U) - lambda * U

      ! d^* V^* and dU^*
      dV = matmul(gaps(si+1:si+T,si+1:si+T), V)
      dU = matmul(gaps(si+1:si+T,si+1:si+T), U)

      ! We reuse the defined symbols to save a matrix multiplication here
      U = transpose(U) ; V = transpose(V)

      ! We can save some effort here in the future, H20 is antisymmetric     
      H11(si+1:si+T, si+1:si+T)  = matmul(U, hU) + matmul(U, dV) &
                              &  - matmul(V, dU) - matmul(V, hV)

      si = si +  T
      sb = sb +2*T
    enddo
  end function calcH11

  function calcN20(Bogo, blocks) result(N20)
    !---------------------------------------------------------------------------
    ! Calculate the 20-component of the particle number operator
    !
    !   N^{20} = U^dagger V - V^{\dagger}, U
    !
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

  function particle_number_bogo(bogo, blocks) result(part)
      !-------------------------------------------------------------------------
      ! Calculate the particle number associated with a given Bogoliubov 
      ! transformation
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in) :: bogo(:,:)
      real(KIND=dp)             :: part
      real(KIND=dp), allocatable:: V(:,:)
      integer, intent(in)       :: blocks(4)

      integer :: B,  N, N2, T, sb
      
      part = 0.0d0
      sb = 0
      do B=1,4,2
        N = blocks(B) ; if(N.eq.0) cycle
        N2= blocks(B+1)
        T = N + N2

        V = Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)

        part = part + sum(V**2)
        sb = sb + 2 * T
     enddo

     ! Time*reversal factor 2
$TR  part = 2* part
  end function particle_number_bogo

  subroutine diagonalise_H11(Bogo, prev, H11, blocks, Eqp)
    !---------------------------------------------------------------------------
    ! Diagonalise a (precalculated) H^11 to obtain quasiparticle energies
    ! but (more importantly)                   
    !
    ! * Transform the Bogoliubov transformation to this basis
    ! * as well as all the previous updates
    ! 
    ! Note that this transformation does not alter the many-body state 
    ! (and hence no observables), but other parts of the code DO depend on the 
    ! HFB Hamiltonian to be actually diagonalised (such as the calculation of
    ! the rotational correction).
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: H11(:,:)
    real(KIND=dp), intent(inout) :: Bogo(:,:), prev(:,:)
    real(KIND=dp), intent(out)   :: Eqp(:)
    real(KIND=dp), allocatable   :: A(:,:), work(:)
    integer, intent(in)          :: blocks(4)
    integer                      :: si, sb, B, N, N2, lwork, ifail, T
      
    si = 0 ; sb = 0
    do B=1,4,2
      N = blocks(B) ; if(N.eq.0) cycle
      N2= blocks(B+1) 
      T = N + N2

      A = H11(si+1:si+T, si+1:si+T)

      ! Diagonalize the first symmetry subblock
      lwork = -1; allocate(work(1))
      call DSYEV( 'V', 'U', N, A(1:N,1:N), N, Eqp(si+1:si+N),work,lwork,ifail)
      lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
      call DSYEV( 'V', 'U', N, A(1:N,1:N), N, Eqp(si+1:si+N),work,lwork,ifail)
      deallocate(work)
      
      ! Diagonalize the second symmetry subblock
      lwork = -1; allocate(work(1))
      call DSYEV( 'V', 'U', N2, A(N+1:T,N+1:T), N2, &
      &                                       Eqp(si+N+1:si+T),work,lwork,ifail)
      lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
      call DSYEV( 'V', 'U', N2, A(N+1:T,N+1:T), N2, &
      &                                       Eqp(si+N+1:si+T),work,lwork,ifail)
      deallocate(work)
        
      !-------------------------------------------------------------------------
      ! Transform the U and V matrices
      ! Commented out for the moment, but this allows us to transform into the 
      ! correct basis for preconditioning
      bogo(sb  +1:sb+  T, sb+T+1:sb+2*T) = &
                                   matmul(Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T),A)
      bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T) = & 
                                   matmul(Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T),A)
      ! But also transform the past updates, such that they stay consistent.                              
      prev(si  +1:si+  T, si+1:si+T) = &
                             matmul(transpose(A),prev(si  +1:si+  T, si+1:si+T))
      prev(si  +1:si+  T, si+1:si+T) = &
                             matmul(prev(si  +1:si+  T, si+1:si+T), A)
      !-------------------------------------------------------------------------

      si = si +   T
      sb = sb + 2*T
    enddo
  end subroutine diagonalise_H11
end module HFB_gradient
