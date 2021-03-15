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

  use geninfo
  use wavefunctions

  implicit none

contains 

  subroutine gradient_step(h,gaps,blocks,targetN,Bogo,Eqp,config,lambda,maxiter)
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------  
    integer, intent(in)          :: blocks(4), maxiter
    real(KIND=dp), intent(inout) :: Bogo(:,:), lambda, targetN
    real(KIND=dp), intent(inout) :: Eqp(:)
    real(KIND=dp), intent(in)    :: h(:,:), gaps(:,:), config(:)
    real(KIND=dp), allocatable   :: H20(:,:), N20(:,:), grad(:,:)
    real(KIND=dp), allocatable   :: H11(:,:)
    real(KIND=dp)                :: alpha, particles, gradnorm
    integer                      :: B, sb, N, N2, iter, fermiiter, T

    logical                      :: converged = .false.  
  
    converged = .false. 

    !---------------------------------------------------------------------------
    do iter=1,maxiter
        ! Calculate the gradient 
        H20 = calcH20(Bogo,h,gaps,blocks)
        N20 = calcN20(Bogo, blocks)

        call find_fermi_secant(Bogo, H20, N20,  Eqp, lambda, gradnorm,         &
        &                            particles, targetN, config, blocks)

        H11 = calcH11(Bogo, h, gaps, lambda, blocks)
        call diagonalise_H11(bogo, H11, blocks, Eqp)

        ! Calculate the norm of the gradient and check for convergence
        if(gradnorm .lt. 1d-12) converged = .true.

        if(converged) then
          exit
        endif
    enddo
    print ('(a4,2e12.3,99f10.3)'), 'GRAD',  gradnorm, (Particles - targetN), alpha, maxval(Eqp), minval(abs(Eqp))

    !---------------------------------------------------------------------------
    ! Finally, we transform the rest of the Bogoliubov transformation
    sb = 0
    do B=1,4,2
      N = blocks(B)   ; if(N.eq.0) cycle
      N2= blocks(B+1)
      T = N + N2
      ! Populate the columns of the Bogoliubov transformation that have not 
      ! been evolved. Note the extra minus sign when time-reversal is 
      ! conserved.
$TR   Bogo(sb  +1:sb  +T, sb+1:sb+T) =-Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)   
$NTR  Bogo(sb  +1:sb  +T, sb+1:sb+T) = Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)
      Bogo(sb+T+1:sb+2*T, sb+1:sb+T) = Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T)   

      sb = sb + 2*T
    enddo

  end subroutine gradient_step
  
  function precon_grad(grad, Eqp)  result(Pgrad)
      
    real(KIND=dp), intent(in)  :: grad(:,:), Eqp(:)
    real(KIND=dp)              :: fac
    real(KIND=dp), allocatable :: Pgrad(:,:)
    integer   :: i,j

    Pgrad = grad
    do i=1, size(grad,1)
        do j=1, size(grad,1)
          fac = max(Eqp(i)+Eqp(j), 2.0d0)

          Pgrad(i,j) = Pgrad(i,j)/fac
        enddo
    enddo

  end function precon_grad

  subroutine find_fermi_secant(Bogo, H20, N20,  Eqp, lambda, gradnorm,         &
  &                            particles, targetN, config, blocks)
    !---------------------------------------------------------------------------
    !
    !
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: H20(:,:), N20(:,:), Eqp(:)
    real(KIND=dp), intent(in)    :: targetN, config(:)
    integer, intent(in)          :: blocks(8)

    real(KIND=dp)                :: alpha
    real(KIND=dp)                :: df, dn(2)
    real(KIND=dp), intent(out)   :: particles, gradnorm
    real(KIND=dp), intent(inout) :: lambda, bogo(:,:)
    real(KIND=dp), allocatable   :: grad(:,:), newbogo(:,:)

    integer :: iter
  
    if(.not. all(Eqp.eq.0.0d0)) then
      alpha = 0.005 !10.0/(minval(abs(Eqp)) + maxval(Eqp))
    else
      alpha = 0.005
    endif

    particles = particle_number_bogo(bogo, config, blocks)
    df = 0.0d0; dn = 0.0d0
    do iter=1,100
      grad = H20 - lambda * N20
      !grad = precon_grad(grad, Eqp)
      newbogo = GradUpdate(grad, bogo, alpha, blocks)
      particles = particle_number_bogo(newbogo, config, blocks)

      dn(2) = dn(1)
      dn(1) = particles - targetN

      if(iter.eq.1) then
         ! We try lambda + 0.1 for the first iteration
         lambda = lambda + 0.1
         df     =          0.1
      else
         df     = - dn(1) * df/(dn(1) - dn(2))
  
         if(abs(df).gt.1.0) df = 0.1 * df/abs(df)
        lambda = lambda + df
      endif

      if(abs(Particles-targetN) .lt. 1d-14) exit 
      !lambda = lambda - (Particles - targetN)   
    enddo
    Bogo = newbogo
    gradnorm = sqrt(sum(grad**2))

  end subroutine find_fermi_secant

  function GradUpdate(grad, bogo, alpha, blocks) result(newbogo)
    !---------------------------------------------------------------------------
    ! Update the bogoliubov transformation with a step in the direction of
    ! "grad" with a step size alpha.
    !
    !     alpha   : step size
    !     grad    : direction of the update
    !     bogo    : (on input) current value of the Bogoliubov transformation
    !      
    !     new_bogo: (on output) updated Bogoliubov transformation
    ! 
    !     Note that this routine only deals with part of the Bogoliubov 
    !     transformation:
    ! 
    !    B_complete = ( V^*  U )    B_here = ( U )
    !                 ( U^*  V )             ( V )
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: alpha
    real(KIND=dp), intent(in)    :: grad(:,:)
    integer, intent(in)          :: blocks(4)
    real(KIND=dp), intent(in   ) :: bogo(:,:)
    real(KIND=dp), allocatable   :: tV(:,:), tU(:,:), aux(:,:), newbogo(:,:)
    integer                      :: B, sb, N, N2, T,  i, info, si

    newbogo = bogo
    sb = 0  ; si = 0
    do B=1,4,2
      N  = blocks(B)   ; if(N.eq.0) cycle
      N2 = blocks(B+1)
      T = N + N2

      tU = bogo(sb+  1:sb+  T,sb+T+1:sb+2*T)
      tV = bogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)

      ! U = U - alpha V gradient
      newbogo(sb  +1:sb+T ,sb+T+1:sb+2*T)=newbogo(sb  +1:sb+T ,sb+T+1:sb+2*T)  &
      &                         -alpha *matmul(tV,grad(si+1:si+T,si+1:si+T)) 

      ! V = V - alpha U gradient
      newbogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)=newbogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)&
$NTR      &                         -alpha *matmul(tU,grad(si+1:si+T,si+1:si+T)) 
$TR       &                         +alpha *matmul(tU,grad(si+1:si+T,si+1:si+T)) 

      si = si +   T
      sb = sb + 2*T
    enddo
  
    !Don't forget to orthonormalise (which could also be achieved through a 
    ! Cholesky decomposition)
    call ortho_bogo(newbogo, blocks)
    
  end function GradUpdate

  subroutine ortho_bogo(bogo, blocks)
    !---------------------------------------------------------------------------
    ! Gramm-Schmidt routine to orthogonalise a Bogoliubov matrix
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

  function calcH20(Bogo,h,gaps,blocks) result(H20)
    !---------------------------------------------------------------------------
    !
    ! H20 =   U^{dagger} h   V^* - V^{\dagger} \Delta^* V^*
    !       - V^{dagger} h^t U^* + U^{\dagger} \Delta   U^*
    !
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)               :: H(:,:), bogo(:,:), gaps(:,:)
    integer, intent(in)                     :: blocks(4)


    real(KIND=dp), allocatable :: H20(:,:), U(:,:), V(:,:)
    real(KIND=dp), allocatable :: hV(:,:), hU(:,:), dV(:,:), dU(:,:)
    integer                    :: B, N, N2, si, sb, i, j, T

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
$TR      H20(si+1:si+T, si+1:si+T)  =-matmul(U, hV) + matmul(U, dU) &
$TR                              &  - matmul(V, hU) - matmul(V, dV)

      si = si +  T
      sb = sb +2*T
    enddo
  end function calcH20

  function calcH11(Bogo,h,gaps,lambda,blocks) result(H11)
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: h(:,:), bogo(:,:), gaps(:,:)
    real(KIND=dp), intent(in)  :: lambda
    real(KIND=dp), allocatable :: H11(:,:), U(:,:), V(:,:)
    real(KIND=dp), allocatable :: hV(:,:), hU(:,:), dV(:,:), dU(:,:)
    integer, intent(in)        :: blocks(4)
    integer                    :: B, N, N2, si, sb, i, j, T

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
                              &  + matmul(V, dU) - matmul(V, hV)

      !U = transpose(U) ; V = transpose(V)

      !print *, 'H'
      !do i=1, T
      !  print ('(99f10.3)'), h(si+i, si+1:si+T)
      !enddo
      !print *!
      !print *, 'U'
      !do i=1, T
      !  print ('(99f10.3)'), U(i, 1:T)
      !enddo
      !print *
      !print *, 'V'
      !do i=1, T
      !  print ('(99f10.3)'), V(i, 1:T)
      !enddo
      !print *

      !print *, 'HV'
      !do i=1, T
      !  print ('(99f10.3)'), hV(i, 1:T)
      !enddo
      !print *
      !print *, 'HU'
      !do i=1, T
      !  print ('(99f10.3)'), hU(i, 1:T)
      !enddo
      !print *!!!

      !print *, 'H11'
      !do i=1, T
      !  print ('(99f10.3)'), H11(si+i, si+1:si+T)
      !enddo
      !print *

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

    integer :: B, si, N, N2, T, sb,i

    allocate(N20(sum(blocks), sum(blocks))) ; N20 = 0

    sb = 0 ; si = 0
    do B=1,4,2
      N = blocks(B) ; if(N.eq.0) cycle
      N2= blocks(B+1)
      T = N + N2

      U = Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T)
      V = Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)

$TR   N20(si+1:si+T,si+1:si+T) =-matmul(transpose(U),V)-matmul(transpose(V),U)
$NTR  N20(si+1:si+T,si+1:si+T) = matmul(transpose(U),V)-matmul(transpose(V),U)
        
      si = si +     T
      sb = sb + 2 * T
    enddo
  end function calcN20

  function particle_number_bogo(bogo, config, blocks) result(part)
      !-------------------------------------------------------------------------
      !
      !
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in) :: bogo(:,:), config(:)
      real(KIND=dp)             :: part
      real(KIND=dp), allocatable:: V(:,:), U(:,:)
      integer, intent(in)       :: blocks(4)

      integer :: B, si, N, N2, T, sb,i,j
      part = 0.0d0
      sb = 0
      do B=1,4,2
        N = blocks(B) ; if(N.eq.0) cycle
        N2= blocks(B+1)
        T = N + N2

        U = Bogo(sb+  1:sb+  T, sb+T+1:sb+2*T)
        V = Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)
        !print *, 'CONFIG', config(sb+1:sb+2*T)        
        do i=1,T  
          do j=1,T
            part = part +     V(i,j) **2  !&config(sb+T+j) *
!            &           + (1-config(sb+T+j))* U(i,j) **2
          enddo
        enddo        
        sb = sb + 2 * T
        !print ('(90f10.3)'), config(sb+T+1:sb+2*T)
     enddo

     ! Time*reversal factor 2
$TR  part = 2* part
  end function particle_number_bogo

  subroutine diagonalise_H11(Bogo, H11, blocks, Eqp)
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: H11(:,:)
    real(KIND=dp), intent(inout) :: Bogo(:,:)
    real(KIND=dp), intent(out)   :: Eqp(:)
    real(KIND=dp), allocatable   :: A(:,:), work(:)
    integer, intent(in)          :: blocks(4)
    integer                      :: si, sb, i,j, B, N, N2, lwork, ifail, T
      
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
      call DSYEV( 'V', 'U', N2, A(N+1:T,N+1:T), N2, Eqp(si+N+1:si+T),work,lwork,ifail)
      lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
      call DSYEV( 'V', 'U', N2, A(N+1:T,N+1:T), N2, Eqp(si+N+1:si+T),work,lwork,ifail)
      deallocate(work)
      
      !lwork = -1; allocate(work(1))
      !call DSYEV( 'V', 'U', T, A, T, Eqp(si+1:si+T),work,lwork,ifail)
      !lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
      !call DSYEV( 'V', 'U', T, A, T, Eqp(si+1:si+T),work,lwork,ifail)
      !deallocate(work)

      !print ('(a3,99f10.3)') , 'EQP', EQP(si+1:si+T)

      ! Transform the U and V matrices
      Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T) = matmul(Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T),A)
      Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T) = matmul(Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T),A)

      si = si +   T
      sb = sb + 2*T
    enddo
  end subroutine diagonalise_H11
end module HFB_gradient
