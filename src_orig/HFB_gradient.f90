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

  subroutine gradient_step(h,gaps,blocks,targetN,Bogo,Eqp,lambda,maxiter)
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------  
    integer, intent(in)          :: blocks(4), maxiter
    real(KIND=dp), intent(inout) :: Bogo(:,:), lambda, targetN
    real(KIND=dp), intent(inout) :: Eqp(:)
    real(KIND=dp), intent(in)    :: h(:,:), gaps(:,:)
    real(KIND=dp), allocatable   :: H20(:,:), N20(:,:), newbogo(:,:), grad(:,:)
    real(KIND=dp), allocatable   :: H11(:,:)
    real(KIND=dp)                :: alpha, particles, gradnorm
    integer                      :: B, si, sb, N, N2, iter, fermiiter,i

    logical                      :: converged = .false.  
  
    converged = .false. 
    alpha     = 0.02   ! fixed step size for now

    do iter=1, maxiter
        ! Calculate the gradient 
        H20 = calcH20(Bogo,h,gaps,blocks)
        N20 = calcN20(Bogo, blocks)

        particles = particle_number_bogo(bogo, blocks)
        !print *, 0, lambda, particles, particles-targetN
        do fermiiter = 1,100
          ! Try this particular value of Lambda
          grad = H20 - lambda * N20
          ! Make a step in the right direction
          newbogo = GradUpdate(grad, bogo, alpha, blocks)

          ! Check the new particle number
          particles = particle_number_bogo(newbogo, blocks)
          !print *, fermiiter, lambda, particles, particles-targetN
          if(abs(Particles-targetN) .lt. 1d-16) exit 
          lambda = lambda - (Particles - targetN)                      
        enddo
        bogo = newbogo

        gradnorm = sqrt(sum(grad**2))

        if(gradnorm .lt. 1d-16) converged = .true.
        if(converged) then
          exit
        endif
    enddo
    H11 = calcH11(Bogo, h, gaps, lambda, blocks)
    call diagonalise_H11(Bogo, H11, blocks, Eqp)
 !   call calc_qp_energies(h, gaps, lambda, bogo, blocks, Eqp, dispqp) 
!    print ('(99f10.3)'), Eqp
    print *, iter, gradnorm, maxval(abs(grad)), maxval(abs(bogo-newbogo)),  &
    &         fermiiter, particles-targetN

    !stop

  end subroutine gradient_step

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

      !print *, 'GRAD'
      !do i=1,T
      !  print ('(99f10.3)'), grad(si+i, si+1:si+T)
      !enddo

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
    
      !allocate(aux(T,T)) 
      !aux =  - alpha**2 * matmul(grad, grad)
      !do i=1, T
      !  aux(i,i) = 1 + aux(i,i) 
      !enddo

      ! Compute the Cholesky factorization
      !call DPOTRF('L', T, aux, T, info )
      !if(info .ne. 0) then
      !    print *, 'Problem with Cholesky decomposition in HFB_gradient.'
      !    print *, 'INFO = ', info
      !    stop
      !endif
      ! Compute the inverse square root
      !call DTRSM('R', 'L', 'N', 'N', T, T, 1.0d0, aux, T, &
      !&                                     newbogo(sb  +1:sb+  T,sb+1:sb+T), T)
      !call DTRSM('R', 'L', 'N', 'N', T, T, 1.0d0, aux, T, &
      !&                                     newbogo(sb+T+1:sb+2*T,sb+1:sb+T), T)

      !deallocate(aux)
   !   sb = sb + 2*T
   ! enddo
  
    ! No need to orthonormalize if using a Cholesky factorization
    !!Don't forget to orthonormalise
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
        Overlap    = 1.0_dp/sqrt(overlap)
        bogo(sb+1:sb+2*T,sb+T+j) = bogo(sb+1:sb+2*T,sb+T+j)*overlap
        ! Orthogonalise all the rest against vector j
        do i=j+1,T
          Overlap = 0.0_dp
          do k=1,2*T
              Overlap = Overlap + bogo(sb+k,sb+T+j) * bogo(sb+k,sb+T+i)
          enddo
          Bogo(:,sb+T+i) = Bogo(:,sb+T+i) - Overlap * Bogo(:,sb+T+j)
        enddo
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

  subroutine calc_qp_energies(h,gaps,lambda,Bogo,blocks,Eqp,dispqp) 
     !--------------------------------------------------------------------------
     !
     !
     !--------------------------------------------------------------------------
     integer, intent(in)        :: blocks(4)
     real(KIND=dp), intent(in)  :: h(:,:), gaps(:,:),  lambda 
     real(KIND=dp), intent(inout) :: Bogo(:,:)
     real(KIND=dp), allocatable :: U(:,:), V(:,:), hV(:,:), hU(:,:), Htot(:,:)
     integer                    :: B, sb, si, i, N, N2, T, j

     real(KIND=dp), intent(out) :: Eqp(:), dispqp(:)
  
     Eqp    = 0.0d0
     dispqp = 0.0d0
     
     si = 0 ; sb = 0
     do B=1,4,2
      N  = blocks(B)    ; if(N.eq.0) cycle 
      N2 = blocks(B+1)

      T = N + N2

      allocate(Htot(2*T, 2*T))
      Htot(1:T,1:T)         =  h(si+1:si+T,si+1:si+T)
      Htot(T+1:2*T,T+1:2*T) = -h(si+1:si+T,si+1:si+T)
      do i=1, T
        Htot(i,i)     = Htot(i,i)     - lambda
        Htot(i+T,i+T) = Htot(i+T,i+T) + lambda
      enddo
      Htot(1:T,T+1:2*T) = gaps(si+1:si+T,si+1:si+T)
      Htot(T+1:2*T,1:T) = gaps(si+1:si+T,si+1:si+T)

      do i=1,2*T
          print('(99f10.3)'), Htot(i,1:2*T)
      enddo
      print*

      Bogo(sb  +1:sb  +T, sb+1:sb+T) =-Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)   
      Bogo(sb+T+1:sb+2*T, sb+1:sb+T) = Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T)   

      do i=1,2*T
          print('(99f10.3)'), Bogo(sb+i,sb+1:sb+2*T)
      enddo
      print*



      Htot = matmul(Htot,(Bogo(sb+1:sb+2*T, sb+1:sb+2*T)))
      Htot = matmul(transpose(Bogo(sb+1:sb+2*T, sb+1:sb+2*T)), Htot)
      do i=1,2*T
          print('(99f10.3)'), Htot(i,1:2*T)
      enddo
      print*
      deallocate(Htot)

      ! Getting the U and V out to make the formulas explicit
      ! and the matrix multiplications memory-local
      U = Bogo(sb  +1:sb+  T,sb+T+1:sb+2*T)
      V = Bogo(sb+T+1:sb+2*T,sb+T+1:sb+2*T)

      !  h V^* and h^t U^*
      hV = matmul(h(si+1:si+T,si+1:si+T), V) - lambda * V
      hU = matmul(h(si+1:si+T,si+1:si+T), U) - lambda * U
      
      ! d^* V^* and dU^*
      hV =-hV +  matmul(gaps(si+1:si+T,si+1:si+T), U)
      hU = hU +  matmul(gaps(si+1:si+T,si+1:si+T), V)

      do i=1, T
        do j=1, T
          Eqp(si+i)   = Eqp(si+i)     +  U(j,i) * hU(j,i) +  V(j,i) * hV(j,i)
          dispqp(si+i)= dispqp(si+i)  + hU(j,i) * hU(j,i) + hV(j,i) * hV(j,i)
        enddo
        dispqp(si+i) = dispqp(si+i) - Eqp(si+i)**2
      enddo

      si = si +  T
      sb = sb +2*T
     enddo
  end subroutine calc_qp_energies

  function particle_number_bogo(bogo, blocks) result(part)
      !-------------------------------------------------------------------------
      !
      !
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in) :: bogo(:,:)
      real(KIND=dp)             :: part
      real(KIND=dp), allocatable:: V(:,:)
      integer, intent(in)       :: blocks(4)

      integer :: B, si, N, N2, T, sb,i,j
      part = 0.0d0
      sb = 0
      do B=1,4,2
        N = blocks(B) ; if(N.eq.0) cycle
        N2= blocks(B+1)
        T = N + N2

        V = Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)
        
        do i=1,T  
          do j=1,T
            part = part + V(i,j) **2
          enddo
        enddo        
        sb = sb + 2 * T
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

      ! Diagonalize this symmetry-subblock
      lwork = -1; allocate(work(1))
      call DSYEV( 'V', 'U', T, A, T, Eqp(si+1:si+T),work,lwork,ifail)
      lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
      call DSYEV( 'V', 'U', T, A, T, Eqp(si+1:si+T),work,lwork,ifail)
      deallocate(work)
      
      !print ('(a3,99f10.3)') , 'EQP', EQP(si+1:si+T)

      ! Transform the U and V matrices
      Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T) = matmul(Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T),A)
      Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T) = matmul(Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T),A)

      si = si +   T
      sb = sb + 2*T
    enddo
  end subroutine diagonalise_H11
end module HFB_gradient
