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

  subroutine gradient_step(h,gaps,blocks,targetparticles,Bogo,Eqp,lambda,maxiter)
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------  
    integer, intent(in)          :: blocks(4), maxiter
    real(KIND=dp), intent(inout) :: Bogo(:,:), lambda, targetparticles
    real(KIND=dp), intent(inout) :: Eqp(:)
    real(KIND=dp), intent(in)    :: h(:,:), gaps(:,:)
    real(KIND=dp), allocatable   :: H20(:,:), N20(:,:), newbogo(:,:), grad(:,:)
    real(KIND=dp)                :: alpha, particles
    integer                      :: B, si, sb, N, N2, iter

    logical                      :: converged = .false.  
  
    converged = .false. 
    alpha     = 0.02   ! fixed step size for now

    do iter=1, maxiter
    
        ! Calculate the gradient 
        H20 = calcH20(Bogo,h,gaps,blocks)
        N20 = calcN20(Bogo, blocks)

        ! Try this particular value of Lambda
        grad = H20 - lambda * N20
        ! Make a step in the right direction
        newbogo = GradUpdate(grad, bogo, alpha, blocks)
        ! Check the new particle number
        particles = particle_number_bogo(bogo, blocks)
      
    enddo

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
    integer                      :: B, sb, N, N2, T,  i, info

    sb = 0

    newbogo = bogo
    do B=1,4,2
      N  = blocks(B)   ; if(N.eq.0) cycle
      N2 = blocks(B+1)

      T = N + N2

      tU = bogo(sb+  1:  T,sb+1:sb+T)
      tV = bogo(sb+T+1:2*T,sb+1:sb+T)

      ! U = U - alpha V gradient
      newbogo(sb  +1:T     ,sb+1:sb+T) = newbogo(sb  +1:T     ,sb+1:sb+T)      & 
      &                         -alpha *matmul(tV,grad(sb  +1:T     ,sb+1:sb+T)) 
      ! V = V - alpha U gradient
      newbogo(sb+T+1:sb+2*T,sb+1:sb+T) = newbogo(sb+T+1:sb+2*T,sb+1:sb+T)      &
      &                         -alpha *matmul(tU,grad(sb  +1:T     ,sb+1:sb+T)) 


      allocate(aux(T,T)) 
      aux =  - alpha**2 * matmul(grad, grad)
      do i=1, T
        aux(i,i) = 1 + aux(i,i) 
      enddo

      ! Compute the Cholesky factorization
      call DPOTRF('L', T, aux, T, info )
      if(info .ne. 0) then
          print *, 'Problem with Cholesky decomposition in HFB_gradient.'
          print *, 'INFO = ', info
          stop
      endif
      ! Compute the inverse square root
      call DTRSM('R', 'L', 'N', 'N', T, T, 1.0d0, aux, T, &
      &                                     newbogo(sb  +1:sb+  T,sb+1:sb+T), T)
      call DTRSM('R', 'L', 'N', 'N', T, T, 1.0d0, aux, T, &
      &                                     newbogo(sb+T+1:sb+2*T,sb+1:sb+T), T)

      deallocate(aux)
      sb = sb + 2*T
    enddo
  
    ! No need to orthonormalize if using a Cholesky factorization
    !!Don't forget to orthonormalise
    !call ortho_bogo(newbogo, blocks)
    
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
            Overlap = Overlap + Bogo(i,j)**2
        enddo
        Overlap    = 1.0_dp/sqrt(overlap)
        bogo(:,j) = bogo(:,j)*overlap
        ! Orthogonalise all the rest against vector j
        do i=j+1,T
          Overlap = 0.0_dp
          do k=1,2*T
              Overlap = Overlap + bogo(k,j) * bogo(k,i)
          enddo
          Bogo(:,i) = Bogo(:,i) - Overlap * Bogo(:,j)
        enddo
      enddo
      sb = sb+ 2*N+2*N2
    enddo

  end subroutine ortho_bogo

  function calcH20(Bogo,h,gaps,blocks) result(H20)
    !---------------------------------------------------------------------------
    !
    ! H20 =   U^{dagger} h   V^* - V^{\dagger} \Delta^* V^*
    !       - V^{dagger} h^t U^* - U^{\dagger} \Delta   U^*
    !
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: H(:,:), bogo(:,:), gaps(:,:)
    integer, intent(in)        :: blocks(4)
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
      U = Bogo(sb  +1:sb+  T,si+1:si+T)
      V = Bogo(sb+T+1:sb+2*T,si+1:si+T)

      !  h V^* and h^t U^*
      hV = matmul(   h(si+1:si+T,si+1:si+T), V)
      hU = matmul(   h(si+1:si+T,si+1:si+T), U)

      ! d^* V^* and dU^*
      dV = matmul(gaps(si+1:si+T,si+1:si+T), V)
      dU = matmul(gaps(si+1:si+T,si+1:si+T), U)

      ! We reuse the defined symbols to save a matrix multiplication here
      U = transpose(U) ; V = transpose(V)
      hV = hV + dU
      hU = hU + dV
      ! We can save some effort here in the future, H20 is antisymmetric     
      H20(si+1:si+T, si+1:si+T) = matmul(U, hV) - matmul(V,hU)

      si = si +  T
      sb = sb +2*T
    enddo

  end function calcH20

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

    sb = 0
    do B=1,4,2
      N = blocks(B) ; if(N.eq.0) cycle
      N2= blocks(B+1)
      T = N + N2

      U = Bogo(sb  +1:sb+  T, sb+1:sb+T)
      V = Bogo(sb+T+1:sb+2*T, sb+1:sb+T)

      N20(sb+1:sb+T,sb+1:sb+T) = matmul(transpose(U),V)-matmul(transpose(V),U)
        
      sb = sb + 2 * T
    enddo
  end function calcN20

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

        V = Bogo(sb+T+1:sb+2*T, sb+1:sb+T)
        
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

end module HFB_gradient
