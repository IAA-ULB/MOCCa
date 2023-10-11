module basis_transform
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
 ! Module containing the routines to change from one spwf-basis to the next.
 !
 ! Convention: Here an orthonormal transformation means an orthonormal matrix 
 !             where each column represents a new basis vector, expressed in 
 !             the old basis. 
 !
 !============================================================================== 
 use geninfo
 use wavefunctions, only : HFblocks, nwt, spwf_map, nwt_local, hfblocks_global
 
 implicit none

 !------------------------------------------------------------------------------
 ! Constructing the canonical basis can be demanding, so we include a cutoff
 ! to save some CPU cycles for spwfs with small matrix elements.
 real(KIND=dp), parameter :: basis_cut = 1d-12

contains

 subroutine transform_spwfs_inplace(psi, transfo)
  !-----------------------------------------------------------------------------
  ! Perform a linear transformation of the spwfs, in-place in memory. 
  !
  ! This routine attempts to have the smallest memory-cost possible, performing
  ! the transformation symmetry-block by symmetry-block. This requires a temp
  ! matrix with a non-negligible size. I'm (=W.R.) sure there exists truly
  ! 'in-place' approaches where this cost can be avoided, but I don't know them.
  !
  ! Input:
  !  psi     : input set of spwfs, to be transformed
  !  transfo : unitary transformation C
  !
  ! Output:
  !  psi     : transformed set of spwfs
  !               psi' = C^T psi 
  !-----------------------------------------------------------------------------
  integer                      :: wave1, wave2, B, N, si
  integer                      :: wave1_global, wave2_global
  real(KIND=dp), intent(inout) :: psi(mv,4,nwt_local)
  real(KIND=dp), intent(in)    :: transfo(nwt,nwt)
  real(KIND=dp), allocatable   :: temp(:,:,:)

  si  = 0
  do B=1,8
    N = HFBlocks(B)  ;  if(N .eq. 0) cycle 

    allocate(temp(mv,4,N)) 
    temp = 0.0
    do wave1=1,N    ! The local index of this spwf is si+wave1
      do wave2=1,N  ! The local index of this spwf is si+wave2
        wave1_global = spwf_map(si+wave1) ! Global index
        wave2_global = spwf_map(si+wave2) ! Global index      

        ! Don't bother if the wavefunction is not important enough
        if(abs(Transfo(wave2_global,wave1_global)).lt.basis_cut) cycle
        temp(:,:,wave1) = temp(:,:,wave1) +                                    &
        &                 Transfo(wave2_global,wave1_global) * psi(:,:,si+wave2) 
      enddo 
    enddo
    psi(:,:,si+1:si+N) =  temp
    deallocate(temp)

    si = si +  N
  enddo
 end subroutine transform_spwfs_inplace

 subroutine transform_spwfs(psi_in, psi_out, transfo)
  !-----------------------------------------------------------------------------
  ! Perform a linear transformation of the spwfs in psi_in, producing psi_out.
  !
  ! Input:
  !  psi_in  : input set of spwfs, to be transformed
  !  transfo : unitary transformation C
  !
  ! Output:
  !  psi_out : transformed set of spwfs
  !               psi' = C^T psi 
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in)                :: psi_in(mv,4,nwt_local)
  real(KIND=dp), intent(in)                :: transfo(nwt,nwt)
  real(KIND=dp), intent(out), allocatable  :: psi_out(:,:,:)
  integer                                  :: wave1, wave2, B, N, si
  integer                                  :: wave1_global, wave2_global
 
  if(.not.allocated(psi_out)) then
    allocate(psi_out(mv,4,nwt_local))
  endif

  si      = 0
  psi_out = 0.0
  do B=1,8
    N = HFBlocks(B)  ;  if(N .eq. 0) cycle 
    do wave1=1,N      ! The local index of this spwf is si + wave1
      do wave2=1,N    ! The local index of this spwf is si + wave2

        wave1_global = spwf_map(si+wave1) ! Global index
        wave2_global = spwf_map(si+wave2) ! Global index
        ! Don't bother if the wavefunction is not important enough
        if(abs(Transfo(wave2_global,wave1_global)).lt.basis_cut) cycle

        psi_out(:,:,si+wave1)  = psi_out(:,:,si+wave1) +                       &
        &              Transfo(wave2_global,wave1_global) * psi_in(:,:,si+wave2) 
      enddo 
    enddo
!    do wave1=1,N      ! The local index of this spwf is si + wave1
!      do wave2=1,N    ! The local index of this spwf is si + wave2
!        print *, wave1, wave2, sum(psi_out(:,:,si+wave1)*psi_out(:,:,si+wave2))*dv
!      enddo
!    enddo
!    print *
!    call stp('Canbasis built')
    si = si +  N
  enddo
 end subroutine transform_spwfs

 function transform_mat(M, transfo) result(Mc)
  !-----------------------------------------------------------------------------
  ! Transform the matrix M with the orthonormal transformation C = transfo.
  !
  ! This routine is a bit wasteful: it constructs the new matrix in both 
  ! signature blocks at the same time. When time-reversal is broken, this means
  ! that some 0's are multiplied and added together. Since N/N2/T is generally
  ! small, I don't care enough to do better.
  ! 
  ! Input:
  !     M    : matrix to transform
  !  transfo : unitary transformation C to employ 
  !            (in the conventions of this module)
  !
  ! Output:
  !     Mc = C^T M C
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in)    :: M(nwt,nwt)
  real(KIND=dp), intent(in)    :: transfo(nwt,nwt)
  real(KIND=dp)                :: Mc(nwt,nwt)
  integer                      :: B, N, N2, si, T

  si      = 0   
  Mc = 0.0d0
  do B=1,8,2
    N = HFBlocks_global(B)  ;  if(N .eq. 0) cycle 
    N2= HFBlocks_global(B+1)
    T = N + N2

    Mc(si+1:si+T, si+1:si+T) =&
    &             matmul(M(si+1:si+T, si+1:si+T), transfo(si+1:si+T, si+1:si+T))
    Mc(si+1:si+T, si+1:si+T) =&
    & matmul( transpose(transfo(si+1:si+T, si+1:si+T)),Mc(si+1:si+T, si+1:si+T))

    si = si +  T
  enddo  

 end function transform_mat 

 function transform_vec(V, transfo) result(Vc)
  !-----------------------------------------------------------------------------
  ! Transform a set of vectors V with the orthonormal transformation C=transfo.
  !
  ! Input: 
  !    V     : set of vectors to transform
  !  transfo : unitary transformation C to employ 
  !            (in the conventions of this module)
  ! Output:
  !     Vc = C^T V
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in)    :: V(nwt,nwt)
  real(KIND=dp)                :: Vc(nwt,nwt)
  real(KIND=dp), intent(in)    :: transfo(nwt,nwt)
  integer                      :: B, N, si
  
  si = 0   
  Vc = 0.0d0
  do B=1,8
    N = HFBlocks_global(B)  ;  if(N .eq. 0) cycle 

    Vc(si+1:si+N, si+1:si+N) =&
    & matmul( transpose(transfo(si+1:si+N, si+1:si+N)),V(si+1:si+N, si+1:si+N))

    si = si +  N
  enddo  
  
 end function transform_vec
 
 function transform_diag(H, transfo) result(Hc)
  !-----------------------------------------------------------------------------
  ! 
  !
  ! Input:
  !   H    : matrix 
  ! transfo: unitary transformation C
  !
  ! Output:
  !   Hc   : diagonal matrix elements of 
  !            Hc = C^T H 
  !          (not that Hc is necessarily diagonal)
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in)    :: H(nwt,nwt)
  real(KIND=dp)                :: Hc(nwt)
  real(KIND=dp), intent(in)    :: transfo(nwt,nwt)
  integer                      :: B, N, si, wave1, wave2, wave3
  
  si = 0   
  do B=1,8
    N = HFBlocks_global(B)  ;  if(N .eq. 0) cycle 
    do wave1=1,N
      Hc(si+wave1) = 0
      do wave2=1,N
        do wave3=1,N
         Hc(si+wave1) = Hc(si+wave1) +  Transfo(si+wave2,si+wave1) &
         &                            * Transfo(si+wave3,si+wave1) & 
         &                            * H(si+wave2,si+wave3) 
        enddo
      enddo
    enddo
    si = si +  N
  enddo
 end function transform_diag
 
 function transform_bogo(Bogo, transfo) result(Bc)
  !-----------------------------------------------------------------------------
  ! Input: 
  !    B     : Bogoliubov transformation to transform
  !  transfo : unitary transformation C to employ 
  !            (in the conventions of this module)
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in)    :: Bogo(2*nwt,2*nwt)
  real(KIND=dp)                :: Bc  (2*nwt,2*nwt)
  real(KIND=dp), allocatable   :: A   (    :,    :)
  real(KIND=dp), intent(in)    :: transfo(nwt,nwt)

  integer                      :: B, N, N2, si, sb, T
    
  si = 0
  sb = 0
  Bc = 0.0
  do B=1,8,2
    N = HFBlocks_global(B)  ;  if(N .eq. 0) cycle 
    N2= HFblocks_global(B+1) 

    T = N+N2 
    
    A = transpose(transfo(si  +1:si+  T, si  +1:si  +T))
    
    Bc(sb  +1:sb+  T, sb+T+1:sb+2*T) = &
    &            matmul(A,Bogo(sb  +1:sb+  T, sb+T+1:sb+2*T))
    Bc(sb+T+1:sb+2*T, sb+T+1:sb+2*T) = & 
    &            matmul(A,Bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T))

    Bc(sb  +1:sb+  T, sb+  1:sb+  T) = &
    &            matmul(A, Bogo(sb  +1:sb+  T, sb+  1:sb+  T))
    Bc(sb+T+1:sb+2*T, sb+  1:sb+  T) = & 
    &            matmul(A, Bogo(sb+T+1:sb+2*T, sb+  1:sb+  T))

    si = si +   T
    sb = sb + 2*T
  enddo

 end function transform_bogo
 
end module basis_transform
