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
 use compilation
 use geninfo
 use wavefunctions, only : HFblocks, nwt
 
 implicit none

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
  real(KIND=dp), intent(inout) :: psi(mv,4,nwt)
  real(KIND=dp), intent(in)    :: transfo(nwt,nwt)
  real(KIND=dp), allocatable   :: temp(:,:,:)

  si      = 0   
  do B=1,8
    N = HFBlocks(B)  ;  if(N .eq. 0) cycle 
    
    allocate(temp(mv,4,N)) 
    temp = 0.0
    do wave1=1,N 
      do wave2=1,N
        temp(:,:,wave1) = temp(:,:,wave1) +                                    &
        &                     Transfo(si+wave2,si+wave1) * psi(:,:,si+wave2) 
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
  real(KIND=dp), intent(in)                :: psi_in(mv,4,nwt),transfo(nwt,nwt)
  real(KIND=dp), intent(out), allocatable  :: psi_out(:,:,:)
  integer                                  :: wave1, wave2, B, N, si
  
  if(.not.allocated(psi_out)) then
    allocate(psi_out(mv,4,nwt))
  endif

  si      = 0   
  psi_out = 0.0
  do B=1,8
    N = HFBlocks(B)  ;  if(N .eq. 0) cycle 
    do wave1=1, N 
      do wave2=1,N
        psi_out(:,:,si+wave1)  = psi_out(:,:,si+wave1) +                       &
        &                     Transfo(si+wave2,si+wave1) * psi_in(:,:,si+wave2) 
      enddo 
    enddo
    si = si +  N
  enddo
 end subroutine transform_spwfs

 function transform_mat(M, transfo) result(Mc)
  !-----------------------------------------------------------------------------
  ! Transform the matrix M with the orthonormal transformation C = transfo.
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
  integer                      :: B, N, si

  si      = 0   
  Mc = 0.0d0
  do B=1,8
    N = HFBlocks(B)  ;  if(N .eq. 0) cycle 

    Mc(si+1:si+N, si+1:si+N) =&
    &             matmul(M(si+1:si+N, si+1:si+N), transfo(si+1:si+N, si+1:si+N))
    Mc(si+1:si+N, si+1:si+N) =&
    & matmul( transpose(transfo(si+1:si+N, si+1:si+N)),Mc(si+1:si+N, si+1:si+N))

    si = si +  N
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
    N = HFBlocks(B)  ;  if(N .eq. 0) cycle 

    Vc(si+1:si+N, si+1:si+N) =&
    & matmul( transpose(transfo(si+1:si+N, si+1:si+N)),V(si+1:si+N, si+1:si+N))

    si = si +  N
  enddo  
  
 end function transform_vec
 
end module basis_transform
