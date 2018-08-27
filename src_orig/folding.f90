module folding
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
 ! Module that is useful for folding various functions on the mesh. 
 !==============================================================================
 
 use geninfo

 implicit none
 
contains

 pure function Gaussian(r1,r2, r0) result(G)
  !-----------------------------------------------------------------------------
  ! Returns the 1D Gaussian function given by
  ! 
  !     G(|r1 - r2|) = (r0 * sqrt(pi))**(-1) * exp[-(r/r0)**2]
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in) :: r1, r2, r0
  real(KIND=dp) :: dr, G
  
  dr = abs(r1-r2)/(r0)
  G = 1.0/(r0 * sqrt(pi)) * exp(-dr**2)
  return
 end function Gaussian
 
 subroutine Gauss_1D(G, mesh, m, r0, p)
    !---------------------------------------------------------------------------
    ! Function that constructs a matrix to fold in 1-D with a Gaussian of
    ! parameter r0 and symmetry sign p. (To be filled in by Hephaestos later.)
    !---------------------------------------------------------------------------

    real(KIND=dp)             :: G(m,m)
    integer, intent(in)       ::  m, p
    real(KIND=dp), intent(in) :: r0, mesh(m)
    
    integer :: i,j
    
    ! Elements actually represented on the mesh
    do i=1,m
        do j=1,m          
            G(i,j) = Gaussian(mesh(i), mesh(j), r0)
        enddo
    enddo
    !---------------------------------------------------------------------------
    ! Elements to be gotten by symmetry.
    do i=1,m
        do j=1,m          
            G(i,j) = G(i,j) + p*Gaussian(-mesh(i), mesh(j), r0)
        enddo
    enddo

    ! Normalize, to avoid the numerical errors due to the mesh discretization.
    do i=1,m    
        G(:,i) = G(:,i)/(sum(G(:,i)*dx))
    enddo
 end subroutine Gauss_1D
 
 function FoldGaussian(f,Gx,Gy,Gz, mx, my, mz) result(Folded)
  !-----------------------------------------------------------------------------
  ! Returns the folded function 
  !
  !  Folded(x2,y2,z2) = int dx1 dy1 dz1 f(x1,y1,z1) 
  !                                Gx(|x1-x2|)Gx(|y1-y2|)Gx(|z1-z2|)
  !
  ! Where the folding functions Gx/Gy/Gz should be constructed beforehand and
  ! passed in. They are in often Gaussians.  
  !-----------------------------------------------------------------------------
    
  real(KIND=dp), intent(in), target :: f(:,:,:)
  real(KIND=dp), intent(in)         :: Gx(:,:), Gy(:,:), Gz(:,:)
  real(KIND=dp), target             :: folded(mx,my,mz)
  integer                           :: i,j,k
  integer, intent(in)               :: mx,my,mz
  !-----------------------------------------------------------------------------
  ! First convolute/fold along the X-direction
  do k=1,mz
    do j=1,my
      folded(:,j,k) = matmul(Gx, f(1:mx,j,k))*dx
    enddo
  enddo
  ! Then the Y-direction
  do k=1,mz
    do i=1,mx
      folded(i,:,k) = matmul(Gy, folded(i,1:my,k))*dx
    enddo
  enddo
  ! Then the Z-direction
  do j=1,my
    do i=1,mx
      folded(i,j,:) = matmul(Gz, folded(i,j,:))*dx
    enddo
  enddo
  !-----------------------------------------------------------------------------
 end function FoldGaussian
 
end module folding
