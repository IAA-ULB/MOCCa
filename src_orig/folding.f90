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
 
contains

 pure function Gaussian(r1,r2, r0) result(G)
  !-----------------------------------------------------------------------------
  ! Returns the 1D Gaussian function given by
  ! 
  !     G(|r1 - r2|) = (r0 * sqrt(pi))**(-1) * exp[-(r/r0)**2]
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in) :: r1, r2, r0
  real(KIND=dp) :: dr 
  
  dr = abs(r1-r2)/r0
  G = 1.0/(r0 * sqrt(pi)) * exp(-dr**2)
  return
 end function Gaussian
 
 function FoldGaussian(f,Gx,Gy,Gz) result(Folded)
  !-----------------------------------------------------------------------------
  ! Returns the folded function 
  !
  !  Folded(x2,y2,z2) = int dx1 dy1 dz1 f(x1,y1,z1) 
  !                                Gx(|x1-x2|)Gx(|y1-y2|)Gx(|z1-z2|)
  !
  ! Where the folding functions Gx/Gy/Gz should be constructed beforehand and
  ! passed in. They are in often Gaussians.  
  !-----------------------------------------------------------------------------
    
  real(KIND=dp), intent(in), target :: f(mv)
  real(KIND=dp), intent(in)         :: Gx(nx), Gy(ny), Gz(nz)
  real(KIND=dp), target             :: folded(mv)
  
  real(KIND=dp), pointer    :: f3(:,:,:), folded3(:,:,:)
  !-----------------------------------------------------------------------------
  ! Put the problem in 3D, to make things easier. 
  f3(1:nx,1:ny,1:nz)       => f
  folded3(1:nx,1:ny,1:nz)  => folded
  
  !-----------------------------------------------------------------------------
  ! First convolute/fold along the X-direction
  do k=1,nz
    do j=1,ny
      folded3(:,j,k) = matmul(Gx, f3(1:nx,j,k))
    enddo
  enddo
  ! Then the Y-direction
  do k=1,nz
    do i=1,nx
      folded3(i,:,k) = matmul(Gy, folded3(i,1:ny,k))
    enddo
  enddo
  ! Then the Z-direction
  do j=1,ny
    do i=1,nx
      folded3(i,j,:) = matmul(Gz, folded3(i,j,:))
    enddo
  enddo
  !-----------------------------------------------------------------------------
 end function FoldGaussian
 
end module folding
