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
    ! parameter r0 and symmetry sign p. 
    !
    ! Input:
    !   mesh : mesh coordinates along this direction
    !   m    : number of mesh points along this direction
    !   p    : symmetry sign for reflection along this direction.
    !   r0   : width of the Gaussian
    !
    ! Output:
    !   G  :  Constructed folding matrix
    !
    !          G(i,j) =     (r0 * sqrt(pi))**(-1) * exp[-(|+r_i - r_j|/r0)**2]
    !                 + p * (r0 * sqrt(pi))**(-1) * exp[-(|-r_i - r_j|/r0)**2]
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! A few notes:
    !
    !  *) After construction, we normalize the folding matrix to offset 
    !     numerical errors due to the mesh discretisation. In this way, we 
    !     can guarantee that 
    !
    !         int dx int dx' G(x,x') f(x') = sum_ij G(i,j) f(j) = int dx f(x)  
    !                                                           = sum_i  f(j)
    !     i.e. that we don't change integrals on the mesh.
    ! 
    !  *) If the relevant Cartesian axis is completely stored in the code, 
    !     meaning that there is no reflection symmetry, it is up to the user
    !     to pass in p = 0, removing trivially the contribution of any 
    !     reflection symmetry along the axis.
    !---------------------------------------------------------------------------

    real(KIND=dp)             :: G(m,m)
    integer, intent(in)       ::  m, p
    real(KIND=dp), intent(in) :: r0, mesh(m)
    
    integer :: i,j, ind
    
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
    !---------------------------------------------------------------------------
    !NS: for periodic boundary conditions add contrubution from 2 (symmetric)
    !neighboors. Should be enough for realistic box sizes due to rapid fall down
    !of the exponent.
#if(USE_Periodic==1) 
    do i=1,m
        do j=1,m          
            G(i,j) = G(i,j) + Gaussian(mesh(i)-(1+p)*m*dx, mesh(j), r0)        &
            &      +  Gaussian((1-2*p)*mesh(i)+(1+p)*m*dx, mesh(j), r0)
        enddo
    enddo
#endif
    !---------------------------------------------------------------------------
    ! Normalize, to avoid the numerical errors due to the mesh discretization.
    ! Technical note: we normalize all columns with the norm of ONE PARTICULAR
    !                 column, chosen "sufficiently far away" from the boundary
    !                 of the mesh. If we would normalize G for j = m, on the 
    !                 boundary, we would divide by too small a number, as the 
    !                 Gaussian should extend BEYOND the mesh. 
    !                 Naively, we could choose ind = 1 for this, 
    !                 but this is ON the boundary of the mesh when this axis
    !                 is not reduced through a conserved symmetry. For 
    !                 reasonable meshes and reasonable folding sizes, m/2+1
    !                 is several points away from either boundary of the mesh.
    ind = m/2 + 1 
    
    G(:,:) = G(:,:)/(sum(G(:,ind)*dx))
 end subroutine Gauss_1D
 
 function FoldGaussian(f,Gx,Gy,Gz, mx, my, mz) result(Folded)
  !-----------------------------------------------------------------------------
  ! Returns the folded function 
  !
  !  Folded(x2,y2,z2) = int dx1 dy1 dz1 f(x1,y1,z1) 
  !                                Gx(|x1-x2|)Gx(|y1-y2|)Gx(|z1-z2|)
  !
  ! Where the folding functions Gx/Gy/Gz should be constructed beforehand and
  ! passed in. 
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
      folded(i,j,:) = matmul(Gz, folded(i,j,1:mz))*dx
    enddo
  enddo
  !-----------------------------------------------------------------------------
 end function FoldGaussian
 
end module folding
