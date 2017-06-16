module derivatives
 !=======================================================================
 !  #######   ##   #    # #####   ##   #      #    #  ####
 !     #     #  #  ##   #   #    #  #  #      #    # #
 !     #    #    # # #  #   #   #    # #      #    #  ####
 !     #    ###### #  # #   #   ###### #      #    #      #
 !     #    #    # #   ##   #   #    # #      #    # #    #
 !     #    #    # #    #   #   #    # ######  ####   ####
 !
 !  Copyright W. Ryssens & M. Bender
 !
 !=======================================================================
 !
 ! Module that defines the derivatives of functions on the mesh. 
 !
 !
 !=======================================================================
 !
 ! Important CAVEAT: at the moment Tantalus will only allow you to use
 ! symmetry combinations that give rise to 'local' derivatives, i.e. 
 ! symmetry combinations that will let you relate 
 !  
 !   f(-i, j, k) => f(i,j,k)
 !   f( i,-j, k) => f(i,j,k)
 !   f( i, j,-k) => f(i,j,k)
 !
 ! If this is not the case, the matrix multiplications become 'nonlocal'
 ! in the memory-storage meaning of the word, and a naive implementation 
 ! of the matrix multiplication would result in exploding matrix sizes.
 ! 
 !=======================================================================
 !------------------------------------------------------------------------------
 ! To be replaced by Hephaestos  
 !
 ! LINESIZE   : total box size as a function of nx/ny/nz
 ! LINESIZEX $LINESIZEX
 ! LINESIZEY $LINESIZEY
 ! LINESIZEZ $LINESIZEZ
 !
 !------------------------------------------------------------------------------
 
 use compilation
 use geninfo
 
 implicit none
 
 !------------------------------------------------------------------------------
 ! Contains the matrix elements to perform a derivation on the mesh, in either 
 ! the X-, Y- or Z-direction.
 ! 
 ! There are two matrices for every direction. 
 !   When the direction is not affected by any symmetry
 !       A(:,:,1) => derivation matrix
 !       A(:,:,2) => 0, and never used
 ! 
 !   When symmetry is relevant, but derivatives are 'local'
 !       A(:,:,1) => derivation matrix when reflection quantum number is -1
 !       A(:,:,2) => derivation matrix when reflection quantum number is +1
 !             
 !   When symmetry is relevant, but derivatives are 'non-local'
 !       A(:,:,1) => matrix to multiply the vector with
 !       A(:,:,2) => matrix to multiply the symmetric variant with
 !
 !------------------------------------------------------------------------------
 real*8, allocatable ::   derX(:,:,:),  derY(:,:,:),  derZ(:,:,:)
 !------------------------------------------------------------------------------
 ! Contains the matrix elements to perform the operation of second order
 ! (diagonal) derivatives on the mesh. 
 !------------------------------------------------------------------------------
 real*8, allocatable :: laplaX(:,:,:),laplaY(:,:,:),laplaZ(:,:,:)
 
 interface derive_grad
    module procedure derive_grad_1D
    module procedure derive_grad_3D
 end interface
 
 interface derive_lap
    module procedure derive_lap_1D
    module procedure derive_lap_3D
 end interface
 
 interface derive
    module procedure derive_1D
    module procedure derive_3D
 end interface
 
contains 
    
 subroutine inilag
    !---------------------------------------------------------------------------
    ! Computes the Lagrange derivative coefficients for this particular 
    ! symmetry combination.
    !---------------------------------------------------------------------------

    integer       :: i,j, linX, linY, linZ
    real(KIND=dp) :: sinA, A, B, sinB, C, D

    ! Allocate the arrays
    allocate(derX(nx,nx,4), laplaX(nx,nx,4))
    allocate(derY(ny,ny,4), laplaY(ny,ny,4))
    allocate(derZ(nz,nz,4), laplaZ(nz,nz,4))
    
    derX   = 0.0d0 ; derY   = 0.0d0 ; derZ   = 0.0d0
    laplaX = 0.0d0 ; laplaY = 0.0d0 ; laplaZ = 0.0d0
    
    linX = $LINESIZEX
    linY = $LINESIZEY
    linZ = $LINESIZEZ
    
    do i=1,nx
        do j=1,nx
            A           = (pi * (i - j))/linX
            sinA        = sin(A)
            B           = (pi * (i - linX + j-1))/linX 
            sinB        = sin(B)
            
            C = (-1)**(i-j)       *pi/(linX*dx*sinA)
            D = (-1)**(i-linX+j-1)*pi/(linX*dx*sinB)
            
            if(i.eq.j) C = 0
            
            derX(i,j,1) = $DERX_ONE
            derX(i,j,2) = $DERX_TWO
        enddo
    enddo
    
    LaplaX(:,:,1) = matmul(derX(:,:,2),derX(:,:,1))
    LaplaX(:,:,2) = matmul(derX(:,:,1),derX(:,:,2))    
   
    do i=1,ny
        do j=1,ny
            
            A           = (pi * (i - j))/linY
            sinA        = sin(A)
            B           = (pi * (i - linY + j-1))/linY 
            sinB        = sin(B)
            
            C = (-1)**(i-j)       *pi/(linY*dx*sinA)
            D = (-1)**(i-linY+j-1)*pi/(linY*dx*sinB)
                      
            if(i.eq.j) C = 0
            
            derY(i,j,1) = $DERY_ONE
            derY(i,j,2) = $DERY_TWO
        enddo
    enddo
        
    LaplaY(:,:,1) = matmul(derY(:,:,2),derY(:,:,1))
    LaplaY(:,:,2) = matmul(derY(:,:,1),derY(:,:,2))    
    
    do i=1,nz
        do j=1,nz
            A           = (pi * (i - j))/linZ
            sinA        = sin(A)
            B           = (pi * (i - linZ + j-1))/linZ 
            sinB        = sin(B)
            
            C = (-1)**(i-j)       *pi/(linZ*dx*sinA)
            D = (-1)**(i-linZ+j-1)*pi/(linZ*dx*sinB)
            
          
            if(i.eq.j) C = 0
          
            derZ(i,j,1) = $DERZ_ONE
            derZ(i,j,2) = $DERZ_TWO
        enddo
    enddo
        
    LaplaZ(:,:,1) = matmul(derZ(:,:,2),derZ(:,:,1))
    LaplaZ(:,:,2) = matmul(derZ(:,:,1),derZ(:,:,2))    
 end subroutine inilag   
 
 subroutine Derive_3d(f, px, py, pz, fx, fy, fz, df)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the derivative of a function on the mesh.
    ! This routine exists, and combines Derive_grad and Derive_lap for because
    ! many compilers optimize the combination better than both routines 
    ! separately (likely due to cache reusing).
    !
    ! fx = First order derivative in the x direction
    ! fy = First order derivative in the y direction
    ! fz = First order derivative in the z direction
    ! df = Laplacien of the function.
    !
    ! px = sign of the symmetry transformation in the x-direction
    ! py = sign of the symmetry transformation in the y-direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)  :: f(:,:,:)
    real(KIND=dp), intent(out) :: fx(:,:,:), fy(:,:,:), fz(:,:,:), df(:,:,:)
    integer, intent(in)        :: px,py,pz
    
    integer                    :: i,k, sx, sy,sz
    
    sx = (px + 3)/2 ! These are equal to 
    sy = (py + 3)/2 !    1    if pi =   -1  or 0
    sz = (pz + 3)/2 !    2    if pi =   +1 
    
    do i=1,ny*nz
        fx(:,i,1) =                 matmul(derX  (:,:,sx),f(:,i,1))
        df(:,i,1) =                 matmul(laplaX(:,:,sx),f(:,i,1))
    enddo   
    do k=1,nz
        do i=1,nx
            fy(i,:,k) =             matmul(derY  (:,:,sy),f(i,:,k))
            df(i,:,k) = df(i,:,k) + matmul(laplaX(:,:,sy),f(i,:,k))
        enddo
    enddo
    do i=1,nx*ny
        fz(i,1,:) =                 matmul(derZ  (:,:,sz),f(i,1,:))
        df(i,1,:) = df(i,1,:) +     matmul(laplaX(:,:,sz),f(i,1,:))
    enddo
    
 end subroutine Derive_3d

$N2DIAG subroutine Derive_tot(f, px, py, pz, df, ddf)
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG    ! Subroutine that computes the following derivatives on the mesh
$N2DIAG    !
$N2DIAG    ! first order derivatives: x,y,z 
$N2DIAG    ! diagonal second order derivatives :: xx, yy, zz
$N2DIAG    !
$N2DIAG    ! df(:,1)    = First order derivative in the x direction
$N2DIAG    ! df(:,1)    = First order derivative in the y direction
$N2DIAG    ! df(:,1)    = First order derivative in the z direction
$N2DIAG    ! ddf(:,i,j) = Second order derivative in the (i,j) direction. 
$N2DIAG    !
$N2DIAG    ! px = sign of the symmetry transformation in the x-direction
$N2DIAG    ! py = sign of the symmetry transformation in the y-direction
$N2DIAG    ! pz = sign of the symmetry transformation in the z-direction
$N2DIAG    !
$N2DIAG    ! diag = 0   All of the second order derivatives are calculated.
$N2DIAG    ! diag = 1   Only the xx,yy and zz derivatives are calculated.
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG    
$N2DIAG    real(KIND=dp), intent(in)  :: f(:,:,:)
$N2DIAG    real(KIND=dp), intent(out) :: df(:,:,:,:), ddf(:,:,:,:,:)
$N2DIAG    integer, intent(in)        :: px,py,pz
$N2DIAG    
$N2DIAG    integer                    :: i,k, sx, sy,sz
$N2DIAG    
$N2DIAG    sx = (px + 3)/2 ! These are equal to 
$N2DIAG    sy = (py + 3)/2 !    1    if pi =   -1  or 0
$N2DIAG    sz = (pz + 3)/2 !    2    if pi =   +1 
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG    !  First order derivatives and diagonal second-order ones
$N2DIAG    do i=1,ny*nz
$N2DIAG         df(:,i,1,1)   =        matmul(derX  (:,:,sx),f(:,i,1))
$N2DIAG        ddf(:,i,1,1,1) =        matmul(laplaX(:,:,sx),f(:,i,1)) 
$N2DIAG    enddo   
$N2DIAG    do k=1,nz
$N2DIAG        do i=1,nx
$N2DIAG            df(i,:,k,2)    =    matmul(derY  (:,:,sy),f(i,:,k))
$N2DIAG            ddf(i,:,k,2,2) =    matmul(laplaY(:,:,sy),f(i,:,k))                        
$N2DIAG        enddo
$N2DIAG    enddo
$N2DIAG    do i=1,nx*ny
$N2DIAG        df(i,1,:,3)    =        matmul(derZ  (:,:,sz),f(i,1,:))
$N2DIAG        ddf(i,1,:,3,3) =        matmul(laplaZ(:,:,sz),f(i,1,:))
$N2DIAG    enddo
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG end subroutine Derive_tot

$N2ALL subroutine Derive_tot(f, px, py, pz, df, ddf)
$N2ALL    !---------------------------------------------------------------------------
$N2ALL    ! Subroutine that computes the following derivatives on the mesh
$N2ALL    !
$N2ALL    ! first order derivatives: x,y,z 
$N2ALL    ! all second order derivatives :: xx, xy, xz, yx, yy, yz, zx, zy, zz
$N2ALL    !
$N2ALL    ! df(:,1)    = First order derivative in the x direction
$N2ALL    ! df(:,1)    = First order derivative in the y direction
$N2ALL    ! df(:,1)    = First order derivative in the z direction
$N2ALL    ! ddf(:,i,j) = Second order derivative in the (i,j) direction. 
$N2ALL    !
$N2ALL    ! px = sign of the symmetry transformation in the x-direction
$N2ALL    ! py = sign of the symmetry transformation in the y-direction
$N2ALL    ! pz = sign of the symmetry transformation in the z-direction
$N2ALL    !
$N2ALL    ! diag = 0   All of the second order derivatives are calculated.
$N2ALL    ! diag = 1   Only the xx,yy and zz derivatives are calculated.
$N2ALL    !---------------------------------------------------------------------------
$N2ALL    
$N2ALL    real(KIND=dp), intent(in)  :: f(:,:,:)
$N2ALL    real(KIND=dp), intent(out) :: df(:,:,:,:), ddf(:,:,:,:,:)
$N2ALL    integer, intent(in)        :: px,py,pz
$N2ALL    
$N2ALL    integer                    :: i,k, sx, sy,sz
$N2ALL    
$N2ALL    sx = (px + 3)/2 ! These are equal to 
$N2ALL    sy = (py + 3)/2 !    1    if pi =   -1  or 0
$N2ALL    sz = (pz + 3)/2 !    2    if pi =   +1 
$N2ALL    !---------------------------------------------------------------------------
$N2ALL    !  First order derivatives and diagonal second-order ones
$N2ALL    do i=1,ny*nz
$N2ALL         df(:,i,1,1)   =        matmul(derX  (:,:,sx),f(:,i,1))
$N2ALL        ddf(:,i,1,1,1) =        matmul(laplaX(:,:,sx),f(:,i,1)) 
$N2ALL    enddo   
$N2ALL    do k=1,nz
$N2ALL        do i=1,nx
$N2ALL            df(i,:,k,2)    =    matmul(derY  (:,:,sy),f(i,:,k))
$N2ALL            ddf(i,:,k,2,2) =    matmul(laplaY(:,:,sy),f(i,:,k))                        
$N2ALL        enddo
$N2ALL    enddo
$N2ALL    do i=1,nx*ny
$N2ALL        df(i,1,:,3)    =        matmul(derZ  (:,:,sz),f(i,1,:))
$N2ALL        ddf(i,1,:,3,3) =        matmul(laplaZ(:,:,sz),f(i,1,:))
$N2ALL    enddo
$N2ALL    !---------------------------------------------------------------------------
$N2ALL    ! Off-diagonal second order derivatives
$N2ALL    do k=1,nz
$N2ALL      do i=1,nx
$N2ALL          ddf(i,:,k,2,1) =      matmul(derY  (:,:,sy),df(i,:,k,1))
$N2ALL          ddf(i,:,k,1,2) =      ddf(i,:,k,2,1) 
$N2ALL      enddo
$N2ALL    enddo
$N2ALL    
$N2ALL    do i=1,nx*ny
$N2ALL      ddf(i,1,:,3,1) =      matmul(derZ  (:,:,sz),df(i,1,:,1))
$N2ALL      ddf(i,1,:,1,3) =      ddf(i,1,:,3,1) 
$N2ALL      ddf(i,1,:,3,2) =      matmul(derZ  (:,:,sz),df(i,1,:,2))
$N2ALL      ddf(i,1,:,2,3) =      ddf(i,1,:,3,2) 
$N2ALL    enddo
$N2ALL
$N2ALL end subroutine Derive_tot
 
 subroutine Derive_grad_3d(f, px, py, pz, fx, fy, fz)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the gradient of a function on the mesh.
    !
    !
    ! fx = First order derivative in the x direction
    ! fy = First order derivative in the y direction
    ! fz = First order derivative in the z direction

    ! px = sign of the symmetry transformation in the x-direction
    ! py = sign of the symmetry transformation in the y-direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)  :: f(:,:,:)
    real(KIND=dp), intent(out) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    integer, intent(in)        :: px,py,pz
    
    integer                    :: i,j,k, sx, sy,sz
    
    sx = (px + 3)/2 ! These are equal to 
    sy = (py + 3)/2 !    1    if pi =   -1  or 0
    sz = (pz + 3)/2 !    2    if pi =   +1 
    
    do k=1,nz
        do j=1,ny
        fx(:,j,k) =             matmul(derX  (:,:,sx),f(:,j,k))
$DERSYMX        fx(:,j,k) = fx(:,j,k) + matmul(derX  (:,:,sx),f($SYMPARTNERX))
        enddo
    enddo   
    do k=1,nz
        do i=1,nx
            fy(i,:,k) =             matmul(derY  (:,:,sy),f(i,:,k))
$DERSYMY            fy(i,:,k) = fy(i,:,k) + matmul(derX  (:,:,sy),f($SYMPARTNERY))
        enddo
    enddo
    do j=1,ny
        do i=1,nx
        fz(i,j,:) =             matmul(derZ  (:,:,sz),f(i,j,:))
$DERSYMZ        fz(i,j,:) = fz(i,j,:) + matmul(derZ  (:,:,sz),f($SYMPARTNERZ))
        enddo
    enddo
    
 end subroutine Derive_grad_3d
 
 subroutine Derive_X(f, px, fx)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the gradient of a function on the mesh.
    !
    !
    ! fx = First order derivative in the x direction
    ! px = sign of the symmetry transformation in the x-direction
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)  :: f(:,:,:)
    real(KIND=dp), intent(out) :: fx(:,:,:)
    integer, intent(in)        :: px
    
    integer                    :: j,k, sx
    
    sx = (px + 3)/2 
    
    do k=1,nz
        do j=1,ny
        fx(:,j,k) =             matmul(derX  (:,:,sx),f(:,j,k))
$DERSYMX        fx(:,j,k) = fx(:,j,k) + matmul(derX  (:,:,sx),f($SYMPARTNERX))
        enddo
    enddo   
    
 end subroutine Derive_X
 
  subroutine Derive_Y(f, py, fy)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the gradient of a function on the mesh.
    !
    ! fy = First order derivative in the y direction
    ! py = sign of the symmetry transformation in the y-direction
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)  :: f(:,:,:)
    real(KIND=dp), intent(out) :: fy(:,:,:)
    integer, intent(in)        :: py
    integer                    :: i,k, sy
    
    sy = (py + 3)/2 !    1    if pi =   -1  or 0
    
    do k=1,nz
        do i=1,nx
            fy(i,:,k) =             matmul(derY  (:,:,sy),f(i,:,k))
$DERSYMY            fy(i,:,k) = fy(i,:,k) + matmul(derX  (:,:,sy),f($SYMPARTNERY))
        enddo
    enddo
    
 end subroutine Derive_Y
 
  subroutine Derive_Z(f, pz, fz)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the gradient of a function on the mesh.
    ! fz = First order derivative in the z direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)  :: f(:,:,:)
    real(KIND=dp), intent(out) :: fz(:,:,:)
    integer, intent(in)        :: pz
    
    integer                    :: i,j,sz
    
    sz = (pz + 3)/2 !    2    if pi =   +1 
    
    do j=1,ny
        do i=1,nx
        fz(i,j,:) =             matmul(derZ  (:,:,sz),f(i,j,:))
$DERSYMZ        fz(i,j,:) = fz(i,j,:) + matmul(derZ  (:,:,sz),f($SYMPARTNERZ))
        enddo
    enddo
    
 end subroutine Derive_Z
 
 subroutine Derive_lap_3D(f, px, py, pz, df)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the laplacian of a function on the mesh.
    !
    ! fx = First order derivative in the x direction
    ! fy = First order derivative in the y direction
    ! fz = First order derivative in the z direction
    ! df = Laplacien of the function.
    !
    ! px = sign of the symmetry transformation in the x-direction
    ! py = sign of the symmetry transformation in the y-direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)  :: f(:,:,:)
    real(KIND=dp), intent(out) ::  df(:,:,:)
    integer, intent(in)        :: px,py,pz
    
    integer                    :: i,k, sx, sy,sz
    
    sx = (px + 3)/2 ! These are equal to 
    sy = (py + 3)/2 !    1    if pi =   -1  or 0
    sz = (pz + 3)/2 !    2    if pi =   +1 
    
    do i=1,ny*nz
        df(:,i,1) =                 matmul(laplaX(:,:,sx),f(:,i,1))
    enddo   
    do k=1,nz
        do i=1,nx
            df(i,:,k) = df(i,:,k) + matmul(laplaY(:,:,sy),f(i,:,k))
        enddo
    enddo
    do i=1,nx*ny
        df(i,1,:) = df(i,1,:) +     matmul(laplaZ(:,:,sz),f(i,1,:))
    enddo
    
 end subroutine Derive_lap_3D
 
 subroutine Derive_lap_1d(f, px, py, pz, df)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the gradient of a function on the mesh, but on
    ! one that is stored as a vector of nx*ny*nz points.
    !
    ! We use a dirty trick here, by simply reshaping with pointers, which should
    ! avoid copying matrices and not impact the speed. (Let's see in practice.)
    !
    ! fx = First order derivative in the x direction
    ! fy = First order derivative in the y direction
    ! fz = First order derivative in the z direction

    ! px = sign of the symmetry transformation in the x-direction
    ! py = sign of the symmetry transformation in the y-direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in), target  :: f(:)
    real(KIND=dp), intent(out), target :: df(:)
    integer, intent(in)        :: px,py,pz
    real(KIND=dp), pointer     :: f3(:,:,:), df3(:,:,:)
    
    f3(1:nx,1:ny,1:nz)   => f
    df3(1:nx,1:ny,1:nz)  => df
    
    call Derive_lap_3d(f3, px,py,pz,df3)
    
 end subroutine Derive_lap_1d
 
 subroutine Derive_grad_1d(f, px, py, pz, fx, fy, fz)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the gradient of a function on the mesh, but on
    ! one that is stored as a vector of nx*ny*nz points.
    !
    ! We use a dirty trick here, by simply reshaping with pointers, which should
    ! avoid copying matrices and not impact the speed. (Let's see in practice.)
    !
    ! fx = First order derivative in the x direction
    ! fy = First order derivative in the y direction
    ! fz = First order derivative in the z direction

    ! px = sign of the symmetry transformation in the x-direction
    ! py = sign of the symmetry transformation in the y-direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in), target  :: f(:)
    real(KIND=dp), intent(out), target :: fx(:), fy(:), fz(:)
    integer, intent(in)        :: px,py,pz
    real(KIND=dp), pointer     :: f3(:,:,:), fx3(:,:,:), fy3(:,:,:), fz3(:,:,:)
    
    f3 (1:nx,1:ny,1:nz)  => f
    fx3(1:nx,1:ny,1:nz)  => fx
    fy3(1:nx,1:ny,1:nz)  => fy
    fz3(1:nx,1:ny,1:nz)  => fz
    
    call Derive_grad_3d(f3, px,py,pz,fx3, fy3, fz3)
    
 end subroutine Derive_grad_1d
 
 subroutine Derive_1d(f, px, py, pz, fx, fy, fz, df)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the first and second order derivative on the mesh
    !
    ! We use a dirty trick here, by simply reshaping with pointers, which should
    ! avoid copying matrices and not impact the speed. (Let's see in practice.)
    !
    ! fx = First order derivative in the x direction
    ! fy = First order derivative in the y direction
    ! fz = First order derivative in the z direction

    ! px = sign of the symmetry transformation in the x-direction
    ! py = sign of the symmetry transformation in the y-direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in), target  :: f(:)
    real(KIND=dp), intent(out), target :: df(:), fx(:), fy(:), fz(:)
    integer, intent(in)        :: px,py,pz
    real(KIND=dp), pointer     :: f3(:,:,:), df3(:,:,:)
    real(KIND=dp), pointer     :: fx3(:,:,:), fy3(:,:,:), fz3(:,:,:)
    
    f3(1:nx,1:ny,1:nz)   => f
    df3(1:nx,1:ny,1:nz)  => df
    fx3(1:nx,1:ny,1:nz)  => fx
    fy3(1:nx,1:ny,1:nz)  => fy
    fz3(1:nx,1:ny,1:nz)  => fz
    call Derive_3d(f3, px,py,pz,fx3,fy3,fz3,df3)
    
 end subroutine Derive_1d

end module derivatives
