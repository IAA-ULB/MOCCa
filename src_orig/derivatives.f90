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
 !
 ! D   : can be used to set variable D to zero: derivatives will not use
 !       symmetries
 ! DX $DX
 ! DY $DY
 ! DZ $DZ
 !------------------------------------------------------------------------------
 
 use compilation
 use geninfo
 
 implicit none
 
 !------------------------------------------------------------------------------
 ! Contains the matrix elements to perform a derivation on the mesh, in either 
 ! the X-, Y- or Z-direction.
 !------------------------------------------------------------------------------
 real*8, allocatable ::   derX(:,:,:),  derY(:,:,:),  derZ(:,:,:)
 !------------------------------------------------------------------------------
 ! Contains the matrix elements to perform the operation of the Laplacian on the
 ! mesh. 
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
    allocate(derX(nx,nx,2), laplaX(nx,nx,2))
    allocate(derY(ny,ny,2), laplaY(ny,ny,2))
    allocate(derZ(nz,nz,2), laplaZ(nz,nz,2))
    
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
            
            D=$DX
                     
            if(i.eq.j) C = 0
            
            derX(i,j,1) = C + D 
            derX(i,j,2) = C - D
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
            
            D=$DY
            
            if(i.eq.j) C = 0
            
            derY(i,j,1) = C + D
            derY(i,j,2) = C - D
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
            
            ! D needs to be set to zero when there is no symmetry in the Z
            ! direction
            D=$DZ
            
            if(i.eq.j) C = 0
          
            derZ(i,j,1) = C + D
            derZ(i,j,2) = C - D 
        enddo
    enddo
        
    LaplaZ(:,:,1) = matmul(derZ(:,:,2),derZ(:,:,1))
    LaplaZ(:,:,2) = matmul(derZ(:,:,1),derZ(:,:,2))    
   
    
!    do j=1,ny
!        do i=1,nx
!            print *, derX(i,j,1), derY(i,j,1), derZ(i,j,1)
!        enddo
!        print *
!    enddo
!    print *, '-----------------------------------------------------'
!    do j=1,ny
!    do i=1,nx
!        print *, derX(i,j,2), derY(i,j,2), derZ(i,j,2)
!    enddo
!    enddo
!    print *

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

 subroutine Derive_tot(f, px, py, pz, df, ddf, diag)
    !---------------------------------------------------------------------------
    ! Subroutine that computes all of the derivatives of a function on the mesh.
    !
    ! df(:,1)    = First order derivative in the x direction
    ! df(:,1)    = First order derivative in the y direction
    ! df(:,1)    = First order derivative in the z direction
    ! ddf(:,i,j) = Second order derivative in the (i,j) direction. 
    !
    ! px = sign of the symmetry transformation in the x-direction
    ! py = sign of the symmetry transformation in the y-direction
    ! pz = sign of the symmetry transformation in the z-direction
    !
    ! diag = 0   All of the second order derivatives are calculated.
    ! diag = 1   Only the xx,yy and zz derivatives are calculated.
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)  :: f(:,:,:)
    real(KIND=dp), intent(out) :: df(:,:,:,:), ddf(:,:,:,:,:)
    integer, intent(in)        :: px,py,pz, diag
    
    integer                    :: i,k, sx, sy,sz
    
    sx = (px + 3)/2 ! These are equal to 
    sy = (py + 3)/2 !    1    if pi =   -1  or 0
    sz = (pz + 3)/2 !    2    if pi =   +1 

    !---------------------------------------------------------------------------
    !  First order derivatives and diagonal second-order ones
    do i=1,ny*nz
         df(:,i,1,1)   =        matmul(derX  (:,:,sx),f(:,i,1))
        ddf(:,i,1,1,1) =        matmul(laplaX(:,:,sx),f(:,i,1)) 
    enddo   
    do k=1,nz
        do i=1,nx
            df(i,:,k,2)    =    matmul(derY  (:,:,sy),f(i,:,k))
            ddf(i,:,k,2,2) =    matmul(laplaY(:,:,sy),f(i,:,k))                        
        enddo
    enddo
    do i=1,nx*ny
        df(i,1,:,3)    =        matmul(derZ  (:,:,sz),f(i,1,:))
        ddf(i,1,:,3,3) =        matmul(laplaZ(:,:,sz),f(i,1,:))
    enddo
    !---------------------------------------------------------------------------
    if(diag.eq.0) then
                ! Other second-order derivatives
    endif   
 end subroutine Derive_tot
 
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
    
    integer                    :: i,k, sx, sy,sz
    
    sx = (px + 3)/2 ! These are equal to 
    sy = (py + 3)/2 !    1    if pi =   -1  or 0
    sz = (pz + 3)/2 !    2    if pi =   +1 
    
    do i=1,ny*nz
        fx(:,i,1) =                 matmul(derX  (:,:,sx),f(:,i,1))
    enddo   
    do k=1,nz
        do i=1,nx
            fy(i,:,k) =             matmul(derY  (:,:,sy),f(i,:,k))
        enddo
    enddo
    do i=1,nx*ny
        fz(i,1,:) =                 matmul(derZ  (:,:,sz),f(i,1,:))
    enddo
    
 end subroutine Derive_grad_3d
 
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
            df(i,:,k) = df(i,:,k) + matmul(laplaX(:,:,sy),f(i,:,k))
        enddo
    enddo
    do i=1,nx*ny
        df(i,1,:) = df(i,1,:) +     matmul(laplaX(:,:,sz),f(i,1,:))
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
