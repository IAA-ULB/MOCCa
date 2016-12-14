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
            
            derX(i,j,1) = C - D 
            derX(i,j,2) = C + D
        enddo
    enddo
   
    do i=1,ny
        do j=1,ny
            if(i .eq. j) cycle
            A           = (pi * (i - j))/linY
            sinA        = sin(A)
            B           = (pi * (i - linZ + j))/linZ 
            sinB        = sin(B)
            
            C = (-1)**(i-j)     *pi/(linY*dx*sinA)
            D = (-1)**(i-linY+j)*pi/(linY*dx*sinB)
            
            D=$DY
            
            if(i.eq.j) C = 0
            
            derY(i,j,1) = C - D
            derY(i,j,2) = C + D
        enddo
    enddo
    
    do i=1,nz
        do j=1,nz
            if(i .eq. j) cycle
            A           = (pi * (i - j))/linZ
            sinA        = sin(A)
            B           = (pi * (i - linZ + j))/linZ 
            sinB        = sin(B)
            
            C = (-1)**(i-j)     *pi/(linZ*dx*sinA)
            D = (-1)**(i-linZ+j)*pi/(linZ*dx*sinB)
            
            ! D needs to be set to zero when there is no symmetry in the Z
            ! direction
            D=$DZ
            
            if(i.eq.j) C = 0
          
            derZ(i,j,1) = C - D
            derZ(i,j,2) = C + D 
        enddo
    enddo

 end subroutine inilag   
 
 
 subroutine Derive(f, px, py, pz, fx, fy, fz, df)
    !---------------------------------------------------------------------------
    ! Subroutine that computes the derivative of a function on the mesh.
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
    do k=1,ny
        do i=1,nx
            fy(i,:,k) =             matmul(derY  (:,:,sy),f(i,:,k))
            df(i,:,k) = df(i,:,k) + matmul(laplaX(:,:,sy),f(i,:,k))
        enddo
    enddo
    do i=1,nx*ny
        fz(i,1,:) =                 matmul(derZ  (:,:,sz),f(i,1,:))
        df(i,1,:) = df(i,1,:) +     matmul(laplaX(:,:,sz),f(i,1,:))
    enddo
    
 end subroutine Derive
 
 
 
end module derivatives
