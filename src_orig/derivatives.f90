module derivatives
 !==============================================================================
 !  #######   ##   #    # #####   ##   #      #    #  ####
 !     #     #  #  ##   #   #    #  #  #      #    # #
 !     #    #    # # #  #   #   #    # #      #    #  ####
 !     #    ###### #  # #   #   ###### #      #    #      #
 !     #    #    # #   ##   #   #    # #      #    # #    #
 !     #    #    # #    #   #   #    # ######  ####   ####
 !
 !  Copyright W. Ryssens & M. Bender
 !
 !==============================================================================
 !
 ! Module that defines the derivatives of functions on the mesh.
 !
 ! Many changes are made by Hephaestos to this source code.
 !
 ! LINESIZE   : total box size as a function of nx/ny/nz
 ! LINESIZEX  = $LINESIZEX
 ! LINESIZEY  = $LINESIZEY
 ! LINESIZEZ  = $LINESIZEZ
 !
 ! Flags to decide which derivatives of the spwfs get calculated
 !  exclamation marks means they are commented out
 !      Only diagonal second order derivatives
 !      N2DIAG  = $N2DIAG
 !
 !      All second order derivatives (and diagonal third order ones)
 !      N2ALL   = $N2ALL
 !
 !      All third order derivatives
 !      N3ALL   = $N3ALL
 !==============================================================================
 ! Technical notes:
 ! - - - - - - - - - - -
 !
 ! * At the moment Tantalus will only allow you to use symmetry combinations 
 !   that give rise to 'local' derivatives, i.e. symmetry combinations that will
 !   let you relate
 !
 !   f(-i, j, k) => f(i,j,k)
 !   f( i,-j, k) => f(i,j,k)
 !   f( i, j,-k) => f(i,j,k)
 !
 !   If this is not the case, the matrix multiplications become 'nonlocal'
 !   in the memory-storage meaning of the word.
 !
 ! * For historical and readability reasons, the derivatives are all implemented
 !   with respect to 3D functions, i.e. f(nx,ny,nz). However, all wavefunctions
 !   and densities are defined as vectors on the mesh (i.e. f(nx*ny*nz)) for 
 !   speed reasons and to keep the number of indices down. This disparity is 
 !   currently solved using pointer remapping.
 !
 ! * However, the previous point places some constraints on the structure of 
 !   routines in terms of the vectorisation that can be achieved by compilers.
 !   In particular, pointer remapping can "break" intent statements. 
 !   Consider for example:
 ! 
 !     function example(f)
 !        real(KIND=dp),intent(in), target:: f(:)
 !        real(KIND=dp),pointer           :: f3(:,:,:)
 !        f3(1:nx,1:ny,1:nz) => f(1:nx*ny*nz)
 !
 !        [some loop over f3]
 !     end function
 ! 
 !   Compilers will not know that f3 will not be changed during the execution
 !   of the function, even if f was declared as "intent(in)". The gotcha is 
 !   of course that the pointer can still be reassigned ....
 !   
 !   For CRAY compilers at least, and possibly for other compilers as well, 
 !   this kind of structure makes vectorization of the loop impossible. 
 !    
 !   One can avoid this issue by adding an extra layer, i.e. do the pointer
 !   remapping in one routine, which then calls a second routine to actually
 !   perform the loop. This can be vectorized, since now we can declare the 
 !   inputs to the second routine (=the remapped pointers) to be "intent(in)"
 !   themselves, allowing for easy vectorisation! This is the reason the
 !   derive_X/Y/Z functions come in triplets.
 !  
 !
 ! * Another optimisation trick is to avoid any explicit allocatable arrays 
 !   inside "small" routines. If there are such explicit allocates, then 
 !   CRAY compilers (and likely others too) can refuse to inline such routines.
 !   Don't write
 !       A = derX(:,:,sx)
 !       [Some loop involving A]
 !   But rather write the loop explicitly with derX; the compiler will be more
 !   free to inline.
 !
 !==============================================================================
 ! Further thoughts on optimisation, not implemented yet
 ! 
 ! 1. add explicit contiguous statements in more places, perhaps it will help
 !    the compilers optimize.
 ! 2. add explicit call to BLAS routines (DAXPY notably, possibly DGEM)
 !
 !==============================================================================

 use geninfo

 implicit none

 !------------------------------------------------------------------------------
 ! Contains the matrix elements to perform a derivation on the mesh, in either
 ! the X-, Y- or Z-direction.
 !
 ! There are two matrices for every direction.
 !   When the direction is not affected by any symmetry
 !       A(:,:,1) => derivative matrix
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

 interface derive_tot
    module procedure derive_tot_1D
    module procedure derive_tot_3D
 end interface

 interface derive_tot_periodic
    module procedure Derive_tot_periodic_1D
    module procedure Derive_tot_periodic_3D
 end interface

 interface derive_lap
    module procedure derive_lap_1D
    module procedure derive_lap_3D
 end interface

 interface derive_X
    module procedure derive_X_single_1D
    module procedure derive_X_single_3D
    module procedure derive_X_spwf
 end interface

 interface derive_Y
    module procedure derive_Y_single_1D
    module procedure derive_Y_single_3D
    module procedure derive_Y_spwf
 end interface

 interface derive_Z
    module procedure derive_Z_single_1D
    module procedure derive_Z_single_3D
    module procedure derive_Z_spwf
 end interface

contains

 subroutine inilag
    !---------------------------------------------------------------------------
    ! Computes the Lagrange derivative coefficients for this particular
    ! symmetry combination.
    ! NS: modified for correct sign
    !---------------------------------------------------------------------------

    integer       :: i,j, linX, linY, linZ
    real(KIND=dp) :: sinA, A, B, sinB, C, D

#if(USE_Periodic>0)
    real(KIND=dp) :: E, F
#endif

    ! Allocate the arrays
    allocate(derX(nx,nx,4), laplaX(nx,nx,4))
    allocate(derY(ny,ny,4), laplaY(ny,ny,4))
    allocate(derZ(nz,nz,4), laplaZ(nz,nz,4))

    derX   = 0.0d0 ; derY   = 0.0d0 ; derZ   = 0.0d0
    laplaX = 0.0d0 ; laplaY = 0.0d0 ; laplaZ = 0.0d0

    linX = $LINESIZEX
    linY = $LINESIZEY
    linZ = $LINESIZEZ

#if(USE_Periodic==0)
    do i=1,nx
        do j=1,nx
            A           = (pi * (i - j))/linX
            sinA        = sin(A)
            B           = (pi * (i + j-1))/linX
            sinB        = sin(B)

            C = (-1)**(i-j)  *pi/(linX*dx*sinA)
            D = (-1)**(i+j-1)*pi/(linX*dx*sinB) !check the sign

            if(i.eq.j) C = 0

            derX(i,j,1) = $DERX_ONE
            derX(i,j,2) = $DERX_TWO
        enddo
    enddo

    !---------------------------------------------------------------------------
    ! Attention: this construction for the laplacian needs doublechecking for
    ! non-local derivative combinations
    if(linX .eq. 2*nx) then
      ! X-axis is symmetry reduced:
      !     Delta^+_xx = Nabla^-_x \cdot Nabla^+_x
      !     Delta^-_xx = Nabla^+_x \cdot Nabla^-_x
      LaplaX(:,:,1) = matmul(derX(:,:,2),derX(:,:,1))
      LaplaX(:,:,2) = matmul(derX(:,:,1),derX(:,:,2))
    else
      ! X-axis is not symmetry reduced: Delta_xx = Nabla_x * Nabla_x
      LaplaX(:,:,1) = matmul(derX(:,:,1),derX(:,:,1))
      LaplaX(:,:,2) = 0.0
    endif

    do i=1,ny
        do j=1,ny

            A           = (pi * (i - j))/linY
            sinA        = sin(A)
            B           = (pi * (i + j-1))/linY
            sinB        = sin(B)

            C = (-1)**(i-j)       *pi/(linY*dx*sinA)
            D = (-1)**(i+j-1)*pi/(linY*dx*sinB)

            if(i.eq.j) C = 0

            derY(i,j,1) = $DERY_ONE
            derY(i,j,2) = $DERY_TWO
        enddo
    enddo

    if(linY .eq. 2*ny) then
      ! Y-axis is symmetry reduced:
      !     Delta^+_yy = Nabla^-_y \cdot Nabla^+_y
      !     Delta^-_yy = Nabla^+_y \cdot Nabla^-_y
      LaplaY(:,:,1) = matmul(derY(:,:,2),derY(:,:,1))
      LaplaY(:,:,2) = matmul(derY(:,:,1),derY(:,:,2))
    else
      ! Y-axis is not symmetry reduced: Delta_yy = Nabla_y * Nabla_y
      LaplaY(:,:,1) = matmul(derY(:,:,1),derY(:,:,1))
      LaplaY(:,:,2) = 0.0
    endif

    do i=1,nz
        do j=1,nz
            A           = (pi * (i - j))/linZ
            sinA        = sin(A)
            B           = (pi * (i + j-1))/linZ
            sinB        = sin(B)

            C = (-1)**(i-j)       *pi/(linZ*dx*sinA)
            D = (-1)**(i+j-1)*pi/(linZ*dx*sinB)


            if(i.eq.j) C = 0

            derZ(i,j,1) = $DERZ_ONE
            derZ(i,j,2) = $DERZ_TWO
        enddo
    enddo

    if(linZ .eq. 2*nz) then
      ! Z-axis is symmetry reduced:
      !     Delta^+_zz = Nabla^-_z \cdot Nabla^+_z
      !     Delta^-_zz = Nabla^+_z \cdot Nabla^-_z
      LaplaZ(:,:,1) = matmul(derZ(:,:,2),derZ(:,:,1))
      LaplaZ(:,:,2) = matmul(derZ(:,:,1),derZ(:,:,2))
    else
      ! Z-axis is not symmetry reduced: Delta_zz = Nabla_z * Nabla_z
      LaplaZ(:,:,1) = matmul(derZ(:,:,1),derZ(:,:,1))
      LaplaZ(:,:,2) = 0.0
    endif

#else
    !---------------------------------------------------------------------------
    !NS:
    ! This part acts in case of periodic boundary conditions
    ! Computes the Lagrange derivative coefficients for this particular
    ! symmetry combination.
    ! Employs phase shift to make even-N LF functions strictly periodic
    ! So now the derivative matrices are complex
    ! Therefore for reduced symmetries:
    !          der(:,:,1)-real part with even proj
    !          der(:,:,2)-imag part with even proj
    !          der(:,:,3)-ireal part with odd proj
    !          der(:,:,4)-imag part with odd proj
    ! else:
    !          der(:,:,1)-real part
    !          der(:,:,2)-imag part
    !---------------------------------------------------------------------------
    do i=1,nx
        do j=1,nx
            A           = (pi * (i-j))  /linX
            sinA        = sin(A)
            B           = (pi * (i+j-1))/linX
            sinB        = sin(B)

            E = (-1)**(i-j)  *pi/(linX*dx*sinA)
            F = (-1)**(i+j-1)*pi/(linX*dx*sinB)

            !adding phase shifts reals
            C=E*cos(k_shx*(i-j)  *dx)
            D=F*cos(k_shx*(i+j-1)*dx)

            if(i.eq.j) C = 0

            derX(i,j,1) = $DERX_ONE
            derX(i,j,3) = $DERX_TWO

            !adding phase shifts imag
            C=E*sin(k_shx*(i-j)  *dx)
            D=F*sin(k_shx*(i+j-1)*dx)

            if(i.eq.j) C = k_shx

            derX(i,j,2) = $DERX_ONE
            derX(i,j,4) = $DERX_TWO

        enddo
    enddo

    !---------------------------------------------------------------------------
    ! Attention: this construction for the laplacian needs doublechecking for
    ! non-local derivative combinations
    if(linX .eq. 2*nx) then
      ! X-axis is symmetry reduced:
      !e.g., for even wf f(x)=f*(-x):
      !     Re{Delta^+_xx} = Re{Nabla^-_x} \cdot Re{Nabla^+_x} -
      !     Im{Nabla^+_x} \cdot Im{ Nabla^+_x}
      !     Im{Delta^+_xx} = Re{Nabla^+_x} \cdot Im{Nabla^+_x} +
      !     Im{Nabla^-_x} \cdot Re{ Nabla^+_x}
      LaplaX(:,:,1) = matmul(derX(:,:,3),derX(:,:,1)) -                        &
      &               matmul(derX(:,:,2),derX(:,:,2))
      LaplaX(:,:,2) = matmul(derX(:,:,1),derX(:,:,2)) +                        &
      &               matmul(derX(:,:,4),derX(:,:,1))
      LaplaX(:,:,3) = matmul(derX(:,:,1),derX(:,:,3)) -                        &
      &               matmul(derX(:,:,4),derX(:,:,4))
      LaplaX(:,:,4) = matmul(derX(:,:,2),derX(:,:,3)) +                        &
      &               matmul(derX(:,:,3),derX(:,:,4))
    else
      ! X-axis is not symmetry reduced: Delta_xx = Nabla_x * Nabla_x
      LaplaX(:,:,1) = matmul(derX(:,:,1),derX(:,:,1)) -                        &
      &               matmul(derX(:,:,2),derX(:,:,2))
      LaplaX(:,:,2) = matmul(derX(:,:,2),derX(:,:,1)) +                        &
      &               matmul(derX(:,:,1),derX(:,:,2))
      LaplaX(:,:,3) = 0.d0
      LaplaX(:,:,4) = 0.d0
    endif

    do i=1,ny
        do j=1,ny
            A           = (pi * (i - j))  /linY
            sinA        = sin(A)
            B           = (pi * (i + j-1))/linY
            sinB        = sin(B)

            E = (-1)**(i-j)  *pi/(linY*dx*sinA)
            F = (-1)**(i+j-1)*pi/(linY*dx*sinB)

            !adding phase shifts reals
            C=E*cos(k_shy*(i-j)  *dx)
            D=F*cos(k_shy*(i+j-1)*dx)

            if(i.eq.j) C = 0

            derY(i,j,1) = $DERY_ONE
            derY(i,j,3) = $DERY_TWO

            !adding phase shifts imag
            C=E*sin(k_shy*(i-j)  *dx)
            D=F*sin(k_shy*(i+j-1)*dx)

            if(i.eq.j) C = k_shy

            derY(i,j,2) = $DERY_ONE
            derY(i,j,4) = $DERY_TWO

        enddo
    enddo

    !---------------------------------------------------------------------------
    ! Attention: this construction for the laplacian needs doublechecking for
    ! non-local derivative combinations
    if(linY .eq. 2*ny) then
      ! Y-axis is symmetry reduced:
      !e.g., for even wf f(y)=f*(-y):
      !     Re{Delta^+_yy} = Re{Nabla^-_y} \cdot Re{Nabla^+_y} -
      !     Im{Nabla^+_y} \cdot Im{ Nabla^+_y}
      !     Im{Delta^+_yy} = Re{Nabla^+_y} \cdot Im{Nabla^+_y} +
      !     Im{Nabla^-_y} \cdot Re{ Nabla^+_y}
      LaplaY(:,:,1) = matmul(derY(:,:,3),derY(:,:,1)) -                        &
      &               matmul(derY(:,:,2),derY(:,:,2))
      LaplaY(:,:,2) = matmul(derY(:,:,1),derY(:,:,2)) +                        &
      &               matmul(derY(:,:,4),derY(:,:,1))
      LaplaY(:,:,3) = matmul(derY(:,:,1),derY(:,:,3)) -                        &
      &               matmul(derY(:,:,4),derY(:,:,4))
      LaplaY(:,:,4) = matmul(derY(:,:,2),derY(:,:,3)) +                        &
      &               matmul(derY(:,:,3),derY(:,:,4))
    else
      ! X-axis is not symmetry reduced: Delta_xx = Nabla_x * Nabla_x
      LaplaY(:,:,1) = matmul(derY(:,:,1),derY(:,:,1)) -                        &
      &               matmul(derY(:,:,2),derY(:,:,2))
      LaplaY(:,:,2) = matmul(derY(:,:,2),derY(:,:,1)) +                        &
      &               matmul(derY(:,:,1),derY(:,:,2))
      LaplaY(:,:,3) = 0.d0
      LaplaY(:,:,4) = 0.d0
    endif


    do i=1,nz
        do j=1,nz
            A           = (pi * (i - j))  /linZ
            sinA        = sin(A)
            B           = (pi * (i + j-1))/linZ
            sinB        = sin(B)

            E = (-1)**(i-j)  *pi/(linZ*dx*sinA)
            F = (-1)**(i+j-1)*pi/(linZ*dx*sinB)

            !adding phase shifts reals
            C=E*cos(k_shz*(i-j)  *dx)
            D=F*cos(k_shz*(i+j-1)*dx)

            if(i.eq.j) C = 0

            derZ(i,j,1) = $DERZ_ONE
            derZ(i,j,3) = $DERZ_TWO

            !adding phase shifts imag
            C=E*sin(k_shz*(i-j)  *dx)
            D=F*sin(k_shz*(i+j-1)*dx)

            if(i.eq.j) C = k_shz

            derZ(i,j,2) = $DERZ_ONE
            derZ(i,j,4) = $DERZ_TWO

        enddo
    enddo

    !---------------------------------------------------------------------------
    ! Attention: this construction for the laplacian needs doublechecking for
    ! non-local derivative combinations
    if(linZ .eq. 2*nz) then
      ! Z-axis is symmetry reduced:
      !e.g., for even wf f(z)=f*(-z):
      !     Re{Delta^+_zz} = Re{Nabla^-_z} \cdot Re{Nabla^+_z} -
      !     Im{Nabla^+_z} \cdot Im{ Nabla^+_z}
      !     Im{Delta^+_zz} = Re{Nabla^+_z} \cdot Im{Nabla^+_z} +
      !     Im{Nabla^-_z} \cdot Re{ Nabla^+_z}
      LaplaZ(:,:,1) = matmul(derZ(:,:,3),derZ(:,:,1)) -                        &
      &               matmul(derZ(:,:,2),derZ(:,:,2))
      LaplaZ(:,:,2) = matmul(derZ(:,:,1),derZ(:,:,2)) +                        &
      &               matmul(derZ(:,:,4),derZ(:,:,1))
      LaplaZ(:,:,3) = matmul(derZ(:,:,1),derZ(:,:,3)) -                        &
      &               matmul(derZ(:,:,4),derZ(:,:,4))
      LaplaZ(:,:,4) = matmul(derZ(:,:,2),derZ(:,:,3)) +                        &
      &               matmul(derZ(:,:,3),derZ(:,:,4))
    else
      ! X-axis is not symmetry reduced: Delta_xx = Nabla_x * Nabla_x
      LaplaZ(:,:,1) = matmul(derZ(:,:,1),derZ(:,:,1)) -                        &
      &               matmul(derZ(:,:,2),derZ(:,:,2))
      LaplaZ(:,:,2) = matmul(derZ(:,:,2),derZ(:,:,1)) +                        &
      &               matmul(derZ(:,:,1),derZ(:,:,2))
      LaplaZ(:,:,3) = 0.d0
      LaplaZ(:,:,4) = 0.d0
    endif
#endif

 end subroutine inilag

$N2DIAG subroutine Derive_tot_3D(f, px, py, pz, df, ddf)
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
$N2DIAG    ! Note that higher-order derivative tensors are stored in lexicographical order
$N2DIAG    ! in order to cut down on the number of indices and wasted computation.
$N2DIAG    !            1    2    3    4    5    6    7    8    9    10
$N2DIAG    ! 1st order: Dx   Dy   Dz
$N2DIAG    ! 2nd order: Dxx  Dxy  Dxz  Dyy  Dyz  Dzz
$N2DIAG    ! 3rd order: Dxxx Dxxy Dxxz Dxyy Dxyz Dxzz Dyyy Dyyz Dyzz Dzz
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG
$N2DIAG    real(KIND=dp), intent(in)  :: f(:,:,:)
$N2DIAG    real(KIND=dp), intent(out) :: df(:,:,:,:), ddf(:,:,:,:)
$N2DIAG    integer, intent(in)        :: px,py,pz
$N2DIAG
$N2DIAG    integer                    :: i,j,k,l, sx, sy,sz
$N2DIAG    real(KIND=dp), allocatable :: A(:,:), B(:,:)
$N2DIAG
$N2DIAG    sx = (-px + 3)/2 ! These are equal to !NS:check the sign
$N2DIAG    sy = (-py + 3)/2 !    1    if pi =   -1  or 0
$N2DIAG    sz = (-pz + 3)/2 !    2    if pi =   +1
$N2DIAG
$N2DIAG    df = 0.0d0 ; ddf = 0.0d0
$N2DIAG
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG    !  First order derivatives and diagonal second-order ones
$N2DIAG    A = derX(:,:,sx) ; B = laplaX(:,:,sx)
$N2DIAG    do k=1,nz
$N2DIAG     do j=1,ny
$N2DIAG      do l=1,nx
$N2DIAG       do i=1,nx
$N2DIAG         df(i,j,k,1) =  df(i,j,k,1) + A(i,l)*f(l,j,k)
$N2DIAG        ddf(i,j,k,1) = ddf(i,j,k,1) + B(i,l)*f(l,j,k)
$N2DIAG       enddo
$N2DIAG      enddo
$N2DIAG     enddo
$N2DIAG    enddo
$N2DIAG
!           call dgemm('N','N',  nx,ny*nz,nx,1.0d0,A,nx,f(1:nx,1:ny*nz,1),nx,0.0d0, df(1:nx,1:ny*nz,1,1),nx)
!           call dsymm('L','U',  nx,ny*nz,   1.0d0,B,nx,f(1:nx,1:ny*nz,1),nx,0.0d0,ddf(1:nx,1:ny*nz,1,1),nx)

$N2DIAG    A = derY(:,:,sy) ; B = laplaY(:,:,sy)
$N2DIAG    do k=1,nz
$N2DIAG     do j=1,ny
$N2DIAG      do l=1,ny
$N2DIAG       do i=1,nx
$N2DIAG             df(i,j,k,2) =  df(i,j,k,2) + A(j,l)*f(i,l,k)
$N2DIAG            ddf(i,j,k,4) = ddf(i,j,k,4) + B(j,l)*f(i,l,k)
$N2DIAG       enddo
$N2DIAG      enddo
$N2DIAG     enddo
$N2DIAG    enddo
$N2DIAG
!
$N2DIAG    A = derZ(:,:,sz) ; B = laplaZ(:,:,sz)
$N2DIAG    do k=1,nz
$N2DIAG     do l=1,nz
$N2DIAG      do j=1,ny
$N2DIAG        do i=1,nx
$N2DIAG             df(i,j,k,3) =  df(i,j,k,3) + A(k,l)*f(i,j,l)
$N2DIAG            ddf(i,j,k,6) = ddf(i,j,k,6) + B(k,l)*f(i,j,l)
$N2DIAG          enddo
$N2DIAG        enddo
$N2DIAG      enddo
$N2DIAG    enddo
!           call dgemm('N','T',  nx*ny,nz,nz,1.0d0,f(1:nx*ny,1,1:nz),nx*ny,A,nz,0.0d0, df(1:nx*ny,1,1:nz,3),nx*ny)
!           call dsymm('R','U',  nx*ny,   nz,1.0d0,B,nz,f(1:nx*ny,1,1:nz),nx*ny,0.0d0,ddf(1:nx*ny,1,1:nz,6),nx*ny)

$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG    deallocate(A,B)
$N2DIAG end subroutine Derive_tot_3D

$N2DIAG subroutine Derive_tot_1d(f, px, py, pz, df, ddf)
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG    ! Subroutine that computes the gradient of a function on the mesh, but on
$N2DIAG    ! one that is stored as a vector of nx*ny*nz points.
$N2DIAG    !----------------------------------------------------------------------------
$N2DIAG
$N2DIAG    real(KIND=dp), intent(in),target  :: f(:)
$N2DIAG    real(KIND=dp), intent(out),target,contiguous :: df(:,:), ddf(:,:)
$N2DIAG    integer, intent(in)        :: px,py,pz
$N2DIAG    real(KIND=dp), pointer     :: f3(:,:,:), df3(:,:,:,:), ddf3(:,:,:,:)
$N2DIAG
$N2DIAG    f3 (1:nx,1:ny,1:nz)      => f
$N2DIAG    df3(1:nx,1:ny,1:nz,1:3)  => df
$N2DIAG    ddf3(1:nx,1:ny,1:nz,1:6) => ddf
$N2DIAG
$N2DIAG    call Derive_tot_3d(f3, px,py,pz,df3, ddf3)
$N2DIAG
$N2DIAG end subroutine Derive_tot_1d

$N2DIAG subroutine Derive_tot_periodic_3D(f, px, py, pz, df, ddf)
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG    ! Subroutine that computes the following derivatives on the mesh
$N2DIAG    !
$N2DIAG    ! NS:
$N2DIAG    ! first order derivatives: x,y,z
$N2DIAG    ! diagonal second order derivatives :: xx, yy, zz
$N2DIAG    !
$N2DIAG    ! df(:,1,1)    = Real first order derivative in the x direction
$N2DIAG    ! df(:,1,2)    = Imag first order derivative in the x direction
$N2DIAG    ! df(:,2,:)    = First order derivative in the y direction
$N2DIAG    ! df(:,3,:)    = First order derivative in the z direction
$N2DIAG    ! ddf(:,1,:) = Second order derivative in the xx direction.
$N2DIAG    ! ddf(:,4,:) = Second order derivative in the yy direction.
$N2DIAG    ! ddf(:,6,:) = Second order derivative in the zz direction.
$N2DIAG    !
$N2DIAG    ! px = sign of the symmetry transformation in the x-direction
$N2DIAG    ! py = sign of the symmetry transformation in the y-direction
$N2DIAG    ! pz = sign of the symmetry transformation in the z-direction
$N2DIAG    !
$N2DIAG    ! Note that higher-order derivative tensors are stored in lexicographical order
$N2DIAG    ! in order to cut down on the number of indices and wasted computation.
$N2DIAG    !            1    2    3    4    5    6    7    8    9    10
$N2DIAG    ! 1st order: Dx   Dy   Dz
$N2DIAG    ! 2nd order: Dxx  Dxy  Dxz  Dyy  Dyz  Dzz
$N2DIAG    ! 3rd order: Dxxx Dxxy Dxxz Dxyy Dxyz Dxzz Dyyy Dyyz Dyzz Dzz
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG
$N2DIAG    real(KIND=dp), intent(in)  :: f(:,:,:,:)
$N2DIAG    real(KIND=dp), intent(out) :: df(:,:,:,:,:), ddf(:,:,:,:,:)
$N2DIAG    integer, intent(in)        :: px(:),py(:),pz(:)
$N2DIAG
$N2DIAG    integer                    :: i,k,j,l
$N2DIAG    integer, allocatable       :: sx(:), sy(:),sz(:)
$N2DIAG
$N2DIAG    sx = (-px + 1)/2*2 ! These are equal to
$N2DIAG    sy = (-py + 1)/2*2 !    0    if pi =   +1  or 0
$N2DIAG    sz = (-pz + 1)/2*2 !    2    if pi =   -1
$N2DIAG
$N2DIAG    df = 0.0d0 ; ddf = 0.0d0
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG    !  First order derivatives and diagonal second-order ones
$N2DIAG    !NS: Now d_re f=derX_re*f_re-derX_im*f_im
$N2DIAG    !   sx accounts for the parity
$N2DIAG    do k=1,nz
$N2DIAG     do j=1,ny
$N2DIAG      do l=1,nx
$N2DIAG       do i=1,nx
$N2DIAG         df(i,j,k,1,1) = df(i,j,k,1,1) + derX  (i,l,1+sx(1))*f(l,j,k,1) & 
$N2DIAG                       &               - derX  (i,l,2+sx(2))*f(l,j,k,2)
$N2DIAG         df(i,j,k,1,2) = df(i,j,k,1,2) + derX  (i,l,2+sx(1))*f(l,j,k,1) & 
$N2DIAG                       &               + derX  (i,l,1+sx(2))*f(l,j,k,2)
$N2DIAG        ddf(i,j,k,1,1) = ddf(i,j,k,1,1)+ laplaX(i,l,1+sx(1))*f(l,j,k,1) &
$N2DIAG                       &               - laplaX(i,l,2+sx(2))*f(l,j,k,2)
$N2DIAG        ddf(i,j,k,1,2) = ddf(i,j,k,1,2)+ laplaX(i,l,2+sx(1))*f(l,j,k,1) &
$N2DIAG                       &               + laplaX(i,l,1+sx(2))*f(l,j,k,2)
$N2DIAG       enddo
$N2DIAG      enddo
$N2DIAG    enddo
$N2DIAG    enddo
$N2DIAG
$N2DIAG    do k=1,nz
$N2DIAG     do j=1,ny
$N2DIAG      do l=1,ny
$N2DIAG       do i=1,nx
$N2DIAG         df(i,j,k,2,1) = df(i,j,k,2,1) + derY  (j,l,1+sy(1))*f(i,l,k,1) &
$N2DIAG                       &               - derY  (j,l,2+sy(2))*f(i,l,k,2)
$N2DIAG         df(i,j,k,2,2) = df(i,j,k,2,2) + derY  (j,l,2+sy(1))*f(i,l,k,1) &
$N2DIAG                       &               + derY  (j,l,1+sy(2))*f(i,l,k,2)
$N2DIAG        ddf(i,j,k,4,1) = ddf(i,j,k,4,1)+ laplaY(j,l,1+sy(1))*f(i,l,k,1) &
$N2DIAG                       &               - laplaY(j,l,2+sy(2))*f(i,l,k,2)
$N2DIAG        ddf(i,j,k,4,2) = ddf(i,j,k,4,2)+ laplaY(j,l,2+sy(1))*f(i,l,k,1) &
$N2DIAG                       &               + laplaY(j,l,1+sy(2))*f(i,l,k,2)
$N2DIAG       enddo
$N2DIAG      enddo
$N2DIAG     enddo
$N2DIAG    enddo
$N2DIAG
$N2DIAG    do k=1,nz
$N2DIAG     do l=1,nz
$N2DIAG      do j=1,ny
$N2DIAG       do i=1,nx
$N2DIAG         df(i,j,k,3,1) = df(i,j,k,3,1) + derZ  (k,l,1+sz(1))*f(i,j,l,1) &
$N2DIAG                       &               - derZ  (k,l,2+sz(2))*f(i,j,l,2)
$N2DIAG         df(i,j,k,3,2) = df(i,j,k,3,2) + derZ  (k,l,2+sz(1))*f(i,j,l,1) &
$N2DIAG                       &               + derZ  (k,l,1+sz(2))*f(i,j,l,2)
$N2DIAG        ddf(i,j,k,6,1) = ddf(i,j,k,6,1)+ laplaZ(k,l,1+sz(1))*f(i,j,l,1) &
$N2DIAG                       &               - laplaZ(k,l,2+sz(2))*f(i,j,l,2)
$N2DIAG        ddf(i,j,k,6,2) = ddf(i,j,k,6,2)+ laplaZ(k,l,2+sz(1))*f(i,j,l,1) &
$N2DIAG                       &               + laplaZ(k,l,1+sz(2))*f(i,j,l,2)
$N2DIAG       enddo
$N2DIAG      enddo
$N2DIAG     enddo
$N2DIAG    enddo
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG end subroutine Derive_tot_periodic_3D

$N2DIAG subroutine Derive_tot_periodic_1d(f, px, py, pz, df, ddf)
$N2DIAG    !---------------------------------------------------------------------------
$N2DIAG    ! Subroutine that computes the gradient of a function on the mesh, but on
$N2DIAG    ! one that is stored as a vector of nx*ny*nz points.
$N2DIAG    !----------------------------------------------------------------------------
$N2DIAG
$N2DIAG    real(KIND=dp), intent(in) ,target,contiguous  :: f(:,:)
$N2DIAG    real(KIND=dp), intent(out),target,contiguous  :: df(:,:,:)
$N2DIAG    real(KIND=dp), intent(out),target,contiguous  :: ddf(:,:,:)
$N2DIAG    integer, intent(in)        :: px(:),py(:),pz(:)
$N2DIAG    real(KIND=dp), pointer     :: f3(:,:,:,:), df3(:,:,:,:,:)
$N2DIAG    real(KIND=dp), pointer     :: ddf3(:,:,:,:,:)
$N2DIAG
$N2DIAG    f3(1:nx,1:ny,1:nz,1:2)       => f(:,:)
$N2DIAG    df3(1:nx,1:ny,1:nz,1:3,1:2)  => df(:,:,:)
$N2DIAG    ddf3(1:nx,1:ny,1:nz,1:6,1:2) => ddf(:,:,:)
$N2DIAG
$N2DIAG    call Derive_tot_periodic_3d(f3, px,py,pz,df3, ddf3)
$N2DIAG
$N2DIAG end subroutine Derive_tot_periodic_1d

$N2ALL subroutine Derive_tot_1d(f, px, py, pz, df, ddf)
$N2ALL    !---------------------------------------------------------------------------
$N2ALL    ! Subroutine that computes the gradient of a function on the mesh, but on
$N2ALL    ! one that is stored as a vector of nx*ny*nz points.
$N2ALL    !----------------------------------------------------------------------------
$N2ALL    real(KIND=dp), intent(in),target             :: f(:)
$N2ALL    real(KIND=dp), intent(out),target,contiguous :: df(:,:), ddf(:,:)
$N2ALL    integer, intent(in)        :: px,py,pz
$N2ALL    real(KIND=dp), pointer     :: f3(:,:,:), df3(:,:,:,:), ddf3(:,:,:,:)
$N2ALL
$N2ALL    f3 (1:nx,1:ny,1:nz)      => f
$N2ALL    df3(1:nx,1:ny,1:nz,1:3)  => df
$N2ALL    ddf3(1:nx,1:ny,1:nz,1:6) => ddf
$N2ALL
$N2ALL    call Derive_tot_3d(f3, px,py,pz,df3, ddf3)
$N2ALL
$N2ALL end subroutine Derive_tot_1d

$N2ALL subroutine Derive_tot_3D(f, px, py, pz, df, ddf)
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
$N2ALL    ! Note that higher-order derivative tensors are stored in lexicographical order
$N2ALL    ! in order to cut down on the number of indices and wasted computation.
$N2ALL    !            1    2    3    4    5    6    7    8    9    10
$N2ALL    ! 1st order: Dx   Dy   Dz
$N2ALL    ! 2nd order: Dxx  Dxy  Dxz  Dyy  Dyz  Dzz
$N2ALL    ! 3rd order: Dxxx Dxxy Dxxz Dxyy Dxyz Dxzz Dyyy Dyyz Dyzz Dzz
$N2ALL    !---------------------------------------------------------------------------
$N2ALL
$N2ALL    real(KIND=dp), intent(in)  :: f(:,:,:)
$N2ALL    real(KIND=dp), intent(out) :: df(:,:,:,:), ddf(:,:,:,:)
$N2ALL    integer, intent(in)        :: px,py,pz
$N2ALL    real(KIND=dp), allocatable :: A(:,:), B(:,:)
$N2ALL    integer                    :: i,j,k,l, sx, sy,sz
$N2ALL
$N2ALL    sx = (-px + 3)/2 ! These are equal to
$N2ALL    sy = (-py + 3)/2 !    1    if pi =   -1  or 0
$N2ALL    sz = (-pz + 3)/2 !    2    if pi =   +1
$N2ALL    
$N2ALL    df = 0.0d0 ; ddf = 0.0d0
$N2ALL    !---------------------------------------------------------------------------
$N2ALL    !  First order derivatives and diagonal second-order ones
$N2ALL    A = derX(:,:,sx) ; B = laplaX(:,:,sx)
$N2ALL    do k=1,nz
$N2ALL     do j=1,ny
$N2ALL      do l=1,nx
$N2ALL       do i=1,nx
$N2ALL         df(i,j,k,1) =  df(i,j,k,1) + A(i,l)*f(l,j,k)
$N2ALL        ddf(i,j,k,1) = ddf(i,j,k,1) + B(i,l)*f(l,j,k)
$N2ALL       enddo
$N2ALL      enddo
$N2ALL     enddo
$N2ALL    enddo
$N2ALL
$N2ALL
$N2ALL    A = derY(:,:,sy) ; B = laplaY(:,:,sy)
$N2ALL    do k=1,nz
$N2ALL     do j=1,ny
$N2ALL      do l=1,ny
$N2ALL       do i=1,nx
$N2ALL             df(i,j,k,2) =  df(i,j,k,2) + A(j,l)*f(i,l,k)
$N2ALL            ddf(i,j,k,4) = ddf(i,j,k,4) + B(j,l)*f(i,l,k)
$N2ALL       enddo
$N2ALL      enddo
$N2ALL     enddo
$N2ALL    enddo
$N2ALL
$N2ALL
$N2ALL    A = derZ(:,:,sz) ; B = laplaZ(:,:,sz)
$N2ALL    do k=1,nz
$N2ALL     do l=1,nz
$N2ALL      do j=1,ny
$N2ALL        do i=1,nx
$N2ALL             df(i,j,k,3) =  df(i,j,k,3) + A(k,l)*f(i,j,l)
$N2ALL            ddf(i,j,k,6) = ddf(i,j,k,6) + B(k,l)*f(i,j,l)
$N2ALL          enddo
$N2ALL        enddo
$N2ALL      enddo
$N2ALL    enddo
$N2ALL
$N2ALL    !---------------------------------------------------------------------------
$N2ALL    ! Off-diagonal second order derivatives
$N2ALL    A = derY  (:,:,sy)
$N2ALL    do k=1,nz
$N2ALL     do j=1,ny
$N2ALL      do l=1,ny
$N2ALL       do i=1,nx
$N2ALL          ddf(i,j,k,2) = ddf(i,j,k,2) + A(j,l) * df(i,l,k,1)
$N2ALL       enddo
$N2ALL      enddo
$N2ALL     enddo
$N2ALL    enddo
$N2ALL
$N2ALL    A = derZ  (:,:,sz)
$N2ALL    do k=1,nz
$N2ALL     do l=1,nz
$N2ALL      do j=1,ny
$N2ALL       do i=1,nx
$N2ALL          ddf(i,j,k,3) = ddf(i,j,k,3) + A(k,l) * df(i,j,l,1)
$N2ALL          ddf(i,j,k,5) = ddf(i,j,k,3) + A(k,l) * df(i,j,l,2)
$N2ALL       enddo
$N2ALL      enddo
$N2ALL     enddo
$N2ALL    enddo
$N2ALL    deallocate(A,B)
$N2ALL end subroutine Derive_tot_3D

$N3ALL subroutine Derive_tot_3D(f, px, py, pz, df, ddf, dddf)
$N3ALL    !---------------------------------------------------------------------------
$N3ALL    ! Subroutine that computes the following derivatives on the mesh
$N3ALL    !
$N3ALL    ! first order derivatives: x,y,z
$N3ALL    ! all second order derivatives :: xx, xy, xz, yx, yy, yz, zx, zy, zz
$N3ALL    ! all third order derivatives  :: .....
$N3ALL    !
$N3ALL    ! df(:,1)    = First order derivative in the x direction
$N3ALL    ! df(:,1)    = First order derivative in the y direction
$N3ALL    ! df(:,1)    = First order derivative in the z direction
$N3ALL    ! ddf(:,i,j) = Second order derivative in the (i,j) direction.
$N3ALL    ! dddf(:,i,j,k) =  third order derivative (i,j,k)
$N3ALL    !
$N3ALL    ! px = sign of the symmetry transformation in the x-direction
$N3ALL    ! py = sign of the symmetry transformation in the y-direction
$N3ALL    ! pz = sign of the symmetry transformation in the z-direction
$N3ALL    !
$N3ALL    ! Note that higher-order derivative tensors are stored in lexicographical order
$N3ALL    ! in order to cut down on the number of indices and wasted computation.
$N3ALL    !            1    2    3    4    5    6    7    8    9    10
$N3ALL    ! 1st order: Dx   Dy   Dz
$N3ALL    ! 2nd order: Dxx  Dxy  Dxz  Dyy  Dyz  Dzz
$N3ALL    ! 3rd order: Dxxx Dxxy Dxxz Dxyy Dxyz Dxzz Dyyy Dyyz Dyzz Dzz
$N3ALL    !---------------------------------------------------------------------------
$N3ALL
$N3ALL    real(KIND=dp), intent(in)  :: f(:,:,:)
$N3ALL    real(KIND=dp), intent(out) :: df(:,:,:,:), ddf(:,:,:,:)
$N3ALL    real(KIND=dp), intent(out) :: dddf(:,:,:,:)
$N3ALL    integer, intent(in)        :: px,py,pz
$N3ALL
$N3ALL    integer                    :: i,j,k, sx, sy,sz, ax, ay, az
$N3ALL
$N3ALL    sx = (-px + 3)/2 ! These are equal to
$N3ALL    sy = (-py + 3)/2 !    1    if pi =   -1  or 0
$N3ALL    sz = (-pz + 3)/2 !    2    if pi =   +1
$N3ALL    ax = 3 - sx ! These are equal to
$N3ALL    ay = 3 - sy !    2    if si =    1
$N3ALL    az = 3 - sz !    1    if si =    2
$N3ALL    !---------------------------------------------------------------------------
$N3ALL    !  First order derivatives and diagonal second-order ones
$N3ALL    do i=1,ny*nz
$N3ALL               df(:,i,1,1) =    matmul(derX  (:,:,sx),f(:,i,1))
$N3ALL              ddf(:,i,1,1) =    matmul(laplaX(:,:,sx),f(:,i,1))
$N3ALL    enddo
$N3ALL    do k=1,nz
$N3ALL        do i=1,nx
$N3ALL               df(i,:,k,2) =    matmul(derY  (:,:,sy),f(i,:,k))
$N3ALL              ddf(i,:,k,4) =    matmul(laplaY(:,:,sy),f(i,:,k))
$N3ALL        enddo
$N3ALL    enddo
$N3ALL    do j=1,ny
$N3ALL      do i=1,nx
$N3ALL               df(i,j,:,3) =    matmul(derZ  (:,:,sz),f(i,j,:))
$N3ALL              ddf(i,j,:,6) =    matmul(laplaZ(:,:,sz),f(i,j,:))
$N3ALL      enddo
$N3ALL    enddo
$N3ALL    !---------------------------------------------------------------------------
$N3ALL    ! Off-diagonal second order derivatives
$N3ALL    do k=1,nz
$N3ALL      do i=1,nx
$N3ALL              ddf(i,:,k,2) =    matmul(derY  (:,:,sy),df(i,:,k,1))
$N3ALL      enddo
$N3ALL    enddo
$N3ALL
$N3ALL    do j=1,ny
$N3ALL      do i=1,nx
$N3ALL              ddf(i,j,:,3) =    matmul(derZ  (:,:,sz),df(i,j,:,1))
$N3ALL              ddf(i,j,:,5) =    matmul(derZ  (:,:,sz),df(i,j,:,2))
$N3ALL      enddo
$N3ALL    enddo
$N3ALL    !---------------------------------------------------------------------------
$N3ALL    ! Third order derivatives
$N3ALL    do k=1,nz
$N3ALL      do j=1,ny
$N3ALL             dddf(:,j,k, 1)=    matmul(derX  (:,:,sx),ddf(:,j,k,1))
$N3ALL             dddf(:,j,k, 2)=    matmul(derX  (:,:,ax),ddf(:,j,k,2))
$N3ALL             dddf(:,j,k, 3)=    matmul(derX  (:,:,ax),ddf(:,j,k,3))
$N3ALL
$N3ALL             dddf(:,j,k, 4)=    matmul(derX  (:,:,sx),ddf(:,j,k,4))
$N3ALL             dddf(:,j,k, 5)=    matmul(derX  (:,:,sx),ddf(:,j,k,5))
$N3ALL             dddf(:,j,k, 6)=    matmul(derX  (:,:,sx),ddf(:,j,k,6))
$N3ALL      enddo
$N3ALL    enddo
$N3ALL    do k=1,nz
$N3ALL      do i=1,nx
$N3ALL             dddf(i,:,k, 7)=    matmul(derY  (:,:,sy),ddf(i,:,k,4))
$N3ALL             dddf(i,:,k, 8)=    matmul(derY  (:,:,ay),ddf(i,:,k,5))
$N3ALL             dddf(i,:,k, 9)=    matmul(derY  (:,:,sy),ddf(i,:,k,6))
$N3ALL      enddo
$N3ALL    enddo
$N3ALL    do j=1,ny
$N3ALL      do i=1,nx
$N3ALL             dddf(i,j,:,10) =   matmul(derZ  (:,:,sz),ddf(i,j,:,6))
$N3ALL      enddo
$N3ALL    enddo
$N3ALL
$N3ALL end subroutine Derive_tot_3D

$N3ALL subroutine Derive_tot_1D(f, px, py, pz, df, ddf, dddf)
$N3ALL    real(KIND=dp), intent(in),  target, contiguous :: f(:)
$N3ALL    real(KIND=dp), intent(out), target, contiguous :: df(:,:), ddf(:,:)
$N3ALL    real(KIND=dp), intent(out), target, contiguous :: dddf(:,:)
$N3ALL    integer, intent(in)        :: px,py,pz
$N3ALL
$N3ALL    real(KIND=dp), pointer :: f3(:,:,:), df3(:,:,:,:)
$N3ALL    real(KIND=dp), pointer :: ddf3(:,:,:,:), dddf3(:,:,:,:)
$N3ALL
$N3ALL       f3(1:nx, 1:ny, 1:nz)      => f
$N3ALL      df3(1:nx, 1:ny, 1:nz,1:3)  => df
$N3ALL     ddf3(1:nx, 1:ny, 1:nz,1:6)  => ddf
$N3ALL    dddf3(1:nx, 1:ny, 1:nz,1:10) => dddf
$N3ALL
$N3ALL    call Derive_tot_3D(f3, px, py, pz, df3, ddf3, dddf3)
$N3ALL
$N3ALL end subroutine Derive_tot_1D

 subroutine Derive_X_single_3D(f, px, fx)
    !---------------------------------------------------------------------------
    ! Compute the X-derivative of a single function on the mesh.
    !
    ! fx = First order derivative in the x direction
    ! px = sign of the symmetry transformation in the x-direction
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: f(:,:,:)
    real(KIND=dp), intent(out) :: fx(:,:,:)
    integer, intent(in)        :: px
    integer                    :: i,j,k,l,sx

    sx = (-px + 3)/2
    fx  = 0.0d0
    do k=1,nz
     do j=1,ny
      do l=1,nx
       do i=1,nx
        fx(i,j,k) = fx(i,j,k) + derX(i,l,sx) * f(l,j,k)
       enddo
      enddo
     enddo
    enddo
    ! - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Old code line, not updated when loops got rewritten.
    ! $DERSYMX      fx3(:,j,k) = fx3(:,j,k) + matmul(A,f3($SYMPARTNERX))
    ! - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 end subroutine Derive_X_single_3D
 
 subroutine Derive_X_single_1D(f, px, fx)
    !---------------------------------------------------------------------------
    ! Compute the X-derivative of a single 1D-function on the mesh.
    !
    ! fx = First order derivative in the x direction
    ! px = sign of the symmetry transformation in the x-direction
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in), target  :: f(:)
    real(KIND=dp), intent(out),target  :: fx(:)
    integer, intent(in)                :: px
    real(KIND=dp), pointer             :: f3(:,:,:), fx3(:,:,:)
  
    f3 (1:nx,1:ny,1:nz) => f (1:nx*ny*nz)
    fx3(1:nx,1:ny,1:nz) => fx(1:nx*ny*nz)

    call Derive_X_single_3D(f3, px, fx3)
 end subroutine Derive_X_single_1D
 
 subroutine Derive_X_spwf(f, px, fx)
    !---------------------------------------------------------------------------
    ! Compute the X-derivative of four 1D-functions on the mesh.
    !
    ! fx = First order derivative in the x direction
    ! px = sign of the symmetry transformation in the x-direction
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: f(:,:)
    real(KIND=dp), intent(out) :: fx(:,:)
    integer, intent(in)        :: px(4)
    integer                    :: m

    do m=1,4
     call Derive_X_single_1D(f(:,m), px(m), fx(:,m))
    enddo
 end subroutine Derive_X_spwf

  subroutine Derive_Y_single_3D(f, py, fy)
    !---------------------------------------------------------------------------
    ! Compute the Y-derivative of a single 3D-function on the mesh.
    !
    ! fy = First order derivative in the y direction
    ! py = sign of the symmetry transformation in the y-direction
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in), target  :: f(:,:,:)
    real(KIND=dp), intent(out), target :: fy(:,:,:)
    integer, intent(in)                :: py
    integer                            :: i,j,k,l,sy

    sy = (-py + 3)/2 !    1    if pi =   -1  or 0

    fy = 0.0d0
    do k=1,nz
     do j=1,ny
      do l=1,ny
       do i=1,nx
        fy(i,j,k) = fy(i,j,k) + derY(j,l,sy) * f(i,l,k)
       enddo
      enddo
     enddo
    enddo
    ! - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Old code line, not updated when loops got rewritten.
    ! $DERSYMY       fy3(i,:,k) = fy3(i,:,k) + matmul(A,f3($SYMPARTNERY))
    ! - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 end subroutine Derive_Y_single_3D

 subroutine Derive_Y_single_1D(f, py, fy)
    !---------------------------------------------------------------------------
    ! Compute the Y-derivative of a single 1D-function on the mesh.
    !
    ! fy = First order derivative in the y direction
    ! py = sign of the symmetry transformation in the y-direction
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in), target  :: f(:)
    real(KIND=dp), intent(out),target  :: fy(:)
    integer, intent(in)                :: py
    real(KIND=dp), pointer             :: f3(:,:,:), fy3(:,:,:)
  
    f3 (1:nx,1:ny,1:nz) => f (1:nx*ny*nz)
    fy3(1:nx,1:ny,1:nz) => fy(1:nx*ny*nz)

    call Derive_y_single_3D(f3, py, fy3)
 end subroutine Derive_Y_single_1D
 
 subroutine Derive_Y_spwf(f, py, fy)
    !---------------------------------------------------------------------------
    ! Compute the Y-derivative of four 1D-functions on the mesh.
    !
    ! fy = First order derivative in the y direction
    ! py = sign of the symmetry transformation in the y-direction
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in), target  :: f(:,:)
    real(KIND=dp), intent(out), target :: fy(:,:)
    integer, intent(in)                :: py(4)
    integer                            :: m

    do m=1,4
     call Derive_Y_single_1D(f(:,m), py(m), fy(:,m))
    enddo
 end subroutine Derive_Y_spwf

  subroutine Derive_Z_single_3D(f, pz, fz)
    !---------------------------------------------------------------------------
    ! Compute the Z-derivative of a single 1D-function on the mesh.
    ! fz = First order derivative in the z direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: f(:,:,:)
    real(KIND=dp), intent(out) :: fz(:,:,:)
    integer, intent(in)        :: pz
    integer                    :: i,j,k,l,sz

    sz = (-pz + 3)/2 !    2    if pi =   +1
    fz = 0.0d0
    
    do k=1,nz
     do l=1,nz
      do j=1,ny
       do i=1,nx
         fz(i,j,k) = fz(i,j,k) + derZ(k,l,sz) * f(i,j,l)
       enddo
      enddo
     enddo
    enddo
    ! - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Old code line, not updated when loops got rewritten.
    ! $DERSYMZ      fz3(i,j,:) = fz3(i,j,:) + matmul(A,f3($SYMPARTNERZ))
    ! - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 end subroutine Derive_Z_single_3D

 subroutine Derive_Z_single_1D(f, pz, fz)
    !---------------------------------------------------------------------------
    ! Compute the Z-derivative of a single 1D-function on the mesh.
    !
    ! fz = First order derivative in the z direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in), target  :: f(:)
    real(KIND=dp), intent(out),target  :: fz(:)
    integer, intent(in)                :: pz
    real(KIND=dp), pointer             :: f3(:,:,:), fz3(:,:,:)
  
    f3 (1:nx,1:ny,1:nz) => f (1:nx*ny*nz)
    fz3(1:nx,1:ny,1:nz) => fz(1:nx*ny*nz)

    call Derive_z_single_3D(f3, pz, fz3)
 end subroutine Derive_Z_single_1D

 subroutine Derive_Z_spwf(f, pz, fz)
    !---------------------------------------------------------------------------
    ! Compute the Z-derivative of four functions functions on the mesh.
    ! fz = First order derivative in the z direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: f(:,:)
    real(KIND=dp), intent(out) :: fz(:,:)
    integer, intent(in)        :: pz(4)
    integer                    :: m

    do m=1,4
      call Derive_Z_single_1D(f(:,m), pz(m), fz(:,m))
    enddo
    ! - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Old code line, not updated when loops got rewritten.
    ! $DERSYMZ      fz3(i,j,:) = fz3(i,j,:) + matmul(A,f3($SYMPARTNERZ))
    ! - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 end subroutine Derive_Z_spwf

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

    integer                    :: i,j,k,l, sx, sy,sz

    sx = (-px + 3)/2 ! These are equal to
    sy = (-py + 3)/2 !    1    if pi =   -1  or 0
    sz = (-pz + 3)/2 !    2    if pi =   +1

    df = 0.0d0

    do k=1,nz
     do j=1,ny
      do i=1,nx
       do l=1,nx
        df(i,j,k) = df(i,j,k) + laplaX(i,l,sx) * f(l,j,k)
       enddo
      enddo
     enddo
    enddo

    do k=1,nz
     do j=1,ny
      do l=1,ny
       do i=1,nx
        df(i,j,k) = df(i,j,k) + laplaY(j,l,sy) * f(i,l,k)
       enddo
      enddo
     enddo
    enddo

    do k=1,nz
     do l=1,nz
      do j=1,ny
       do i=1,nx
        df(i,j,k) = df(i,j,k) + laplaZ(k,l,sz) * f(i,j,l)
       enddo
      enddo
     enddo
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

 subroutine clean_derivatives()
      if (allocated(derX)) then
        deallocate(derX)
      end if

      if (allocated(derY)) then
        deallocate(derY)
      end if

      if (allocated(derZ)) then
        deallocate(derZ)
      end if

      if (allocated(laplaX)) then
        deallocate(laplaX)
      end if

      if (allocated(laplaY)) then
        deallocate(laplaY)
      end if

      if (allocated(laplaZ)) then
        deallocate(laplaZ)
      end if
 end subroutine clean_derivatives
end module derivatives
