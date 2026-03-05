module sphamil

  implicit none

  real*8, private, allocatable :: derX(:,:,:), derY(:,:,:), derZ(:,:,:)

contains

  subroutine apply_sphamil( psi, dpsi_x, dpsi_y, dpsi_z, ddpsi_xx, ddpsi_yy, ddpsi_zz, &
       &                    sx, sy, sz, iso, F_I_I, F_Nm_Nm, G_I_NS, hbm, nx, ny, nz, dx, ngrid, hpsi )
    !-----------------------------------------------------------------------------------
    ! Apply the single-particle hamiltonian to a single-particle wavefunction.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    ! Input:
    !    psi :       real*8 array, dimension (ngrid,4)
    !                a single-particle wavefunction
    !
    !    dpsi_x/y/z: real*8 arrays, dimension ( ngrid, 4 )
    !                first order derivatives of \psi
    !
    !    ddpsi_xx/yy/zz: real*8 arrays, dimension (ngrid, 4)
    !                    diagonal second order derivatives of \psi
    !
    !    sx/sz/sz      : integer array, dimension (4)
    !                    the symmetry properties (i.e. +-1) of the four components of psi
    !
    !    iso           : integer, -1 (neutron) or +1 (proton)
    !                    the isospin quantum number of psi
    !
    !    F_I_I         : real*8 array, dimension (ngrid,2)
    !                    the mean-field potential F_I_I
    !                    ATTENTION: last index is neutron (1), proton(2)
    !
    !   F_Nm_Nm        : real*8 array, dimension (ngrid,2)
    !                    the kinetic mean-field potential F_Nm_Nm
    !                    ATTENTION: last index is neutron (1), proton(2)
    !
    !   G_I_NS         : real*8 array, dimension (ngrid,3,3,2)
    !                    the spin-current mean-field potential G_I_NS
    !                    ATTENTION: last index is neutron (1), proton(2)
    !
    !   hbm            : real*8 array, dimension (2)
    !                    numerical values of \hbar^2/( 2 m) for the kinetic energy
    !                    the index is neutron (1), proton(2)
    !
    !   nx,ny,nz       : integer
    !                    the number of grid points in each direction
    !   dx             : real*8
    !                    mesh spacing in units of fm
    !
    !   ngrid          : integer
    !                    the number of grid points
    !
    ! Output:
    !   hpsi           : real*8 array, dimension (ngrid,4)
    !                    the application of the single-particle hamiltonian to psi, i.e.
    !                     [ h \psi ] (r, \sigma)
    !------------------------------------------------------------------------------------

    implicit none

    ! Mesh quantities
    integer, intent(in) :: nx,ny,nz, ngrid
    real*8, intent(in)  :: dx
    ! Wavefunction quantities
    real*8, intent(in)  :: psi(ngrid,4)
    real*8, intent(in)  :: dpsi_x(ngrid,4), dpsi_y(ngrid,4), dpsi_z(ngrid,4)
    real*8, intent(in)  :: ddpsi_xx(ngrid,4), ddpsi_yy(ngrid,4), ddpsi_zz(ngrid,4)
    real*8, intent(out) :: hpsi(ngrid,4)

    integer, intent(in) :: sx(4),sy(4),sz(4)
    integer, intent(in) :: iso

    ! Potentials
    real*8, intent(in)  :: F_I_I(ngrid,2), F_Nm_Nm(ngrid,2), G_I_NS(ngrid,3,3,2)

    ! Interaction options
    real*8, intent(in)  :: hbm(2)

    ! Local variables
    integer :: k, i, it
    real*8  :: temp(ngrid,4)
    integer :: sym(4)
    ! Technical note: dtemp declared with the spinor indices (1-4)
    !                 BEFORE the derivative indices (3). This is to aid the
    !                 memory locality of operations in this particular function
    !                 and is OPPOSITE the conventions of the rest of the code.
    real*8  :: dtemp(ngrid,4,3)

    if(.not. allocated(derX)) call inilag(nx,ny,nz, dx)

    ! Go from isospin  = -1 / +1 to array indices 1, 2
    it = (iso + 3) / 2

    ! Action of the kinetic energy operator
    do k=1,4
       do i=1,ngrid
          hpsi(i,k) = -hbm(it) * (ddpsi_xx(i,k) + ddpsi_yy(i,k) + ddpsi_zz(i,k))
       enddo
    enddo

    !----------------------------------------------------------------------------
    ! Action of the operator corresponding to F_I_I
    temp = 0.0d0
    do i=1,ngrid
       temp(i,1) =  temp(i,1) +  F_I_I(i,it) * psi(i,1)
       temp(i,2) =  temp(i,2) +  F_I_I(i,it) * psi(i,2)
       temp(i,3) =  temp(i,3) +  F_I_I(i,it) * psi(i,3)
       temp(i,4) =  temp(i,4) +  F_I_I(i,it) * psi(i,4)
    enddo
    hpsi(:,1) =  hpsi(:,1) +  temp(:,1)
    hpsi(:,2) =  hpsi(:,2) +  temp(:,2)
    hpsi(:,3) =  hpsi(:,3) +  temp(:,3)
    hpsi(:,4) =  hpsi(:,4) +  temp(:,4)

    !----------------------------------------------------------------------------
    ! Action of F_Nm_Nm symmetrized: 0
    temp = 0.0d0
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  F_Nm_Nm(i,it) * dpsi_x(i,1)
       temp(i,2) =  temp(i,2) -  F_Nm_Nm(i,it) * dpsi_x(i,2)
       temp(i,3) =  temp(i,3) -  F_Nm_Nm(i,it) * dpsi_x(i,3)
       temp(i,4) =  temp(i,4) -  F_Nm_Nm(i,it) * dpsi_x(i,4)
    enddo
    sym(1) = -sX(1)
    sym(2) = -sX(2)
    sym(3) = -sX(3)
    sym(4) = -sX(4)
    call Derive_X(temp(:,:), sym, dtemp(:,:,1),nx,ny,nz)
    hpsi(:,1) =  hpsi(:,1) +  dtemp(:,1,1)
    hpsi(:,2) =  hpsi(:,2) +  dtemp(:,2,1)
    hpsi(:,3) =  hpsi(:,3) +  dtemp(:,3,1)
    hpsi(:,4) =  hpsi(:,4) +  dtemp(:,4,1)

    temp = 0.0
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  F_Nm_Nm(i,it) * dpsi_y(i,1)
       temp(i,2) =  temp(i,2) -  F_Nm_Nm(i,it) * dpsi_y(i,2)
       temp(i,3) =  temp(i,3) -  F_Nm_Nm(i,it) * dpsi_y(i,3)
       temp(i,4) =  temp(i,4) -  F_Nm_Nm(i,it) * dpsi_y(i,4)
    enddo
    sym(1) = -sY(1)
    sym(2) = -sY(2)
    sym(3) = -sY(3)
    sym(4) = -sY(4)
    call Derive_Y(temp(:,:), sym, dtemp(:,:,2),nx,ny,nz)
    hpsi(:,1) =  hpsi(:,1) +  dtemp(:,1,2)
    hpsi(:,2) =  hpsi(:,2) +  dtemp(:,2,2)
    hpsi(:,3) =  hpsi(:,3) +  dtemp(:,3,2)
    hpsi(:,4) =  hpsi(:,4) +  dtemp(:,4,2)

    temp = 0.0
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  F_Nm_Nm(i,it) * dpsi_z(i,1)
       temp(i,2) =  temp(i,2) -  F_Nm_Nm(i,it) * dpsi_z(i,2)
       temp(i,3) =  temp(i,3) -  F_Nm_Nm(i,it) * dpsi_z(i,3)
       temp(i,4) =  temp(i,4) -  F_Nm_Nm(i,it) * dpsi_z(i,4)
    enddo
    sym(1) = -sZ(1)
    sym(2) = -sZ(2)
    sym(3) = -sZ(3)
    sym(4) = -sZ(4)
    call Derive_Z(temp(:,:), sym, dtemp(:,:,3),nx,ny,nz)
    hpsi(:,1) =  hpsi(:,1) +  dtemp(:,1,3)
    hpsi(:,2) =  hpsi(:,2) +  dtemp(:,2,3)
    hpsi(:,3) =  hpsi(:,3) +  dtemp(:,3,3)
    hpsi(:,4) =  hpsi(:,4) +  dtemp(:,4,3)

    !---------------------------------------------------------------------------
    ! First part of the action of G_I_NS
    !   0.5 * i * \sum_{mk}  G_I_NmSk \nabla_m \sigma_k \psi
    temp = 0.0d0
    do i=1,ngrid
       temp(i,1) =  temp(i,1) +  G_I_NS(i,1,1,it) * dpsi_x(i,4)
       temp(i,2) =  temp(i,2) -  G_I_NS(i,1,1,it) * dpsi_x(i,3)
       temp(i,3) =  temp(i,3) +  G_I_NS(i,1,1,it) * dpsi_x(i,2)
       temp(i,4) =  temp(i,4) -  G_I_NS(i,1,1,it) * dpsi_x(i,1)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  G_I_NS(i,1,2,it) * dpsi_x(i,3)
       temp(i,2) =  temp(i,2) -  G_I_NS(i,1,2,it) * dpsi_x(i,4)
       temp(i,3) =  temp(i,3) +  G_I_NS(i,1,2,it) * dpsi_x(i,1)
       temp(i,4) =  temp(i,4) +  G_I_NS(i,1,2,it) * dpsi_x(i,2)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) +  G_I_NS(i,1,3,it) * dpsi_x(i,2)
       temp(i,2) =  temp(i,2) -  G_I_NS(i,1,3,it) * dpsi_x(i,1)
       temp(i,3) =  temp(i,3) -  G_I_NS(i,1,3,it) * dpsi_x(i,4)
       temp(i,4) =  temp(i,4) +  G_I_NS(i,1,3,it) * dpsi_x(i,3)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) +  G_I_NS(i,2,1,it) * dpsi_y(i,4)
       temp(i,2) =  temp(i,2) -  G_I_NS(i,2,1,it) * dpsi_y(i,3)
       temp(i,3) =  temp(i,3) +  G_I_NS(i,2,1,it) * dpsi_y(i,2)
       temp(i,4) =  temp(i,4) -  G_I_NS(i,2,1,it) * dpsi_y(i,1)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  G_I_NS(i,2,2,it) * dpsi_y(i,3)
       temp(i,2) =  temp(i,2) -  G_I_NS(i,2,2,it) * dpsi_y(i,4)
       temp(i,3) =  temp(i,3) +  G_I_NS(i,2,2,it) * dpsi_y(i,1)
       temp(i,4) =  temp(i,4) +  G_I_NS(i,2,2,it) * dpsi_y(i,2)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) +  G_I_NS(i,2,3,it) * dpsi_y(i,2)
       temp(i,2) =  temp(i,2) -  G_I_NS(i,2,3,it) * dpsi_y(i,1)
       temp(i,3) =  temp(i,3) -  G_I_NS(i,2,3,it) * dpsi_y(i,4)
       temp(i,4) =  temp(i,4) +  G_I_NS(i,2,3,it) * dpsi_y(i,3)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) +  G_I_NS(i,3,1,it) * dpsi_z(i,4)
       temp(i,2) =  temp(i,2) -  G_I_NS(i,3,1,it) * dpsi_z(i,3)
       temp(i,3) =  temp(i,3) +  G_I_NS(i,3,1,it) * dpsi_z(i,2)
       temp(i,4) =  temp(i,4) -  G_I_NS(i,3,1,it) * dpsi_z(i,1)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  G_I_NS(i,3,2,it) * dpsi_z(i,3)
       temp(i,2) =  temp(i,2) -  G_I_NS(i,3,2,it) * dpsi_z(i,4)
       temp(i,3) =  temp(i,3) +  G_I_NS(i,3,2,it) * dpsi_z(i,1)
       temp(i,4) =  temp(i,4) +  G_I_NS(i,3,2,it) * dpsi_z(i,2)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) +  G_I_NS(i,3,3,it) * dpsi_z(i,2)
       temp(i,2) =  temp(i,2) -  G_I_NS(i,3,3,it) * dpsi_z(i,1)
       temp(i,3) =  temp(i,3) -  G_I_NS(i,3,3,it) * dpsi_z(i,4)
       temp(i,4) =  temp(i,4) +  G_I_NS(i,3,3,it) * dpsi_z(i,3)
    enddo
    hpsi(:,1) =  hpsi(:,1) +  0.5d0 *  temp(:,1)
    hpsi(:,2) =  hpsi(:,2) +  0.5d0 *  temp(:,2)
    hpsi(:,3) =  hpsi(:,3) +  0.5d0 *  temp(:,3)
    hpsi(:,4) =  hpsi(:,4) +  0.5d0 *  temp(:,4)

    !---------------------------------------------------------------------------
    ! Action of G_I_NS symmetrized: -1
    !  - 0.5 *   \sum_{mk} \nabla_m G_I_NmSk \sigma_k \psi
    temp = 0.0d0
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  G_I_NS(i,1,1,it) * psi(i,4)
       temp(i,2) =  temp(i,2) +  G_I_NS(i,1,1,it) * psi(i,3)
       temp(i,3) =  temp(i,3) -  G_I_NS(i,1,1,it) * psi(i,2)
       temp(i,4) =  temp(i,4) +  G_I_NS(i,1,1,it) * psi(i,1)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) +  G_I_NS(i,1,2,it) * psi(i,3)
       temp(i,2) =  temp(i,2) +  G_I_NS(i,1,2,it) * psi(i,4)
       temp(i,3) =  temp(i,3) -  G_I_NS(i,1,2,it) * psi(i,1)
       temp(i,4) =  temp(i,4) -  G_I_NS(i,1,2,it) * psi(i,2)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  G_I_NS(i,1,3,it) * psi(i,2)
       temp(i,2) =  temp(i,2) +  G_I_NS(i,1,3,it) * psi(i,1)
       temp(i,3) =  temp(i,3) +  G_I_NS(i,1,3,it) * psi(i,4)
       temp(i,4) =  temp(i,4) -  G_I_NS(i,1,3,it) * psi(i,3)
    enddo
    sym(1) = -sX(1)
    sym(2) = -sX(2)
    sym(3) = -sX(3)
    sym(4) = -sX(4)
    call Derive_X(temp(:,:), sym, dtemp(:,:,1),nx,ny,nz)
    hpsi(:,1) =  hpsi(:,1) -  0.5d0 *  dtemp(:,1,1)
    hpsi(:,2) =  hpsi(:,2) -  0.5d0 *  dtemp(:,2,1)
    hpsi(:,3) =  hpsi(:,3) -  0.5d0 *  dtemp(:,3,1)
    hpsi(:,4) =  hpsi(:,4) -  0.5d0 *  dtemp(:,4,1)

    temp = 0.0d0
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  G_I_NS(i,2,1,it) * psi(i,4)
       temp(i,2) =  temp(i,2) +  G_I_NS(i,2,1,it) * psi(i,3)
       temp(i,3) =  temp(i,3) -  G_I_NS(i,2,1,it) * psi(i,2)
       temp(i,4) =  temp(i,4) +  G_I_NS(i,2,1,it) * psi(i,1)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) +  G_I_NS(i,2,2,it) * psi(i,3)
       temp(i,2) =  temp(i,2) +  G_I_NS(i,2,2,it) * psi(i,4)
       temp(i,3) =  temp(i,3) -  G_I_NS(i,2,2,it) * psi(i,1)
       temp(i,4) =  temp(i,4) -  G_I_NS(i,2,2,it) * psi(i,2)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  G_I_NS(i,2,3,it) * psi(i,2)
       temp(i,2) =  temp(i,2) +  G_I_NS(i,2,3,it) * psi(i,1)
       temp(i,3) =  temp(i,3) +  G_I_NS(i,2,3,it) * psi(i,4)
       temp(i,4) =  temp(i,4) -  G_I_NS(i,2,3,it) * psi(i,3)
    enddo
    sym(1) = -sY(1)
    sym(2) = -sY(2)
    sym(3) = -sY(3)
    sym(4) = -sY(4)
    call Derive_Y(temp(:,:), sym, dtemp(:,:,2),nx,ny,nz)
    hpsi(:,1) =  hpsi(:,1) -  0.5d0 *  dtemp(:,1,2)
    hpsi(:,2) =  hpsi(:,2) -  0.5d0 *  dtemp(:,2,2)
    hpsi(:,3) =  hpsi(:,3) -  0.5d0 *  dtemp(:,3,2)
    hpsi(:,4) =  hpsi(:,4) -  0.5d0 *  dtemp(:,4,2)

    temp = 0.0d0
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  G_I_NS(i,3,1,it) * psi(i,4)
       temp(i,2) =  temp(i,2) +  G_I_NS(i,3,1,it) * psi(i,3)
       temp(i,3) =  temp(i,3) -  G_I_NS(i,3,1,it) * psi(i,2)
       temp(i,4) =  temp(i,4) +  G_I_NS(i,3,1,it) * psi(i,1)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) +  G_I_NS(i,3,2,it) * psi(i,3)
       temp(i,2) =  temp(i,2) +  G_I_NS(i,3,2,it) * psi(i,4)
       temp(i,3) =  temp(i,3) -  G_I_NS(i,3,2,it) * psi(i,1)
       temp(i,4) =  temp(i,4) -  G_I_NS(i,3,2,it) * psi(i,2)
    enddo
    do i=1,ngrid
       temp(i,1) =  temp(i,1) -  G_I_NS(i,3,3,it) * psi(i,2)
       temp(i,2) =  temp(i,2) +  G_I_NS(i,3,3,it) * psi(i,1)
       temp(i,3) =  temp(i,3) +  G_I_NS(i,3,3,it) * psi(i,4)
       temp(i,4) =  temp(i,4) -  G_I_NS(i,3,3,it) * psi(i,3)
    enddo
    sym(1) = -sZ(1)
    sym(2) = -sZ(2)
    sym(3) = -sZ(3)
    sym(4) = -sZ(4)
    call Derive_Z(temp(:,:), sym, dtemp(:,:,3),nx,ny,nz)
    hpsi(:,1) =  hpsi(:,1) -  0.5d0 *  dtemp(:,1,3)
    hpsi(:,2) =  hpsi(:,2) -  0.5d0 *  dtemp(:,2,3)
    hpsi(:,3) =  hpsi(:,3) -  0.5d0 *  dtemp(:,3,3)
    hpsi(:,4) =  hpsi(:,4) -  0.5d0 *  dtemp(:,4,3)

end subroutine apply_sphamil

subroutine inilag(nx,ny,nz,dx)
  !---------------------------------------------------------------------------
  ! Computes the Lagrange derivative coefficients for this particular
  ! symmetry combination.
  !---------------------------------------------------------------------------

  integer, intent(in) :: nx,ny,nz
  real*8, intent(in)  :: dx

  integer :: i,j, linX, linY, linZ
  real*8  :: sinA, A, B, sinB, C, D

  real*8, parameter  :: pi=4.0d0*atan2(1.0d0,1.0d0)


  ! Allocate the arrays
  allocate(derX(nx,nx,4), derY(ny,nz,4),derZ(nz,nz,4))
  derX   = 0.0d0 ; derY   = 0.0d0 ; derZ   = 0.0d0

  linX = 2*nx
  linY = 2*ny
  linZ = 2*nz

  do i=1,nx
     do j=1,nx
        A           = (pi * (i - j))/linX
        sinA        = sin(A)
        B           = (pi * (i + j-1))/linX
        sinB        = sin(B)

        C = (-1)**(i-j)  *pi/(linX*dx*sinA)
        D = (-1)**(i+j-1)*pi/(linX*dx*sinB) !check the sign

        if(i.eq.j) C = 0

        derX(i,j,1) = C+D
        derX(i,j,2) = C-D
     enddo
  enddo

  do i=1,ny
     do j=1,ny

        A           = (pi * (i - j))/linY
        sinA        = sin(A)
        B           = (pi * (i + j-1))/linY
        sinB        = sin(B)

        C = (-1)**(i-j)       *pi/(linY*dx*sinA)
        D = (-1)**(i+j-1)*pi/(linY*dx*sinB)

        if(i.eq.j) C = 0

        derY(i,j,1) = C+D
        derY(i,j,2) = C-D
     enddo
  enddo

  do i=1,nz
     do j=1,nz
        A           = (pi * (i - j))/linZ
        sinA        = sin(A)
        B           = (pi * (i + j-1))/linZ
        sinB        = sin(B)

        C = (-1)**(i-j)       *pi/(linZ*dx*sinA)
        D = (-1)**(i+j-1)*pi/(linZ*dx*sinB)


        if(i.eq.j) C = 0

        derZ(i,j,1) = C+D
        derZ(i,j,2) = C-D
     enddo
  enddo

end subroutine inilag

subroutine Derive_X(f, px, fx, nx, ny, nz)
    !---------------------------------------------------------------------------
    ! Compute the X-derivative of four 1D-functions on the mesh.
    !
    ! fx = First order derivative in the x direction
    ! px = sign of the symmetry transformation in the x-direction
    !---------------------------------------------------------------------------
    real*8, intent(in)  :: f(nx*ny*nz,4)
    real*8, intent(out) :: fx(nx*ny*nz,4)
    integer, intent(in) :: nx,ny,nz
    integer, intent(in) :: px(4)
    integer             :: m

    do m=1,4
     call Derive_X_single_1D(f(:,m), px(m), fx(:,m),nx,ny,nz)
    enddo

end subroutine Derive_X

subroutine Derive_X_single_1D(f, px, fx, nx, ny, nz)
    !---------------------------------------------------------------------------
    ! Compute the X-derivative of a single real 1D-function on the mesh.
    !
    ! fx = First order derivative in the x direction
    ! px = sign of the symmetry transformation in the x-direction
    !---------------------------------------------------------------------------
    real*8, intent(in), target  :: f(nx*ny*nz)
    real*8, intent(out),target  :: fx(nx*ny*nz)
    integer, intent(in)         :: px
    integer, intent(in)         :: nx,ny,nz
    real*8, pointer             :: f3(:,:,:), fx3(:,:,:)

    f3 (1:nx,1:ny,1:nz) => f (1:nx*ny*nz)
    fx3(1:nx,1:ny,1:nz) => fx(1:nx*ny*nz)

    call Derive_X_single_3D(f3, px, fx3, nx,ny,nz)
 end subroutine Derive_X_single_1D

 subroutine Derive_X_single_3D(f, px, fx, nx, ny, nz)
    !---------------------------------------------------------------------------
    ! Compute the X-derivative of a single real function on the mesh.
    !
    ! fx = First order derivative in the x direction
    ! px = sign of the symmetry transformation in the x-direction
    !---------------------------------------------------------------------------
    real*8, intent(in)  :: f(nx,ny,nz)
    integer, intent(in) :: nx,ny,nz
    real*8, intent(out) :: fx(nx,ny,nz)
    integer, intent(in) :: px
    integer             :: i,j,k,l,sx

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

 end subroutine Derive_X_single_3D

 subroutine Derive_Y(f, py, fy, nx, ny, nz)
    !---------------------------------------------------------------------------
    ! Compute the Y-derivative of four 1D-functions on the mesh.
    !
    ! fy = First order derivative in the y direction
    ! py = sign of the symmetry transformation in the y-direction
    !---------------------------------------------------------------------------
    real*8, intent(in), target  :: f(nx*ny*nz,4)
    real*8, intent(out), target :: fy(nx*ny*nz,4)
    integer,intent(in)          :: nx,ny,nz
    integer, intent(in)                :: py(4)
    integer                            :: m

    do m=1,4
     call Derive_Y_single_1D(f(:,m), py(m), fy(:,m), nx,ny,nz)
    enddo

 end subroutine Derive_Y

 subroutine Derive_Y_single_3D(f, py, fy, nx, ny, nz)
    !---------------------------------------------------------------------------
    ! Compute the Y-derivative of a single 3D-function on the mesh.
    !
    ! fy = First order derivative in the y direction
    ! py = sign of the symmetry transformation in the y-direction
    !---------------------------------------------------------------------------

    real*8, intent(in), target  :: f(nx,ny,nz)
    real*8, intent(out), target :: fy(nx,ny,nz)
    integer, intent(in)                :: py, nx, ny, nz
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

 end subroutine Derive_Y_single_3D

 subroutine Derive_Y_single_1D(f, py, fy, nx,ny,nz)
    !---------------------------------------------------------------------------
    ! Compute the Y-derivative of a single 1D-function on the mesh.
    !
    ! fy = First order derivative in the y direction
    ! py = sign of the symmetry transformation in the y-direction
    !---------------------------------------------------------------------------
    real*8, intent(in), target  :: f(nx*ny*nz)
    real*8, intent(out),target  :: fy(nx*ny*nz)
    integer, intent(in)         :: nx,ny,nz
    integer, intent(in)                :: py
    real*8, pointer             :: f3(:,:,:), fy3(:,:,:)

    f3 (1:nx,1:ny,1:nz) => f (1:nx*ny*nz)
    fy3(1:nx,1:ny,1:nz) => fy(1:nx*ny*nz)

    call Derive_y_single_3D(f3, py, fy3, nx,ny,nz)
 end subroutine Derive_Y_single_1D

   subroutine Derive_Z_single_3D(f, pz, fz, nx, ny, nz)
    !---------------------------------------------------------------------------
    ! Compute the Z-derivative of a single 1D-function on the mesh.
    ! fz = First order derivative in the z direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    real*8, intent(in)  :: f(nx,ny,nz)
    real*8, intent(out) :: fz(nx,ny,nz)
    integer, intent(in) :: nx,ny,nz
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
 end subroutine Derive_Z_single_3D

 subroutine Derive_Z_single_1D(f, pz, fz, nx,ny,nz)
    !---------------------------------------------------------------------------
    ! Compute the Z-derivative of a single real 1D-function on the mesh.
    !
    ! fz = First order derivative in the z direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    real*8, intent(in), target  :: f(nx*ny*nz)
    integer, intent(in)         :: nx,ny,nz
    real*8, intent(out),target  :: fz(nx*ny*nz)
    integer, intent(in)         :: pz
    real*8, pointer             :: f3(:,:,:), fz3(:,:,:)

    f3 (1:nx,1:ny,1:nz) => f (1:nx*ny*nz)
    fz3(1:nx,1:ny,1:nz) => fz(1:nx*ny*nz)

    call Derive_z_single_3D(f3, pz, fz3, nx,ny,nz)
 end subroutine Derive_Z_single_1D

  subroutine Derive_Z(f, pz, fz, nx,ny,nz)
    !---------------------------------------------------------------------------
    ! Compute the Z-derivative of four functions functions on the mesh.
    ! fz = First order derivative in the z direction
    ! pz = sign of the symmetry transformation in the z-direction
    !---------------------------------------------------------------------------
    real*8, intent(in)  :: f(nx*ny*nz,4)
    integer, intent(in) :: nx, ny, nz
    real*8, intent(out) :: fz(nx*ny*nz,4)
    integer, intent(in)        :: pz(4)
    integer                    :: m

    do m=1,4
     call Derive_Z_single_1D(f(:,m), pz(m), fz(:,m), nx, ny, nz)
    enddo
 end subroutine Derive_Z

end module sphamil
