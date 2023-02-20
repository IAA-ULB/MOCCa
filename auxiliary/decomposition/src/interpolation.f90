module interpolation 

  use compilation 

  implicit none
  
contains

  subroutine quadrature_points(a,b,n,weight, absc)
    !---------------------------------------------------------------------------
    !
    ! Get the abscissa's and weights for a Gauss-Legendre integration with 
    !  n points. 
    !
    ! Shamelessly stolen from PROMESSE.
    !
    !---------------------------------------------------------------------------
    real*16 , parameter :: zero=0.0q0,one=1.0q0,two=2.0q0,four=4.0q0
    real*16 , parameter :: eps=3.0q-19,half=0.5q0,quart=0.25q0
    real*16             :: xm,xr,pi,z,z1,p1,p2,p3,pp
    real*16             :: ab(1:48),ge(1:48)

    real(kind=dp), intent(in) :: a,b
    integer, intent(in)       :: n
    real(kind=dp), allocatable, intent(out) :: weight(:), absc(:)
    integer       :: i,j
    allocate(absc(n), weight(n))
    ! 
    
    if (mod(n,2) .ne. 0) then
      print *, 'Quadrature points asked for need to be even.'
      stop
    endif
    
    pi = four * atan2(one,one)
    do i=1,n/2
        z=cos(pi*(i-quart)/(n+half))
    1   continue
        p1=one
        p2=zero
        do j=1,n
          p3=p2
          p2=p1
          p1=((two*j-one)*z*p2-(j-one)*p3)/j
        enddo
        pp=n*(z*p1-p2)/(z*z-one)
        z1=z
        z=z1-p1/pp
        if(abs(z-z1).gt.eps) goto 1
        ab(i)=-z
        ge(i)=two/((one-z*z)*pp*pp)
    enddo
    
    !                      scale the abszissas and weights to the desired interval
    !                      (a,b) and compute the symmetric counterparts
    xm=(b+a)*half
    xr=(b-a)*half
    do j=1,n/2,1
      absc  (j)    = xm+xr*ab(j      )
      absc  (n/2+j)= xm-xr*ab(n/2+1-j)
      weight(j)    = xr*ge(j      )
      weight(n/2+j)= xr*ge(n/2+1-j)
    enddo
  end subroutine quadrature_points
  
  function Lagrange_interpolation_function(x, xr, nx, dx) result(f)
  !----------------------------------------------------------------------------- 
  !    Calculate the value of a (one-D) Lagrange interpolation function 
  !    (associated with a reference point x_r) at a (different) point x:
  !                 
  !    f_r(x) = (2 nx)^{-1} sin[ a * (x - x_r) ] / [  sin (a * (x - x_r)/(2nx)]
  !
  !      [Eq. (16)  in W. Ryssens et al., PRC 064318 (2015) ]
  !    
  !    where 
  !        a   = pi/dx
  !        x_r = coordinate of the reference point
  !        dx  = spacing of the mesh
  !        nx  = number of points in the given direction 
  !              (in a symmetry-unrestricted box)
  !-----------------------------------------------------------------------------
  integer, intent(in)       :: nx
  real(KIND=dp), intent(in) :: dx, xr, x

  real(KIND=dp) :: a, f, N

  a = pi/dx
  N = 1.0/(2*nx)

  if(abs(x - xr) .lt. 1e-15 ) then
    f = 1.0d0
  else
    f  = N * sin(a*( x - xr))/(sin(a *N* (x - xr)))
  endif
  end function Lagrange_interpolation_function 
  
  
  function Interpolate(newpoints,fvalues, nx, ny, nz, dx) result (interp)
    !---------------------------------------------------------------------------
    ! Interpolate a function defined on a Cartesian mesh with characteristics
    !   (nx,ny,nz,dx) at a set of newpoints.
    !
    ! this routine assumes EV8 symmetries and that the function being 
    ! interpolated does not incur a sign under any plane reflection symmetry.
    !
    ! In:
    !   newpoints   : Cartesian coordinates (x,y,z) as a function of two angular
    !                 coordinates?
    !   fvalues     : values of the function being interpolated
    !   nx,ny,nz,dx : characteristics of the Lagrange mesh
    !
    ! Out:
    !   interp      : values of the function at the new points
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: newpoints(:,:,:), fvalues(nx,ny,nz,3)
    integer, intent(in)       :: nx,ny,nz
    real(KIND=dp), intent(in) :: dx
    
    real(KIND=dp) :: Lx, Ly, Lz, xi,yi,zi
    real(KIND=dp), allocatable :: interp(:,:,:)
    
    integer :: i,j,k,l, N, M, q
    
     N = size(newpoints(:,1,1))
     M = size(newpoints(1,:,1))
     
     allocate(interp(N,M,3)) ; interp = 0.0d0
  
     do l=1,N
      do q=1,M
       do k =1, nz
          zi= dx/2 + (k-1) * dx      

          Lz =      Lagrange_interpolation_function(newpoints(l,q,3),  zi, 2*nz, dx)
          Lz = Lz + Lagrange_interpolation_function(newpoints(l,q,3), -zi, 2*nz, dx)

          do j=1,ny
            yi= dx/2 + (j-1) * dx      

            Ly =      Lagrange_interpolation_function(newpoints(l,q,2),  yi, 2*ny, dx) 
            Ly = Ly + Lagrange_interpolation_function(newpoints(l,q,2), -yi, 2*ny, dx)

            do i=1,nx
               xi= dx/2 + (i-1) * dx      
      
               Lx =      Lagrange_interpolation_function(newpoints(l,q,1),  xi, 2*nx, dx) 
               Lx = Lx + Lagrange_interpolation_function(newpoints(l,q,1), -xi, 2*nx, dx)

               interp(l,q,:) = interp(l,q,:) + Lx * Ly * Lz * fvalues(i,j,k,:)
            enddo
          enddo
       enddo  
      enddo 
     enddo 
  
  end function Interpolate 
  
  
end module interpolation
