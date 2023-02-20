program decomposition
  !-----------------------------------------------------------------------------
  ! Decomposition of a nuclear density defined on a Lagrange mesh into 
  ! spherical harmonics components, see the notes in this directory. 
  !                                                        Copyright: W. Ryssens
  !-----------------------------------------------------------------------------

  use compilation 
  use sphericalharmonics
  use interpolation

  implicit none

  1 format ('------------------ Density decomposition ------------------------')
  2 format ( ' Mesh parameters' )
  3 format ( '   nx = ', i5 , ' ny = ' , i5 , ' nz = ' , i5, ' mv = ' , i5)
  4 format ( '   dx = ', f20.10,' (fm  ) ')
  5 format ( '   dv = ', f5.2,' (fm^3) ')
 51 format ( ' Quantisation axis  = ', i1)  
  6 format ('  Inputfile  :  ', a100)
  7 format ('  Outputfile :  ', a100)
  8 format ('  Skipping ', i3, ' rows at the start of the density file. ')
 81 format ('  Reading  ', i3, ' density columns from the input file. ')
  9 format ('-----------------------------------------------------------------')
 10 format (' Multipole moments up to l = ', i3)
 11 format (' Radial mesh' )
 12 format ('   rmin = ', f7.3, ' fm')
 13 format ('   rmax = ', f7.3, ' fm')
 14 format ('   dr   = ', f7.3, ' fm')
 15 format ('   nr   = ',   i3)
 16 format (' Angular mesh' )
 17 format ('   ntheta = ',  i3)
 18 format ('   nphi   = ',  i3)
 
 20 format ('---------------  Verification -----------------------------------')
 21 format ('                       Original       Decomposition     Ratio')
 22 format (' Number of particles:', 3es15.7)
 23 format (' RMS radius:         ', 3es15.7)
 24 format (' Q_',2i2' :            ', 3es15.7)

  !-----------------------------------------------------------------------------
  ! Parameters of the original Lagrange mesh
  integer       :: nx,ny,nz, mv
  real(KIND=dp) :: dx, dv
  !-----------------------------------------------------------------------------
  ! Maximum ell to do the decomposition for
  integer       :: maxl = 4
  !-----------------------------------------------------------------------------
  ! Filenames for input and output
  character(len=100)  :: inputfilename, outputfilename
  logical             :: exists
  !-----------------------------------------------------------------------------
  ! Input in Tantalus format contains a few lines that this code should not 
  ! read. Should be adapted if there are no comments at the top of a .den file
  integer :: skiprows = 14
  !-----------------------------------------------------------------------------
  ! Number of densities that should be read from the original file
  integer :: ncol = 3
  !-----------------------------------------------------------------------------
  ! Quantisation axis for the orientation of the multipole moments:
  !   1 => X-axis
  !   2 => Y-axis
  !   3 => Z-axis
  ! See the tantalus documentation for more information.
  integer :: quantisationaxis = 3
  !-----------------------------------------------------------------------------
  ! Values of the spherical harmonics on our (theta,phi) mesh
  !                                        |-phi 
  !                                        v v- theta
  real(KIND=dp), allocatable :: spherharm(:,:,:,:)  
  !                                            ^ ^
  !                                            | - m
  !                                            - - l
  ! Other angle-dependent integration factors
  real(KIND=dp), allocatable :: integ_factors(:,:)  !(phi, theta) 
  ! Values of the spherical harmonics on the Cartesian mesh, for error checking
  real(KIND=dp), allocatable :: sph_old(:,:,:,:,:,:)
  !-----------------------------------------------------------------------------
  ! Density on the original Lagrange mesh
  !                                    nx ny nz
  real(KIND=dp), allocatable :: density(:,:,:,:)
  !                                           | isospin: 1 neutrons
  !                                                      2 protons
  !                                                      3 charge
  !-----------------------------------------------------------------------------
  ! Coordinates of the coordinate mesh, only used for verification purposes
  real(KIND=dp), allocatable :: x(:), y(:), z(:)
  !-----------------------------------------------------------------------------
  ! Values of the decomposition of the density
  !                                        | r
  !                                        v v- l
  real(KIND=dp), allocatable :: decomp_den(:,:,:,:)
  !                                            ^ ^- isospin (as above)
  !                                            | - m 
  !
  !-----------------------------------------------------------------------------
  ! Coordinates of the new mesh in spherical coordinates and interpolated
  real(KIND=dp), allocatable :: r(:), theta(:), phi(:)
  ! Cartesian coordinates of this new mesh
  real(KIND=dp), allocatable :: newpoints(:,:,:)
  ! Interpolated density on this new mesh
  real(KIND=dp), allocatable :: interden(:,:,:)
  ! Parameters of the radial part of the mesh
  real(KIND=dp) :: rmin=0.1, rmax=10.0, dr 
  ! Number of points in each discretisation direction
  integer       :: nr = 20 ,  ntheta = 6, nphi = 6  
  ! Integration weights (to be determined for theta and phi)  
  real(KIND=dp), allocatable :: th_weight(:), phi_weight(:)
  !-----------------------------------------------------------------------------
  ! Looping and temporary variables
  integer :: i,j,k,l,m, q, rr
  real(KIND=dp) :: costheta,cosmphi, factorialquotient, fac, sintheta, rms_old
  real(KIND=dp) :: Anew, Aold, oldQ, newQ
  !-----------------------------------------------------------------------------
  
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Read information from STDIN
  NameList /main/ nx,ny,nz,dx, inputfilename, outputfilename, skiprows, &
  &               rmin, rmax, nr, maxl, ntheta, nphi, ncol, quantisationaxis
  read (unit=*, nml=main)
 
  mv = nx*ny*nz
  dv = 8*dx**3
  dr = (rmax-rmin)/(nr-1)
 
  print 1
  print 2
  print 3, nx,ny,nz,mv
  print 4, dx
  print 5, dv
  print 51, quantisationaxis
  print 6, Inputfilename
  print 7, Outputfilename
  print 8, skiprows
  print 81, ncol
  print 9
  print 10, maxl
  print 11
  print 12, rmin
  print 13, rmax
  print 14, dr
  print 15, nr
  print 16
  print 17, ntheta
  print 18, nphi
  print 9

  allocate(x(nx), y(ny), z(nz))
  allocate(density(nx,ny,nz,ncol))

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Read a density from file in the format of Tantalus, meaning
  !
  !     X   Y   Z    Rho_n   Rho_p  Rho_charge
  !
  ! X,Y,Z are of no relevance and are discarded.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   
  inquire(file=inputfilename, exist=exists)
  if(.not.exists) then
    print *, 'Input file specified does not exist!'
    stop
  endif
  open (unit=10,file=inputfilename)
  
  do k=1,skiprows
    read(unit=10,fmt=*)      
  enddo
  
  do k=1,nz
    do j=1,ny
      do i=1,nx
        read(unit=10,fmt=*) x(i), y(j), z(k), density(i,j,k,:)
      enddo
    enddo
  enddo
  close(unit=10)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Construct the interpolation mesh in terms of 
  ! * r     :  simple linearly spaced points
  ! * theta :  Gauss-Legendre abscissa's 
  ! * phi   :  Gauss-Legendre abscissa's
  allocate(r(nr))
  do i=1,nr
    r(i) = rmin + (i-1) * dr 
  enddo  

  ! We only integrate theta up to pi/2 and multiply by two afterwards
  ! which we can do because of EV8 symmetry!
  call quadrature_points(0.0d0, pi/2,ntheta,th_weight,theta)
  th_weight = 2 *th_weight

  ! We only integrate phi up to pi/2 and multiply by two afterwards
  ! which we can do because of EV8 symmetry and we only integrate even m
  call quadrature_points(0.0d0, pi/2,nphi,phi_weight,phi)
  phi_weight = 4 *phi_weight

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! Generate the relevant spherical harmonic values
  ! Note: 
  !  (a)  we only need the real components ov Ylm with l and m both even
  !       but calculate all real parts for human time efficiency
  !
  !  (b) we multiply below by the factor
  !                  (-1)^{m}[(2*l+1)/4\pi (l-m)!/(l+m)!]^{1/2}
  !                
  
  allocate(spherharm(nphi,ntheta,0:maxl, 0:maxl)) ; spherharm = 0.0d0
  allocate(integ_factors(nphi,ntheta))            ; integ_factors = 0.0d0

  do j=1,ntheta
    costheta = cos(theta(j))
    do m=0,maxl
      do i=1, nphi    
        cosmphi = cos(m*phi(i))
        ! Ylm
        spherharm(i,j,:,m)     = legendre_pnm(maxl, m, cosTheta) * cosmphi
      enddo
    enddo
  enddo
  ! Calculation of position-independent factor
  do m=0, maxl
    do l=m,maxl ! Note that l>= m
      factorialquotient=1.0_dp
      do q=1,2*m
        factorialquotient = factorialquotient/( l + m - q + 1)
      enddo
      fac=(-1)**m*sqrt(factorialquotient)
      fac = fac * sqrt((2*l +1)/(4 * pi))
      spherharm(:,:,l,m) = spherharm(:,:,l,m) * fac 
    enddo
  enddo
  do j=1,ntheta
    sintheta = sin(theta(j))
    do i=1,nphi
      ! All other angle-dependent factors in the integral
      integ_factors(i,j) = sinTheta * th_weight(j) * phi_weight(i) 
    enddo
  enddo

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Now we start actual computation
  allocate(decomp_den(nr,0:maxl, 0:maxl,ncol), interden(nphi, ntheta,ncol))
  allocate(newpoints(nphi, ntheta,3)) 
  decomp_den = 0.0d0 ; interden = 0.0d0 ; newpoints = 0.0d0
  
  ! Loop over radial coordinate
  do rr=1,nr
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! First, we build all the points that we are going to interpolate to
    do j=1,ntheta
      do i=1,nphi
        select case(quantisationaxis)
        case(1)
          newpoints(i,j,2) = r(rr) * sin(theta(j)) * cos(phi(i))   ! y = r sin(theta) cos(phi)
          newpoints(i,j,3) = r(rr) * sin(theta(j)) * sin(phi(i))   ! z = r sin(theta) sin(phi)
          newpoints(i,j,1) = r(rr) * cos(theta(j))                 ! x = r cos(theta)
        case(2)
          newpoints(i,j,1) = r(rr) * sin(theta(j)) * cos(phi(i))   ! x = r sin(theta) cos(phi)
          newpoints(i,j,3) = r(rr) * sin(theta(j)) * sin(phi(i))   ! z = r sin(theta) sin(phi)
          newpoints(i,j,2) = r(rr) * cos(theta(j))                 ! y = r cos(theta)
        case(3)
          newpoints(i,j,1) = r(rr) * sin(theta(j)) * cos(phi(i))   ! x = r sin(theta) cos(phi)
          newpoints(i,j,2) = r(rr) * sin(theta(j)) * sin(phi(i))   ! y = r sin(theta) sin(phi)
          newpoints(i,j,3) = r(rr) * cos(theta(j))                 ! Z = r cos(theta)
        end select
      enddo
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! We build the interpolated density
    interden = Interpolate(newpoints,density, nx, ny, nz, dx)
    !    
    ! Looping over all EVEN (and real) multipole moments
    do l=0, maxl,2
      do m=0, l, 2
        do q=1,ncol
          ! The spherharm array takes care of all relevant factors in the 
          ! integrand that are not the density.
          decomp_den(rr, l, m, q) = sum(interden(:,:,q) * spherharm(:,:,l,m) & 
          &                         *   integ_factors(:,:)) 

          if(m .ne. 0) then
              decomp_den(rr,l,m,q) =  2*decomp_den(rr,l,m,q)           
          endif
        enddo
      enddo
    enddo
  enddo
 
  !-----------------------------------------------------------------------------
  ! Start output
  !-----------------------------------------------------------------------------
  open (unit=10,file=outputfilename)

  write(unit=10, fmt='(a10,7x)', advance='no') '#R [fm]'
  do l = 0, maxl,2
    do m = 0, l, 2
      do q=1,ncol
        write(unit=10, fmt='(a5,2i1,a8)', advance='no') ' rho_', l, m, ' [fm^-3]'
      enddo
    enddo
  enddo
  write(unit=10, fmt='()')


  do rr=1,nr
    write(unit=10, fmt='(e15.6)', advance='no') r(rr)  
    do l = 0, maxl,2
      do m = 0, l, 2
        do q=1,ncol
          write(unit=10, fmt='(e15.6)', advance='no') decomp_den(rr,l,m,q)
        enddo
      enddo
    enddo
    write(unit=10, fmt='()')      
  enddo
  close(unit=10)

  !-----------------------------------------------------------------------------
  ! Start verification
  !-----------------------------------------------------------------------------

  ! Generate spherical harmonics on the old mesh
  allocate(sph_old(nx,ny,nz,0:maxl,0:maxl,2))
  call GenSphericalHarmonics(maxl,nx,ny,nz,x,y,z,sph_old,quantisationaxis,1)
 

  print 20
  print 21
  do q=1,ncol
    print *, 'Density ', q
    Aold = sum(density(:,:,:,q))*dv
    Anew = sum(decomp_den(:,0,0,q)*r**2)*dr*sqrt(4*pi)
    print 22, Aold, Anew, Anew/Aold
    
    if(maxl .ge. 2) then
      rms_old = 0

      do k=1,nz
        do j=1,ny
          do i=1,nx
            rms_old = rms_old + (x(i)**2 + y(j)**2 + z(k)**2) * density(i,j,k,q)
          enddo
        enddo
      enddo
      rms_old = rms_old *dv/Aold
    endif  
    print 23, sqrt(rms_old), &
    &         sqrt(sum(decomp_den(:,0,0,q)*r**4)*dr*sqrt(4*pi)/Anew), &
    &         sqrt(sum(decomp_den(:,0,0,q)*r**4)*dr*sqrt(4*pi)/Anew)/sqrt(rms_old)

    do l=0, maxl,2
      do m=0,l,2
        oldQ = sum(sph_old(:,:,:,l,m,1) * density(:,:,:,q))*dv
        newQ = sum(r**(l+2) * decomp_den(:,l,m,q))*dr
        
        ! non-trivial factor two because of our conventions!
        if(m.ne.0) newQ = newQ/2
        print 24, l,m,oldQ,newQ, newQ/oldQ
      enddo
    enddo
    print *
  enddo
  
end program decomposition
