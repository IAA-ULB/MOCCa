module momentsofinertia

  use geninfo
  use densities
  use parameterization

  implicit none

  real(KIND=dp) :: Belyaev(3,3), Rigid(3,3)

contains

  subroutine calcrigid()
    !---------------------------------------------------------------------------
    ! Calculate the rgid-rotor moment of inertia of the density.
    ! 
    !   I^rigid_x = m int d^3r rho * (y^2 + z^2)
    !   I^rigid_y = m int d^3r rho * (x^2 + z^2)
    !   I^rigid_z = m int d^3r rho * (x^2 + y^2)
    !
    ! where m is the mass of the nucleon species under consideration. 
    !
    ! When naively calculated by the code, these numbers have units
    !       MeV * fm^2/c^2
    ! We divide by (hbar c)**2 to obtain the final units of
    !       hbar^2 / MeV
    !---------------------------------------------------------------------------
    real(KIND=dp) :: xs(2), ys(2), zs(2)
    real(KIND=dp), pointer :: rho(:,:,:)
    integer       :: it, i,j,k

    Rigid = 0

    xs = 0 ; ys = 0 ; zs = 0

    do it=1,2
      ! Assigning the storage structure in D_I_I a more readable form
      rho(1:nx, 1:ny, 1:nz) => D_I_I(1:nx*ny*nz,it)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            xs(it) = xs(it) + meshx(i)**2 * rho(i,j,k)
            ys(it) = ys(it) + meshy(j)**2 * rho(i,j,k)
            zs(it) = zs(it) + meshz(k)**2 * rho(i,j,k)
          enddo
        enddo
      enddo
    enddo

    xs = xs * dv ; ys = ys *dv ; zs = zs * dv
      
    Rigid(1,1:2) = nucleonmass * (ys + zs)
    Rigid(2,1:2) = nucleonmass * (xs + zs)
    Rigid(3,1:2) = nucleonmass * (xs + ys)

    Rigid(:,3) = sum(Rigid(:,1:2),2)

    ! Converting to the correct units
    Rigid = Rigid/(hbarclum**2)

  end subroutine calcrigid

  subroutine calcBelyaev()
    
  end subroutine calcBelyaev

  subroutine PrintMomentsofIntertia()
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    1 format (21('-'), 'Rotational properties', 20('-'))
    2 format ('                Rigid rotor (hbar^2/MeV)')
    3 format ('              n', 14x ,'p',14x,'t ')
    4 format (' I_R X ', 3f15.7)
    5 format (' I_R Y ', 3f15.7)
    6 format (' I_R Z ', 3f15.7)
    7 format ('                Belyaev     (hbar^2/MeV) ')
    8 format (' I_B X ', 3f15.7)
    9 format (' I_B Y ', 3f15.7)
   10 format (' I_B Z ', 3f15.7)
  100 format (60('-'))

    call calcrigid()
    call calcBelyaev    

    print 1
    print 2
    print 3
    print 4, Rigid(1,:)
    print 5, Rigid(2,:)
    print 6, Rigid(3,:)
    print 7
    print 3
    print  8, Belyaev(1,:)
    print  9, Belyaev(2,:)
    print 10, Belyaev(3,:)
    print 100
  
  end subroutine PrintMomentsofIntertia
end module momentsofinertia
