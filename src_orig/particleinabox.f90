module particleinbox
  !-----------------------------------------------------------------------------
  ! The volume element is 
  !     dv = 2* dx
  !
  !
  !-----------------------------------------------------------------------------


  use compilation
  use geninfo
  use derivatives
  use parameterization

  implicit none

  integer, parameter:: boxstates = 10

contains
  
    subroutine ortho_1D(psi)
        !-----------------------------------------------------------------------
        ! Orthonormalizes the 1D states.
        !-----------------------------------------------------------------------

        integer       :: i, j, N,p
        real(KIND=dp) :: psi(:,:,:) , ME

        N = size(psi,2)
        do p=1,2
          do i=1,N
            ! Normalize state i 
            ME =  sum(psi(:,i,p)**2) * 2 * dx
            psi(:,i,p) = psi(:,i,p)/sqrt(ME)
            do j=i+1,N
              ! Orthogonalize against all the others (of the same parity!)
              ME =  sum(psi(:,i,p)* psi(:,j,p)) * 2 * dx
              psi(:,j,p) = psi(:,j,p) - ME * psi(:,i,p)
            enddo
          enddo
        enddo

    end subroutine ortho_1D

    subroutine solve_pib
      !-------------------------------------------------------------------------
      ! Solve the one-dimensional particle in a box problem for this mesh. 
      !
      ! The hamiltonian is 
      !              H = - hbar^2/(2m) * Delta
      ! but the hbar^2/(2m) is a simple scaling factor that does not affect the
      ! form/shape of the wavefunctions.  
      !
      ! The solutions are of the form (up to a normalization constant) for
      ! k = 1,2,3,...
      !
      !                    cos( k * pi * Z/(N*dx))    if k is odd
      !    psi_k(x) = 
      !                    sin( k * pi * Z/(N*dx))     if k is even
      !
      ! and they differ on wether they are zero or nonzero at z = 0. 
      !
      ! Note that the Hamiltonian does not involve spin, hence the 
      ! single-particle states can a) be perfectly represented as fully spin 
      ! up (and the spin-down states can be obtained by time-reversal) and b)
      ! be chosen completely real. 
      !
      ! Using these features, the symmetry properties of all the psi_k are 
      ! completely determined by k: the cosines are even under z=>-z and the 
      ! sines are odd.
      !-------------------------------------------------------------------------
      
      real(KIND=dp), allocatable :: psi(:,:,:), df(:)
      real(KIND=dp), allocatable :: boxenergies(:,:), disper(:,:)

      real(KIND=dp) :: stepsize = 0.05, k, x(nz)
      integer       :: iter, s, i,j, p

      1 format (16f7.3)

      allocate(psi(nz, boxstates, 2))
      allocate(boxenergies(boxstates,2)) ; allocate(disper(boxstates,2))
      allocate(df(nz))

      call random_number(psi)
      call ortho_1D(psi)

      do i=1,nz
          x(i) = 0.5 * dx + (i-1) * dx
      enddo

      do iter=1, 400
        do p=1,2
          do i=1,boxstates
            df = - matmul(laplaZ(:,:,p),psi(:,i,p)) !+ 0.05 * x**2 * psi(:,i,p)
            ! Calculate expectation values
            boxenergies(i,p) =   sum(psi(:,i,p) * df)                   * 2 * dx
            disper(i,p)      =   sum(df**2)*2*dx - boxenergies(i,p)**2
            ! Evolve
            psi(:,i,p) = psi(:,i,p) - stepsize*(df-boxenergies(i,p)*psi(:,i,p))
          enddo            
        enddo          
        ! Orthonormalize
        call ortho_1D(psi)      

        print *, iter, sum(disper(:,1))/boxstates, sum(disper(:,2))/boxstates
      enddo
      print *
      do i=1,boxstates
        k = pi * (i-0.5) /(nz*dx)
        print *, i, boxenergies(i,1),boxenergies(i,1)/k**2, disper(i,1)
      enddo
      print *
      do i=1,boxstates
        k = pi * (i-0.5)/(nz*dx)
        print *, i, boxenergies(i,2), boxenergies(i,2)/k**2,disper(i,2)
      enddo

      open(unit=6,file='wfs.-.dat' )
      do i=1, boxstates
        write(unit=6, fmt=1), psi(:,i,1)
      enddo
      close(unit=6)

      open(unit=6,file='wfs.+.dat' )
      do i=1, boxstates
        write(unit=6, fmt=1), psi(:,i,2)
      enddo
      close(unit=6)

      stop
    end subroutine solve_pib



end module particleinbox
