module cranking
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
 ! Module for all things cranking and angular momentum.
 !
 !------------------------------------------------------------------------------
 ! HEPHAESTOS keywords:
 ! 
 ! CRANKDIR : $CRANKDIR
 ! CRANKLEN : $CRANKLEN
 !    The array crankdirections controls the axes along which we calculate 
 !    angular momentum. It has size $CRANKLEN+1 and has entries $CRANKDIR plus 
 !    a trailing zero. The trailing zero is included, because the code has to 
 !    compile when we calculate no angular momentum at all, i.e. when time-
 !    reversal is conserved.
 !==============================================================================

 use compilation
 use geninfo
 use derivatives
 use wavefunctions
 use geninfo
 use nil8
 use pairing

 implicit none
 
 !------------------------------------------------------------------------------
 ! Omega:
 !   Cranking frequency or Lagrange multiplier of the angular momentum in 
 !   the three Cartesian directions.
 !------------------------------------------------------------------------------
 real(KIND=dp) :: Omega(3)      = 0.0_dp
 !------------------------------------------------------------------------------
 ! Crankenergy:
 !   Energy associated with the cranking constraint in each Cartesian direction,
 !   i.e. CE(i) = - omega_i * <J_i>
 !------------------------------------------------------------------------------
 real(KIND=dp) :: crankenergy(3)   = 0.0_dp
 !------------------------------------------------------------------------------
 ! TotalAngMom:
 !    Total angular momentum in the three Cartesian directions, calculated
 !    by summation of the single-particle contributions.
 ! AngMomOld: 
 !    Values of the total angular momentum at the previous iteration, used for
 !    readjustment of the cranking constraints.
 ! J2:
 !    Values of the total angular momentum squared, <J_i^2>, for the three
 !    Cartesian directions. 
 ! AMBlock:
 !    Values of the total angular momentum, split by quantum number block.
 !------------------------------------------------------------------------------
 real(KIND=dp), public :: TotalAngMom(3)= 0.0_dp, AngMomOld(3)  = 0.0_dp
 real(KIND=dp), public :: J2(3)         = 0.0_dp, AMBlock(8,3)  = 0.0_dp

 !------------------------------------------------------------------------------
 integer, parameter                        :: cranklen = $CRANKLEN
 integer, parameter, dimension(cranklen+1) :: crankdirections = (/ $CRANKDIR 0/)

contains

  subroutine readcranking(file_number)
    !---------------------------------------------------------------------------
    ! Read the namelist &cranking/
    !---------------------------------------------------------------------------
    integer(dp), intent(in),optional :: file_number
    real(KIND=dp)       :: OmegaX, OmegaY, OmegaZ
    integer             :: i,j,c
    logical             :: NotFound

    namelist /cranking/ OmegaX, OmegaY, OmegaZ

    OmegaX = 0 ; OmegaY = 0 ; OmegaZ = 0

    if(present(file_number)) then
      read (unit=file_number, nml=cranking)
    else
      read (unit=*, nml=cranking)
    endif    

    Omega = (/ OmegaX, OmegaY, OmegaZ /)
    
    ! Check if the asked for value of the vector omega is allowed by the 
    ! current symmetries
    do i=1,3
      if(Omega(i) .ne. 0.0d0) then
        NotFound = .true.
        do j=1,cranklen
          c = crankdirections(j)
          if( c == i ) NotFound = .false.
        enddo
        if(NotFound) then
          print *, 'Disallowed cranking frequency.'
          stop
        endif
      endif
    enddo
  end subroutine readcranking

  subroutine updateAM
    !---------------------------------------------------------------------------
    ! Calculate the total angular momentum and related observables.
    !---------------------------------------------------------------------------  
    integer :: B, N, wave, si, i, c

    totalangmom = 0.0

    si = 0    
    do B=1,8
      N = HFBlocks(B) ; if(N .eq. 0) cycle
      do wave = 1, N  
        do i = 1, cranklen
          c  = crankdirections(i)
          TotalAngMom(c) = TotalAngMom(c) + rho_can(si+wave) * spwf_J(c,si+wave)
        enddo
      enddo
      si = si + N
    enddo

    crankenergy = - omega * TotalAngMom    
  end subroutine updateAM

  subroutine printcranking_init
    !---------------------------------------------------------------------------
    ! Print information on the cranking constraints imposed on the calculation
    ! at the start of the program.
    !---------------------------------------------------------------------------
    1 format ( 29('-'), ' Cranking information ', 29('-'))
    2 format ( ' Cranking frequencies ', /, & 
    &          '    Omega_X = ', f15.3, /,  &
    &          '    Omega_Y = ', f15.3, /,  &
    &          '    Omega_Z = ', f15.3) 

    print 1
    print 2, Omega
        

  end subroutine printcranking_init

  subroutine PrintCranking
    !---------------------------------------------------------------------------
    ! Prints all kinds of information about the expectation value of the 
    ! angular momentum operator and all kinds of angles.
    !---------------------------------------------------------------------------

    1 format (2x,74('_') )
   10 format (2x,74('-'))
    2 format (25('-'), ' Angular Momentum (hbar) ',26('-') )
    3 format (15x, 'Spwfs(*)  ',2x, 'Desired', 5x, 'Omega', 7x, 'E (MeV)' 6x,'Densit. ')

    4 format (3x,'J_',a1,'   ','|', 5f12.5 )
   31 format (3x,'Size  |', 3f12.5,12x,1f12.5)
!   32 format (1x,'ReJT',a1,'   ','|', 5f12.5 )
!   33 format (1x,'ImJT',a1,'   ','|', 5f12.5 )
!   34 format (2x,'|J|' ,a1,'   ','|', 5f12.5 )

    print 2
    print *
    print 3
    print 1

    print 4, 'x', TotalAngMom(1), 0.0, Omega(1), CrankEnergy(1), 0.0 
    print 4, 'y', TotalAngMom(2), 0.0, Omega(2), CrankEnergy(2), 0.0 
    print 4, 'z', TotalAngMom(3), 0.0, Omega(3), CrankEnergy(3), 0.0
    print 1
    print 31, sqrt(sum(totalangmom(1:3)**2)), 0.0, &
    &         sqrt(sum(omega(1:3)**2))      , 0.0
    print 10
  end subroutine PrintCranking

  function crank_spin_potential() result(spot)
    !---------------------------------------------------------------------------
    ! The cranking constraint contributes to the field F_I_S, associated with 
    ! the spin density D_I_S:
    !
    !     F_I_S(r) => F_I_S(r) - \frac{1}{2} \omega f_cut(r)
    !
    ! where f_cut(r) is the cutoff factor also employed for the multipole 
    ! constraints. The factor 1/2 is present because the J_spin = 1/2 Pauli 
    ! sigma matrix.
    !---------------------------------------------------------------------------
    use Moments, only : cutoff

    real(KIND=dp), allocatable :: spot(:,:,:)
    integer :: i, it, c

    allocate(spot(nx*ny*nz,3,2)) ; spot = 0.0d0

    do i=1, cranklen
      c           = crankdirections(i)
      do it=1,2
        spot(:,c,it) = - 0.5 * omega(c) * cutoff(:,it)
      enddo
    enddo
      
    return
  end function crank_spin_potential

  function crank_current_potential() result(jpot)
    !---------------------------------------------------------------------------
    ! The cranking constraint contributes to the field G_I_N, associated with 
    ! the current density D_I_N:
    !
    !       G_I_N => G_I_N - f_cut(r) \vec{\omega} x \vec{r} 
    !---------------------------------------------------------------------------
    use Moments, only : cutoff

    real(KIND=dp), allocatable :: jpot(:,:,:)
    integer                    :: it, mu, indices(2) , nu, ka

    allocate(jpot(nx*ny*nz,3,2)) ; jpot = 0.0d0
    do mu=1, 3
      indices = vector_product(mu)
      nu = indices(1)
      ka = indices(2)
      do it=1,2
          jpot(:,mu, it) = - cutoff(:,it) * &
          &            (omega(nu) * meshgrid(:,ka) - omega(ka) * meshgrid(:,nu))
      enddo
    enddo
    return
  end function crank_current_potential

end module 
