module cranking
 !==============================================================================
 !_________ _______  _       _________ _______  _                 _______
 !\__   __/(  ___  )( (    /|\__   __/(  ___  )( \      |\     /|(  ____ \
 !   ) (   | (   ) ||  \  ( |   ) (   | (   ) || (      | )   ( || (    \/
 !   | |   | (___) ||   \ | |   | |   | (___) || |      | |   | || (_____
 !   | |   |  ___  || (\ \) |   | |   |  ___  || |      | |   | |(_____  )
 !   | |   | (   ) || | \   |   | |   | (   ) || |      | |   | |      ) |
 !   | |   | )   ( || )  \  |   | |   | )   ( || (____/\| (___) |/\____) |
 !   )_(   |/     \||/    )_)   )_(   |/     \|(_______/(_______)\_______)
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
 !
 !  TR  : $TR
 ! NTR  : $NTR
 !==============================================================================

 use compilation
 use geninfo
 use derivatives
 use wavefunctions
 use geninfo
 use nil8
 use pairing
 use densities

 implicit none

 !------------------------------------------------------------------------------
 ! Omega:
 !   Cranking frequency or Lagrange multiplier of the angular momentum in
 !   the three Cartesian directions.
 !------------------------------------------------------------------------------
 real(KIND=dp) :: Omega(3)       = 0.0_dp
 real(KIND=dp) :: Omega_prev(3)  = 0.0_dp
 !------------------------------------------------------------------------------
 ! Crankvalues:  Values of J_mu to use in the cranking constraint
 !------------------------------------------------------------------------------
 real(KIND=dp) :: CrankValues(3) = 0.0_dp
 !------------------------------------------------------------------------------
 ! CrankIntensity:
 !  Intensity of the cranking constraint if cranktype_mu = 1. The code will use
 !  a simple heuristic if left to zero while a constraint is detected.
 !------------------------------------------------------------------------------
 real(KIND=dp) :: CrankIntensity(3) = 0.0_dp
 !------------------------------------------------------------------------------
 ! CrankScalefactor: scaling factor for the projection on the feasible subspace
 !                   see subroutine FeasibleProject in evolution.f90
 !------------------------------------------------------------------------------
 real(KIND=dp) :: CrankScaleFactor(3) = 0.0_dp
 !------------------------------------------------------------------------------
 ! Crankenergy:
 !   Energy associated with the cranking constraint in each Cartesian direction,
 !   i.e. CE(i) = - omega_i * <J_i>.
 !
 ! Crankenergy_cut:
 !   Energy associated with the cranking constraint in each Cartesian direction,
 !   i.e. CE(i) = - omega_i * <J_i>, but where <J_i> is calculated by
 !   integration over the current and spin densities with the multipole cutoff.
 !------------------------------------------------------------------------------
 real(KIND=dp) :: crankenergy(3)     = 0.0_dp
 real(KIND=dp) :: crankenergy_cut(3) = 0.0_dp
 !------------------------------------------------------------------------------
 ! Cranktype
 !    Determines the type of constraint (every cartesian direction separately)
 integer       :: cranktype(3) = 0
 ! Crank_smooth
 ! Whether the update of the cranking constraint uses
 ! (a) .false. => TotalAngMom     , calculated from the spwfs directly
 ! (b) .true.  => TotalAngMom_dens, calculated by integrating the densities
 logical       :: crank_smooth = .false.
 !-----------------------------------------------------------------------------
 ! Whether or not to use the cranking info from file
 logical               :: ContinueCrank= .false.
 !------------------------------------------------------------------------------
 ! TotalAngMom:
 !    Total angular momentum in the three Cartesian directions, calculated
 !    by summation of the single-particle contributions.
 ! TotalAngMom_dens:
 !    Total angular momentum in the three Cartesian directions, calculated
 !    by integration of the spin and current densities.
 ! TotalAngMom_cut:
 !    Total angular momentum in the three Cartesian directions, calculated
 !    by integration of the spin and current densities, but including the
 !    multipole cutoff.
 ! AngMomOld:
 !    Values of the total angular momentum at the previous iteration, used for
 !    readjustment of the cranking constraints.
 ! J2_sp:
 !    Values of the single-particle part of the total angular momentum squared,
 !            <J_i^2>_sp = sum_i rho_ii < i | J_i^2 | i>
 !    where the sum is in the canonical basis and the second J_i on the right
 !    is a SINGLE-PARTICLE operator. This quantity IS NOT EQUAL to <J^2_i>,
 !    the many-body operator.
 !
 ! AMBlock:
 !    Values of the total angular momentum, split by quantum number block.
 !------------------------------------------------------------------------------
 real(KIND=dp) :: TotalAngMom(3)= 0.0_dp, AngMomOld(3)  = 0.0_dp
 real(KIND=dp) :: J2_sp(3)      = 0.0_dp, AMBlock(8,3)  = 0.0_dp
 real(KIND=dp) :: TotalAngMom_dens(3) = 0.0_dp, AngMomOld_dens(3)  = 0.0_dp
 real(KIND=dp) :: TotalAngMom_cut(3)  = 0.0_dp, AngMomOld_cut(3)   = 0.0_dp
 !------------------------------------------------------------------------------
 integer, parameter                        :: cranklen = $CRANKLEN
 integer, parameter, dimension(cranklen+1) :: crankdirections = (/ $CRANKDIR 0/)

contains

  function check_cranking() result(checked)
    !---------------------------------------------------------------------------
    ! Simple function that indicates whether there are active cranking
    ! constraints requiring projection on the feasible subspace.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Output:
    !    checked: .false. if no such constraints present, .true. else.
    !
    !---------------------------------------------------------------------------
    logical :: checked
    integer :: i

    checked = .false.
    do i=1,3
      if(cranktype(i) .eq. 1) then
        checked = .true.
      endif
    enddo

  end function check_cranking

  subroutine readcranking(file_number)
    !---------------------------------------------------------------------------
    ! Read the namelist &cranking/.
    !
    ! Input:
    !     file_number : channel number of opened file where to read from.
    !                   Optional. If not present, read from STDIN.
    !---------------------------------------------------------------------------
    integer(dp), intent(in),optional :: file_number
    real(KIND=dp)       :: OmegaX, OmegaY, OmegaZ
    real(KIND=dp)       :: CrankX, CrankY, CrankZ
    real(KIND=dp)       :: IntensityX, IntensityY, IntensityZ
    real(KIND=dp)       :: ScaleX, ScaleY, ScaleZ
    integer             :: CrankTypeX, CrankTypeY, CrankTypeZ, i
    logical             :: NotFound
$NTR  integer             :: j,c
#if(USE_MPI>0)
    integer             :: mpi_err
#endif

    namelist /cranking/ OmegaX, OmegaY, OmegaZ,             &
    &                   CrankX, CrankY, CrankZ,             &
    &                   CrankTypeX, CrankTypeY, CrankTypeZ, &
    &                   IntensityX, IntensityY, IntensityZ, &
    &                   ScaleX, ScaleY, ScaleZ,             &
    &                   crank_smooth, ContinueCrank

    OmegaX     = 0 ; OmegaY     = 0 ; OmegaZ     = 0
    CrankX     = 0 ; CrankY     = 0 ; CrankZ     = 0
    CrankTypeX = 0 ; CrankTypeY = 0 ; CrankTypeZ = 0
    ScaleX     = 1 ; ScaleY     = 1 ; ScaleZ     = 1
    IntensityX = 0 ; IntensityY = 0 ; IntensityZ = 0

    if(MPI_RANK.eq.0) then
      ! Reading information with the very first MPI rank
      if(present(file_number)) then
        read (unit=file_number, nml=cranking)
      else
        read (unit=*, nml=cranking)
      endif

      !----------------- Assigning Constants based on Input --------------------
      CrankValues      = (/ CrankX,CrankY,CrankZ/)
      Omega            = (/ OmegaX,OmegaY,OmegaZ/)
      CrankType        = (/ CrankTypeX, CrankTypeY, CrankTypeZ/)
      CrankScaleFactor = (/ ScaleX, ScaleY, ScaleZ/)
      CrankIntensity   = (/ IntensityX, IntensityY, IntensityZ/)
    endif

#if(USE_MPI > 0)
    call MPI_Bcast(CrankValues       , 3, MPI_REAL8  ,0,MPI_COMM_WORLD,mpi_err)
    call MPI_Bcast(Omega             , 3, MPI_REAL8  ,0,MPI_COMM_WORLD,mpi_err)
    call MPI_Bcast(CrankType         , 3, MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    call MPI_Bcast(CrankScaleFactor  , 3, MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
    call MPI_Bcast(CrankIntensity    , 3, MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
#endif

    ! Check if the asked for cranking options are allowed by the CONFIG file.
    do i=1,3
      if(Omega(i) .ne. 0.0d0 .or. CrankValues(i) .ne.  0.0d0) then
        NotFound = .true.
$NTR        do j=1,cranklen
$NTR          c = crankdirections(j)
$NTR          if( c == i ) NotFound = .false.
$NTR        enddo
        if(NotFound) then
          print *, 'Disallowed cranking frequency.'
          stop
        endif
      endif
    enddo
  end subroutine readcranking

  subroutine updateAM(R)
    !---------------------------------------------------------------------------
    ! Calculate the total angular momentum and cranking energies.
    !---------------------------------------------------------------------------
$NTR    use Moments, only : cutoff
    ! We only import this if time-reversal is not conserved, otherwise
    ! the compiler complains

    type(DensityVector), intent(in) :: R
$NTR    integer   :: B, N, wave, si, i, c, it
$TR real(KIND=dp) :: trash
    ! Saving all of the history for convergence analysis ...
    angmomold       = totalangmom
    angmomold_dens  = totalangmom_dens
    angmomold_cut   = totalangmom_cut

    ! ... and resetting the current values
    totalangmom = 0.0 ; totalangmom_dens = 0.0d0 ; totalangmom_cut = 0.0d0
    J2_sp       = 0.0

$TR trash = R%D_I_I(1,1) ! To stop compiler complaints when time-reversal is conserved

! $NTR    si = 0
! $NTR    do B=1,8
! $NTR      N = HFBlocks(B) ; if(N .eq. 0) cycle
! $NTR      do wave = 1, N
! $NTR        do i = 1, cranklen
! $NTR          c  = crankdirections(i)
! $NTR          if(pairingtype.ne.2) then
! $NTR            TotalAngMom(c) = TotalAngMom(c) + rho_can(si+wave) * HF_J (c,si+wave)
! $NTR            J2_sp      (c) = J2_sp      (c) + rho_can(si+wave) * HF_J2(c,si+wave)
! $NTR          else
! $NTR            TotalAngMom(c) = TotalAngMom(c) + rho_can(si+wave) * CAN_J (c,si+wave)
! $NTR            J2_sp      (c) = J2_sp(c)       + rho_can(si+wave) * CAN_J2(c,si+wave)
! $NTR          endif
! $NTR        enddo
! $NTR      enddo
! $NTR      si = si + N
! $NTR    enddo
! $NTR
! $NTR    !-------------------------------------------------------------------------
! $NTR    ! And now we integrate the current density and spin density.
! $NTR    do it=1,2
! $NTR      ! Spin part
! $NTR      TotalAngMom_dens(3) = TotalAngMom_dens(3) + &
! $NTR      &                     0.5 *                sum(R%D_I_S(:,3,it))
! $NTR      TotalAngMom_cut(3) = TotalAngMom_cut(3) + &
! $NTR      &                     0.5 * sum( Cutoff(:,it)* R%D_I_S(:,3,it))
! $NTR
! $NTR      do i=1, nx*ny*nz
! $NTR        TotalAngMom_dens(3) = TotalAngMom_dens(3) &
! $NTR        & - meshgrid(i,2) * R%C_I_N(i,1,it) + meshgrid(i,1) * R%C_I_N(i,2,it)
! $NTR
! $NTR        TotalAngMom_cut(3)  = TotalAngMom_cut(3) + cutoff(i,it) * &
! $NTR        & (- meshgrid(i,2) * R%C_I_N(i,1,it) + meshgrid(i,1) * R%C_I_N(i,2,it))
! $NTR      enddo
! $NTR    enddo
! $NTR    TotalAngMom_dens = TotalAngMom_dens * dv
! $NTR    TotalAngMom_cut  = TotalAngMom_cut  * dv

$NTR    !-----------------------------------------------------------------------
$NTR    ! The contribution of the cranking constraint to the total Routhian
$NTR    crankenergy     = - omega * TotalAngMom
$NTR    crankenergy_cut = - omega * TotalAngMom_cut

  end subroutine updateAM

  subroutine printcranking_init
    !---------------------------------------------------------------------------
    ! Print information on the cranking constraints imposed on the calculation
    ! at the start of the program.
    !---------------------------------------------------------------------------
    1 format ( 29('-'), ' Cranking information ', 29('-'))
    2 format ( ' Cranking frequencies:   ', /, &
    &          '    Omega_X = ', f15.3, /,  &
    &          '    Omega_Y = ', f15.3, /,  &
    &          '    Omega_Z = ', f15.3)
    21 format ('    => read from STDIN ')
    22 format ('    => read from FILE  ')
    3 format ( ' Cranking target values: ', /, &
    &          '    J_X     = ', f15.3, /,  &
    &          '    J_Y     = ', f15.3, /,  &
    &          '    J_Z     = ', f15.3)
   31 format ( ' Cranking scale factors: ', /, &
    &          '    J_X     = ', f15.3, /,  &
    &          '    J_Y     = ', f15.3, /,  &
    &          '    J_Z     = ', f15.3)
    4 format ( ' Cranking types: ', /, &
    &          '    J_X     = ', i15, /,  &
    &          '    J_Y     = ', i15, /,  &
    &          '    J_Z     = ', i15 )
    5 format ( ' Cranking on the basis of INTEGRATION OF DENSITIES')
   51 format ( ' Cranking on the basis of SUMMED SPWF ANGULAR MOMENTUM')

    print 1
    print 2, Omega
    if(continuecrank) then
      print 22
    else
      print 21
    endif
    print 3, CrankValues
    print 31, CrankScaleFactor
    print 4, Cranktype
    if(crank_smooth) then
      print 5
    else
      print 51
    endif

  end subroutine printcranking_init

  subroutine PrintCranking(R)
    !---------------------------------------------------------------------------
    ! Prints all kinds of information about the expectation value of the
    ! angular momentum operator and all kinds of angles.
    !---------------------------------------------------------------------------
    character(len=1), parameter     :: dir(3) = (/'x', 'y', 'z'/)
    integer                         :: i
$NTR integer                        :: j
$TR  real(KIND=dp)                  :: trash
    type(DensityVector), intent(in) :: R
$NTR    logical           :: found

    1 format (2x,99('_') )
   10 format (2x,99('-'))
    2 format (30('-'), ' Angular Momentum (hbar) ',46('-') )
    3 format (15x, 'Spwfs(*)  ',7x, 'Desired', 10x, 'Omega', 12x, 'E (MeV)' 12x,'Densit. ')
   31 format (15x, 'Densit.(*)',7x, 'Desired', 10x, 'Omega', 12x, 'E (MeV)' 12x,'Spwfs   ')
    4 format (3x,'J_',a1,'   ','|', 5f17.10 )
   41 format (3x,'Size  |', 3f17.10,17x,1f17.10)
$NTR    5 format (2x,' _______________________________________________________' )
$NTR    6 format (3x,'Open spin')
$NTR    7 format (15x, 'Neutrons', 3x, 'Protons')
$NTR    8 format (3x,a1,1x,'|',3x,'|',4f17.10)

    print 2
    print *

    if(crank_smooth) then
     print 31
    else
     print 3
    endif

    !---------------------------------------------------------------------------
    ! Information on the total angular momentum
    print 1
    do i=1,3
      if(crank_smooth) then
        print 4, dir(i), TotalAngMom_dens(i), CrankValues(i), Omega(i),        &
        &                CrankEnergy(i), TotalAngMom     (i)
      else
        print 4, dir(i), TotalAngMom(i), CrankValues(i), Omega(i),   &
        &                CrankEnergy(i), TotalAngMom_dens(i)
      endif
    enddo
    print 1
    if(crank_smooth) then
      print 41, sqrt(sum(totalangmom_dens**2)), 0.0, &
      &         sqrt(sum(omega**2))      , sqrt(sum(totalangmom**2))
    else
      print 41, sqrt(sum(totalangmom**2)), 0.0, &
      &         sqrt(sum(omega**2))      , sqrt(sum(totalangmom_dens**2))
    endif
    print 10

    !---------------------------------------------------------------------------
    ! Information on the spin density
! $TR trash = R%D_I_I(1,1)
!
! $NTR    print *
! $NTR    print 6
! $NTR    print 7
! $NTR    print 5
! $NTR    do i=1,3
! $NTR      found = .false.
! $NTR      do j=1,cranklen
! $NTR        if(i .eq. crankdirections(j)) found = .true.
! $NTR      enddo
! $NTR      if(found) then
! $NTR        ! There is a possibility for total spin in this Cartesian direction.
! $NTR        print 8     , dir(i), 0.5*sum(R%D_I_S(:,i,1))*dv, &
! $NTR        &                     0.5*sum(R%D_I_S(:,i,2))*dv
! $NTR      else
! $NTR        ! Spin is restricted in this particular direction
! $NTR        print 8     , dir(i), 0.0d0,0.0d0
! $NTR      endif
! $NTR    enddo
! $NTR    print 5
  end subroutine PrintCranking

  function crank_spin_potential() result(spot)
    !---------------------------------------------------------------------------
    ! The cranking constraint contributes to the potential F_I_S, associated 
    ! with the spin density D_I_S:
    !
    !     F_I_S(r) => F_I_S(r) - \frac{1}{2} \omega f_cut(r)
    !
    ! where f_cut(r) is the cutoff factor also employed for the multipole
    ! constraints. The factor 1/2 is present because the J_spin = 1/2 Pauli
    ! sigma matrix.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Output:
    !     spot : the contribution to the spin potential due to cranking
    !---------------------------------------------------------------------------
$NTR    use Moments, only : cutoff

    real(KIND=dp), allocatable :: spot(:,:,:)
$NTR    integer :: i, it, c

    allocate(spot(nx*ny*nz,3,4)) ; spot = 0.0d0
!
! $NTR    do i=1, cranklen
! $NTR      c           = crankdirections(i)
! $NTR      do it=1,2
! $NTR        spot(:,c,it) = - 0.5 * omega(c) * cutoff(:,it)
! $NTR      enddo
! $NTR    enddo
!
! $NTR    spot(:,:,3) = spot(:,:,1) + spot(:,:,2)
! $NTR    spot(:,:,4) = spot(:,:,1) - spot(:,:,2)

    return
  end function crank_spin_potential

  function crank_current_potential() result(jpot)
    !---------------------------------------------------------------------------
    ! The cranking constraint contributes to the potential G_I_N, associated
    ! with the current density D_I_N:
    !
    !       G_I_N => G_I_N - f_cut(r) \vec{\omega} x \vec{r}
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Output:
    !     jpot : the contribution to the current potential due to cranking
    !---------------------------------------------------------------------------
    use Moments, only : cutoff

    real(KIND=dp), allocatable :: jpot(:,:,:)
    integer                    :: it, mu, indices(2) , nu, ka

    allocate(jpot(nx*ny*nz,3,4)) ; jpot = 0.0d0
!     do mu=1, 3
!       indices = vector_product(mu)
!       nu = indices(1)
!       ka = indices(2)
!       do it=1,2
!           jpot(:,mu, it) = -    cutoff(:,it) * &
!           &            (omega(nu) * meshgrid(:,ka) - omega(ka) * meshgrid(:,nu))
!       enddo
!     enddo
!
!     jpot(:,:,3) = jpot(:,:,1) + jpot(:,:,2)
!     jpot(:,:,4) = jpot(:,:,1) - jpot(:,:,2)

    return
  end function crank_current_potential

  subroutine ReadjustCranking()
    !---------------------------------------------------------------------------
    ! Readjust the cranking constraint(s)
    !---------------------------------------------------------------------------
    integer                 :: i
    real(KIND=dp)           :: value
    character(len=1), parameter  :: dir(3) = (/'x', 'y', 'z'/)
   11 format ( '------------------------------------------------------')
   12 format ( ' Constraint on J_', a1,  ' has no intensity parameter.')
   13 format ( ' Heuristic taken from J2_sp:', 1es12.5)

    do i=1,3
        select case(CrankType(i))
        case(0)
          ! Nothing to be done; constant omega-cranking
          Omega_prev(i) = Omega(i)
        case(1)
          if(CrankIntensity(i) .eq. 0.0_dp) then
            !-------------------------------------------------------------------
            ! If left to zero by the user, the code uses this heuristic to guess
            ! an intensity for the cranking constraint. This works very well for
            ! collective cranking. Although still decent, you might want to
            ! manually set intensity_x when cranking around a symmetry axis.
            print 11
            print 12, dir(i)
            CrankIntensity(i) =  1.0d0 / J2_sp(i)
            print 13, CrankIntensity(i)
            print 11
          endif
          ! Actual readjustment of the constraint
          if(crank_smooth) then
            value = TotalAngMom_dens(i)
          else
            value = TotalAngMom(i)
          endif
          ! Save the previous value
          Omega_prev(i) = Omega(i)
          Omega(i) = Omega(i)- CrankIntensity(i)*(value-CrankValues(i))

          ! Debugging printout
!          print '(" ReadjustCranking ",i2,7f12.5,l3)',            &
!           & i, Omega(i), Jtotal, J2_sp(i), CrankIntensity(i),  &
!           & TotalAngMom(i),CrankValues(i),                       &
!           & CrankIntensity(i)*(TotalAngMom(i) - CrankValues(i))
        end select
    enddo
  end subroutine ReadjustCranking

end module
