!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
!    Copyright (C) 2026 W. Ryssens and M. Bender
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Affero General Public License as published
!    by the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!===============================================================================
module BCS
 !==============================================================================
 !  
 ! Module implementing the routines for the solution of the BCS equations.
 !
 !==============================================================================

 use vectors, only: PotentialVector
 use wavefunctions
 use pairingcutoffs

 implicit none
 
 !------------------------------------------------------------------------------
 real(KIND=dp) :: FermiPrec = 1d-9
 integer       :: maxBCSiter = 500
 logical       :: ConstantGap = .false.
 
 real(KIND=dp) :: BCSdispersion(2)
 !------------------------------------------------------------------------------
 ! Values of the BCS gaps Delta.
 real(KIND=dp), allocatable :: BCSGaps(:)
 ! BCS quasiparticle energies
 real(KIND=dp), allocatable :: BCSqps(:)
 !------------------------------------------------------------------------------
 ! BCS occupations, i.e. 2 * v_i^2 (factor 2 due to time-reversal)
 ! - - - - - - - - - - - - - - - --  - - - - - - - - - - - - - - - - - - - - - -
 ! IMPORTANT NOTE:
 !  These are NOT ALWAYS equal to the diagonal matrix elements of the 
 !  density matrix. This is ONLY true at zero temperature. They are ALWAYS
 !  equal to 2 * v_i**2, but NOT always equal to rho_ii !
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable :: BCSoccupations(:)
 !------------------------------------------------------------------------------
 ! BCS occupation factors f, typically defined in the context of 
 ! finite-temperature calculations, but also useful in the context of 
 ! Equal Filling Approximation blocking calculations.
 real(KIND=dp), allocatable :: BCSf(:)
 !------------------------------------------------------------------------------
 procedure(delta_action_dummy), pointer :: delta_action_BCS

  interface
   function delta_action_dummy(psi,&
$N1DELTA                    &      dpsi, &
$N2DELTA                    &            ddpsi, &
$N3DELTA                    &                   dddpsi, &
$SYMDELTA                   &                          sx,sy,sz, &
&                                                            iso, onthefly, F) &
                                                                result(deltapsi)
      !-------------------------------------------------------------------------
      ! Dummy function to allow this module to access the functional.f90 module 
      ! routine to calculate the "action of" Delta.
      !-------------------------------------------------------------------------
      import PotentialVector ! explicit import statement, otherwise the 
                             ! interface would be invalid
      real*8, intent(in)                :: psi(:,:)
      type(PotentialVector), intent(in) :: F
$N1DELTA      real*8, intent(inout) ::   dpsi(:,:,:)
$N2DELTA      real*8, intent(inout) ::  ddpsi(:,:,:)
$N3DELTA      real*8, intent(inout) :: dddpsi(:,:,:)
$SYMDELTA     integer, intent(in)   :: sx(:),sy(:),sz(:)
      integer, intent(in)   :: iso
      real*8, allocatable   :: deltapsi(:,:)
      logical, intent(in)   :: onthefly
   end function
  end interface

  interface
    !-------------------------------------------------------------------------
    ! The following procedures are implemented in submodule BCS_quantities
    ! (src_orig/BCS_quantities.f90). Only their interfaces live here.
    !-------------------------------------------------------------------------
    module subroutine BCSQPEnergies(Fermi)
      real(KIND=dp), intent(in) :: Fermi(2)
    end subroutine BCSQPEnergies

    module subroutine calcBCSOccupations(Fermi)
      real(KIND=dp), intent(in) :: Fermi(2)
    end subroutine calcBCSOccupations

    module subroutine calcBCSdispersion(rho_can, kappa_can)
      real(KIND=dp), intent(in) :: rho_can(:), kappa_can(:)
    end subroutine calcBCSdispersion

    module subroutine clean_BCS()
    end subroutine clean_BCS

    module function average_gap_BCS() result(gap)
      real(KIND=dp) :: gap(2,2)
    end function average_gap_BCS

    module function Order(energies, nwf) result(Indices)
      integer, intent(in)       :: nwf
      real(KIND=dp), intent(in) :: energies(nwf)
      integer                   :: Indices(nwf)
    end function Order
  end interface

contains
 
 subroutine uv_from_occupation(occ, u, v)
  !-----------------------------------------------------------------------------
  ! Simple function that calculates BSC u and v factors in a numerically
  ! safe way from the BCS occupation v^2.
  !
  ! Attention: only VALID if u and v can be assumed to be real.
  !
  ! Input:
  !      occ  : 2 v^2, i.e. what is stored in BCSoccupations
  ! Output:
  !      u, v : bcs factors
  !-----------------------------------------------------------------------------

  real(KIND=dp), intent(in) :: occ
  real(KIND=dp), intent(out) :: u,v

  if(occ/2 .gt. 0.0d0) then
    v = sqrt(occ/2)
  else
    v = 0.0d0
  endif
  
  u = 1 - occ/2
  if(u .gt. 0.0d0) then
    u = sqrt(u)
  else
    u = 0.0d0
  endif
 end subroutine uv_from_occupation
 
 subroutine solvepairing_BCS(fermi, rho_can, kappa_can, qpenergies, gas,       &
 &                             BlockType,Blockindices, blocklowest, blocked_qps)
  !-----------------------------------------------------------------------------
  ! Driver routine for the solving of the BCS equations.
  !
  ! |----
  ! | Do until the Fermi energy is stationary
  ! |   1) Calculate the BCS quasiparticle energies
  ! |        using the matrix elements of h and delta already calculated!
  ! |   2) Construct the correct configuration with either
  ! |       *) ordinary, ground-state BCS
  ! |       *) finite-temperature BCS
  ! |       *) equal-filling blocking BCS
  ! |   3) Calculate the Fermi energy with analytic formula
  ! |----
  !     3) Calculate the entries in rho_pairing and kappa_pairing
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(inout) :: fermi(2)
  real(KIND=dp), intent(inout) :: rho_can(:)
  real(KIND=dp), intent(inout) :: kappa_can(:), qpenergies(:)
  real(KIND=dp)                :: oldfermi(2), fac
  integer                      :: iter, wave
  integer, intent(in)          :: gas

  ! Options for the determination of a blocking configuration
  integer, intent(in)          :: Blockindices(:)
  integer, intent(in)          :: BlockType
  character(len=2), intent(in) :: BlockLowest(:)
  integer, allocatable         :: blocked_qps(:)
  
  1 format('--------------------------------------------------------',/, &
    &        ' BCS calculations did not converge.',/,  &
    &        ' Pairing iterations:  ', i4,/,           &
    &        ' Old Fermi         :  ', 2f15.10,/,      &
    &        ' New Fermi         :  ', 2f15.10,/,      &
    &        ' Difference        :  ', 2f15.10,/,      &
    &        ' Particle deviation:  ', 2f15.10,/,      &
    &        '-----------------------------------------------------------')

  if(.not.allocated(Bcsf)) then
    allocate(bcsf(nwt)) ; bcsf = 0
  endif

  !-----------------------------------------------------------------------------
  ! Iteration start
  do iter =1, maxBCSiter
    oldfermi = fermi
    ! Calculate the quasiparticle energies
    call BCSQPEnergies(Fermi) 
    BCSf = BCSconstructconfiguration(gas, BlockType,Blockindices, blocklowest, & 
                                                                    blocked_qps)
    call BCSFindFermiEnergy(Fermi)
    ! Check for convergence
    if( all(abs(fermi - oldfermi).lt.FermiPrec)) then
      exit
    elseif(iter.eq.maxBCSiter) then
      call calcBCSoccupations(Fermi)
      if(MPI_RANK.eq.0) print 1, iter, oldfermi, fermi, fermi-oldfermi, &
      &                          sum(BCSoccupations(1:nwn)) -neutrons,&
      &                          sum(BCSoccupations(nwn+1:))-protons
    endif          
  enddo
  ! Iteration end  
  call BCSQPEnergies(Fermi)
  call calcBCSoccupations(Fermi)
  !-----------------------------------------------------------------------------
  ! Rho_pairing is diagonal for a BCS calculation
  rho_can = 0.0d0
  do wave=1,nwt
    ! Occupations are 
    !   n_a = f_i + v_i^2 (1 - 2 * f_i)
    fac = BCSf(wave) 
    ! Note that there is already a factor two due to timereversal in    
    ! BCSoccupations, but not in the first term in the formula above.
    rho_can(wave) = 2.0d0*fac + BCSoccupations(wave) * (1 - 2.0d0*fac)

    ! Failsafe
    ! 06/05/'25: disabled because it lead to issues in pasta calculations
    ! if(rho_can(wave).lt.0.0) then
    !   rho_can(wave) = 0.0
    ! endif
  enddo

  !-----------------------------------------------------------------------------
  ! Kappa is in its canonical form for a BCS calculation.
  ! However, we are not storing the time-reversed partners in this case, so
  ! put the matrix elements (i,ibar) on the diagonal anyway.
  !
  ! u * v  = 0.5 * Delta/(sqrt(epsilon**2 + Delta**2))
  !-----------------------------------------------------------------------------
  kappa_can = 0.0d0
  do wave=1,nwt
     ! At finite temperature, the elements of kappa are
     ! kappa_i\bar{i} = u_i v_i ( 1 - 2 * f_i )
     fac             = BCSf(wave)
     kappa_can(wave) = 0.5d0 * BCSgaps(wave)/(BCSqps(wave)) * (1 - 2.0d0*fac)
  enddo

  ! Qpenergies in this case are the BCSqpenergies
  Qpenergies = BCSqps

  ! Side effect: calculate the dispersion 
  call calcBCSdispersion(rho_can, kappa_can)

 end subroutine solvepairing_BCS
 
 subroutine CalcBCSGaps(fermi, stabfactor, F)
    !---------------------------------------------------------------------------
    ! Calculate the BCS pairing gaps.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !  fermi     : chemical potential for both nucleon species
    !  stabfactor: factor to use in the stabilised pairing
    !  F         : set of mean-field potentials
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)         :: fermi(2), stabfactor(2)
    type(PotentialVector), intent(in) :: F 
    integer                      :: wave, iso, wave_global
    real(KIND=dp)                :: trash(2)
    real(KIND=dp), allocatable   :: deltapsi(:,:)
    ! Technical note: deltapsi HAS to be allocatable as opposed to an automatic
    !                 array, because the result of the function delta_action
    !                 is allocatable. This cannot be changed, because delta_action
    !                 is assigned through pointer remapping which requires an
    !                 interface (see delta_action_dummy at the top of this file).
    !                 Said interface cannot be defined in terms of runtime
    !                 variables such as the number of mesh points. Switching to an
    !                 automatic array here will lead to memory leaks with IFORT
    !                 compilers.
#if(USE_MPI>0)
    integer                      :: mpi_err
#endif
    ! trash statement to stop the compiler complaining about unused dummy 
    ! variables
    trash = fermi

    if(ConstantGap) then  
      ! Constantgap pairing
      do wave=1,nwt
          BCSGaps(wave) = 2.0d0 * PCutoffs(wave)**2
      enddo
    else
       ! Use the delta_action function to calculate the matrix elements
       BCSgaps = 0.0d0 ! zeroing to be able to call MPI_ALLREDUCE later
       do wave=1,nwt_local
            wave_global = spwf_map(wave)
            if(wave_global .le. nwn) then
                iso = -1
            else
                iso = +1
            endif

            deltapsi = delta_action_BCS(  hfpsi(:,:,wave) ,             &
$N1DELTA    &                            hfdpsi(:,:,:,wave),            &
$N2DELTA    &                           hfddpsi(:,:,:,wave),            &
$N3DELTA    &                          hfdddpsi(:,:,:,wave),            &
$SYMDELTA   &              sx(:,wave), sy(:,wave), sz(:,wave),          &
            &                                             iso,.false.,F)
 
 
            BCSgaps(wave_global) =  sum(hfpsi(:,:,wave)*deltapsi)*dv*          &
            &             Pcutoffs(wave_global)**2 * (1 + stabfactor((iso+3)/2))
       enddo
#if(USE_MPI>0)
       call MPI_ALLREDUCE(MPI_IN_PLACE, BCSgaps,nwt,MPI_REAL8,MPI_SUM,         &
       &                                                 MPI_COMM_WORLD,mpi_err)
#endif
    endif
  end subroutine CalcBCSGaps

  subroutine BCSFindFermiEnergy (Fermi)
   !----------------------------------------------------------------------------
   ! Implements the analytical formula for the Fermi energy, for fixed 
   ! quasiparticle energies and occupation factors f.
   !  
   !
   ! Starting by demanding that
   !     N     = 2 sum_(i>0) [ f_i  + v_i^2 ( 1 - 2 f_i) ]
   ! and  
   !     v_i^2 = 1/2 ( 1 - (eps_i - mu)/E^qp_i ) 
   !
   ! we obtain 
   ! 
   !  mu = (N - a)/b
   !  
   ! with
   !   a = sum_(i>0) (1  - eps_i/E^qp_i * (1-2f_i))  
   !   b = sum_(i>0) ((1 - 2 f_i)/E^qp_i)  
   !----------------------------------------------------------------------------
    real(KIND=dp),intent(inout) :: Fermi(2)
    real(KIND=dp)               :: LambdaSums(2,2), eqp, nom, Particles(2), fac
    integer                     :: wave, it

    Particles(1) = Neutrons; Particles(2) = Protons

    LambdaSums=0.0_dp 
    !Sum the quasiparticle energies to get the correct lambda.
    do wave=1,nwt
      it = 1
      if(wave .gt. nwn) it = 2 
     
      eqp = BCSqps(wave)
      nom = spenergies(wave)

      fac =  BCSf(wave) 
      lambdasums(1,it)= lambdasums(1,it) +(1.0_dp - nom/eqp*(1-2*fac))
      lambdasums(2,it)= lambdasums(2,it) +       1.0_dp/eqp*(1-2*fac)
    enddo
    do it=1,2
        Fermi(it)=(particles(it)  - lambdasums(1,it))/lambdasums(2,it)
    enddo
  end subroutine BCSFindFermiEnergy

  function BCSConstructConfiguration(gas, blocktype,Blockindices, blocklowest, & 
  &                                  blocked_qps) result(f)

    !---------------------------------------------------------------------------
    ! Construct the correct pairing configuration as a function of various 
    ! user options. 
    !  (a) Zero-temperature, no blocking  : all the BCS-f factors are zero.
    !  (b) Zero-temperature, EFA blocking : f_k = 1/2 for the blocked state
    !  (c) Finite-temperature, no blocking: the f_k are given by a fermi
    !                                       function (modulo different options 
    !                                       for the BCS gas)
    !  (d) Finite-temperature, blocking   : nothing implemented
    !---------------------------------------------------------------------------
    real*8              :: f(nwt), occ, qpmin
    integer             :: wave, NB, i, ind, si, N, B,  qpb, c

    integer, intent(in)          :: Blockindices(:)
    integer, intent(in)          :: BlockType, gas
    character(len=2), intent(in) :: BlockLowest(:)
    integer, allocatable, intent(inout) :: blocked_qps(:)
    integer, allocatable         :: proton_block(:), neutron_block(:)
    integer, allocatable         :: indices(:), toblock(:)
  
    f = 0 ; qpb = 0
    !---------------------------------------------------------------------------
    ! Trash statements to fool CRAY compilers. If this is not here, the compiler
    ! complains about the array being used before being allocated as soon as 
    ! the optimisation level is equal to -O1 or above. I guess this is related 
    ! to the highly-nested nature of this routine ...
    allocate(indices(1))       ; deallocate(indices)
    allocate(proton_block(1))  ; deallocate(proton_block)
    allocate(neutron_block(1)) ; deallocate(neutron_block)
    allocate(toblock(1))       ; deallocate(toblock)
    !--------------------------------------------------------------------------

    if(allocated(blocked_qps)) deallocate(blocked_qps)
    !---------------------------------------------------------------------------
    ! Zero-temperature
    if(inversetemp .lt. 0) then
      occ = 0.0d0
      select case(Blocktype)
      case(0)
        ! No blocking, all the f are zero
        return
      case(1,2,5)
        ! Time-reversal breaking blocking asked for, impossible to do in BCS
        call stp('The code cannot perform true blocking in BCS.')
      case(3,4,6)
        ! Equal filling blocking
        occ = 0.5d0
      end select

      select case(Blocktype)
      case(3)
        !-----------------------------------------------------------------------
        ! We search for a specific configuration, i.e. a quasiparticle with
        ! a specific sp index. In a BCS calculation, this is trivial. 
        NB = size(blockindices)      
        allocate(blocked_qps(NB))
        do i=1, NB
          ind    = blockindices(i)
          f(ind) = occ
          blocked_qps(i) = ind
        enddo
      case(4)
        !-----------------------------------------------------------------------
        ! We search for the lowest qp with specific quantum numbers.
        
        allocate(proton_block(5))  ; proton_block  = 0 
        allocate(neutron_block(5)) ; neutron_block = 0
        allocate(toblock(8))       ; toblock       = 0

        ! Count all qps that were asked for
        do i=1, size(Blocklowest)
          select case (Blocklowest(i))
          case('n+')
              neutron_block(1) = neutron_block(1) + 1 
          case('n-')
              neutron_block(3) = neutron_block(3) + 1 
          case('p+')
              proton_block(1)  = proton_block(1)  + 1
          case('p-')
              proton_block(3)  = proton_block(3)  + 1
          case('n0')
              neutron_block(5) = neutron_block(5) + 1
          case('p0')
              proton_block(5)  = proton_block(5)  + 1
          end select
        enddo
        toblock(1:4) = neutron_block(1:4)
        toblock(5:8) = proton_block(1:4)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! First see where the lowest qp energies are for those qps in block(5)
        if(neutron_block(5).ne.0) then
          do i = 1, neutron_block(5)
            qpmin = 10000000
            si    = 0
            qpb   = 0
            do B=1,4
              N = HFBlocks_global(B) ; if (N.eq.0) cycle
              if(bcsqps(si+toblock(B)+1) .lt. qpmin) then
                qpmin = bcsqps(si+toblock(B)+1)
                qpb   = B
              endif
              si = si +   N
            enddo
            toblock(qpb) = toblock(qpb) + 1
          enddo
        endif

        if(proton_block(5).ne.0) then
          do i = 1, proton_block(5)
            qpmin = 10000000
            si    = sum(HFBlocks_global(1:4)) 
            qpb   = 0
            do B=5,8
              N = HFBlocks_global(B) ; if (N.eq.0) cycle

              if(bcsqps(si+toblock(B)+1) .lt. qpmin) then
                qpmin = bcsqps(si+toblock(B)+1)
                qpb   = B
              endif
              si = si +   N
            enddo
            toblock(qpb) = toblock(qpb) + 1
          enddo
        endif

        allocate(blocked_qps(sum(toblock))) ; blocked_qps = 0
        si = 0
        c  = 0
        do B=1,8
            N = HFBlocks_global(B) ; if(N.eq.0) cycle
            allocate(indices(N))
            indices = Order(BCSqps(si+1:si+N),N)
            do i=1, toblock(B)
              f(si+indices(i)) = occ
              c = c+1
              blocked_qps(c) = si+indices(i)
            enddo            
            si = si + N
            deallocate(indices)
        enddo
        deallocate(proton_block, neutron_block) 
        !-----------------------------------------------------------------------
      end select
    !---------------------------------------------------------------------------
    ! Finite-temperature
    else
      if(blocktype.ne. 0) then
        call stp('Cannot do finite-temperature BCS with blocking.')
      endif
      
      ! Select occupations based on the type of treatment of the gas
      select case(gas)
      case(0)
        ! No special treatment of the gas
        f = 1.d0/(1.d0 + exp(inversetemp * BCSqps))
      case(1)
        ! Not implemented!
        call stp('gastype = 1 is not implemented in the BCS module.')
      case(2)
        ! Only take into account the bound states
        do wave=1,nwt
          if(spenergies(wave) .lt. 0) then
            f(wave) = 1.d0/(1.d0 + exp(inversetemp * BCSqps(wave)))
          else 
            f(wave) = 0d0
          endif
        enddo
      case DEFAULT
        call stp('Unknown value for "gas" in the BCS module.')
      end select
    endif
    return
  end function BCSConstructConfiguration
  
!===============================================================================
!  Never to be used function to define an interface for delta_action
!===============================================================================   
! function delta_action_dummy(psi, dpsi, ddpsi, dddpsi, sx,sy,sz,iso, onthefly) &
!                                                                result(deltapsi)
!    !---------------------------------------------------------------------------
!    ! Dummy function to allow this module to acces the functional.f90 module 
!    ! to acces the information on the acces of deltas.
!    !---------------------------------------------------------------------------
!    
!    real(KIND=dp), intent(in)    :: psi(mv,4)  
!    real(KIND=dp), intent(inout) :: dpsi(mv,3,4),ddpsi(mv,6,4), dddpsi(mv,10,4)
!    integer, intent(in)       :: sx(4),sy(4),sz(4),   iso
!    logical, intent(in)       :: onthefly
! end function
  
end module
