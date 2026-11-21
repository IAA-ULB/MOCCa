submodule (functional) functional_potentials
!==============================================================================
!
! Implementation submodule of module functional: functional_potentials
!  Interface block lives in src_orig/functional.f90.
!
!==============================================================================
implicit none
contains

  module function calcPotentials(R, Fread) result(F)
    !---------------------------------------------------------------------------
    ! Calculate all of the Skyrme potentials and the Coulomb potential
    ! corresponding to the input densityvector R.
    !
    ! Input: 
    !    R        : density-vector for the calculation of the fields
    !    Fread    : optional Potentialvector.
    !               If present, this potentialvector is used to initialize 
    !               the calculated potentialvector; only the fields that 
    !               are detected to be completely zero are then recalculated
    ! Output: 
    !    F        : a fresh potentialvector 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 
    ! Notes
    ! (1)  this function does not do any preconditioning, mixing or anything 
    !      else to the calculated fields. They are calculated "as is" from
    !      the density vector passed in.
    !
    ! (2)  The Coulomb potential is part of the field-vector, and hence 
    !      gets calculated here as well. 
    !---------------------------------------------------------------------------
    use Coulombmod, only : solve_coulomb
    use Coulombmod, only : Coulomb_read_from_file
    use pairing_strengths, only : vmicro, vmicro_stored
    use moments

    type(DensityVector), intent(in)             :: R
    type(PotentialVector), intent(in), optional :: Fread
    type(PotentialVector)                       :: F
    real(KIND=dp)                               :: radius
    integer                                     :: i

    call start_timer(T_potentials)

    ! Copy all the relevant information if this information was presented
    if(present(Fread)) F = Fread

    ! Signal that microscopic pairing strengths have to be recalculated
    vmicro_stored = .false.

    !---------------------------------------------------------------------------
    ! Calculate all contributions to the potentials due to the Skyrme EDF
$CALCPOTENTIALS

    !---------------------------------------------------------------------------
    ! Calculate the Coulomb potentials
    !---------------------------------------------------------------------------
    if((.not. present(Fread)) .or. (.not. Coulomb_read_from_file)) then
        call solve_coulomb(R,F,sx_rho, sy_rho, sz_rho)
    endif
    !---------------------------------------------------------------------------
    ! Simulation of spherical boundary conditions
    if(simulate_spherical_bc) then 
      if(reduZ .eq. 1) then
        if(nx.ne.nz) then
          call stp('Spherical boundary conditions are only meaningful for cubic meshes.')
        endif
      else
        if(2*nx.ne.nz) then
          call stp('Spherical boundary conditions are only meaningful for cubic meshes.')
        endif
      endif
      if(nx.ne.ny) then
          call stp('Spherical boundary conditions are only meaningful for cubic meshes.')
      endif

      do i=1,mv 
        radius = sqrt(sum(meshgrid(i,:)**2))
        if(radius .gt. (nx + 0.5d0) * dx) then 
          F%F_I_I(i,1) = F%F_I_I(i,1) + 1000.0d0
          F%F_I_I(i,2) = F%F_I_I(i,2) + 1000.0d0
          F%F_I_I(i,3) = F%F_I_I(i,1) + F%F_I_I(i,2)
          F%F_I_I(i,4) = F%F_I_I(i,1) - F%F_I_I(i,2)
        endif
      enddo
    endif

    call stop_timer(T_potentials)

  end function calcPotentials

#if($FAM == 1)

  module subroutine calc_perturbed_potentials(  R, dRs, dRa, dR_pp_plus, dR_pp_minus, &
  &                                         dFs, dFa, dF_pp_plus, dF_pp_minus)
    !---------------------------------------------------------------------------
    ! Calculate the linearised response of the potentials (dFs, dFa) around
    ! a set of mean-field densities (R) that are affected by perturbations
    ! (dRs, dRa).
    !
    ! There is a clear separation between so-called symmetric and antisymmetric
    ! quantities, both for the densities and potentials. The TOTAL perturbation
    ! is in both cases the sum of the symmetric and antisymmetric parts, but
    ! these parts have different symmetry properties such that we cannot simply
    ! sum them everywhere.
    !
    ! For the potentials, the expressions for the symmetric and antisymmetric
    ! part are identical; this routine is thus not much more than a wrapper
    ! around two function calls to calc_perturbed_potentials.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   R   : density-vector containing the mean-field densities
    !   dRs : density-vector containing the symmetric perturbing particle-hole densities
    !   dRa : density-vector containing the antisymmetric perturbing particle-hole densities
    !   dR_pp_plus : density-vector containing the perturbing particle-particle densities deduced from kappa
    !   dR_pp_minus: density-vector containing the perturbing particle-particle densities deduced from kappa^*
    !
    ! Output:
    !   dFs : potential-vector containing the symmetric part of the linearised
    !         response of the mean-field potentials
    !   dFa : potential-vector containing the antisymmetric part of the linearised
    !         response of the mean-field potentials
    !   dF_pp_plus : potential-vector containing the linearised response of the pp potentials with kappa
    !   dF_pp_minus: potential-vector containing the linearised response of the pp potentials with kappa^*
    ! 
    !---------------------------------------------------------------------------
    use CoulombMod, only : solve_coulomb_linear_response

    type (DensityVector), intent(in) :: R, dRs, dRa, dR_pp_plus, dR_pp_minus
    type (PotentialVector)           :: dFs, dFa, dF_pp_plus, dF_pp_minus
    
    dFs  = calc_perturbed_potentials_oneoff(R,dRs)  ! perturbation of symmetric ph potentials
    dFa  = calc_perturbed_potentials_oneoff(R,dRa)  ! perturbation of antisymmetric ph potentials

    call solve_coulomb_linear_response(R, dRs,dFs,sx_rho        ,sy_rho        ,sz_rho)          ! Coulomb response to symmetric ph densities
    call solve_coulomb_linear_response(R, dRa,dFa,sx_rho_antisym,sy_rho_antisym,sz_rho_antisym)  ! Coulomb response to antisymmetric ph densities

    dF_pp_plus  = calc_perturbed_potentials_oneoff(R,dR_pp_plus) ! perturbation of pp potentials with kappa
    dF_pp_minus = calc_perturbed_potentials_oneoff(R,dR_pp_minus)! perturbation of pp potentials with kappa^*

    !call print_maxval('F_I_I', dFs%F_I_I, dFa%F_I_I)
    !call print_maxval('F_Nm_Nm', dFs%F_Nm_Nm, dFa%F_Nm_Nm)

    !call print_maxval('F_I_Sx', dFs%F_I_S(:,1,:), dFa%F_I_S(:,1,:))
    !call print_maxval('F_I_Sy', dFs%F_I_S(:,2,:), dFa%F_I_S(:,2,:))
    !call print_maxval('F_I_Sz', dFs%F_I_S(:,3,:), dFa%F_I_S(:,3,:))

    !call print_maxval('G_I_Nx', dFs%G_I_N(:,1,:), dFa%G_I_N(:,1,:))
    !call print_maxval('G_I_Ny', dFs%G_I_N(:,2,:), dFa%G_I_N(:,2,:))
    !call print_maxval('G_I_Nz', dFs%G_I_N(:,3,:), dFa%G_I_N(:,3,:))

    !call print_maxval('G_I_NS', dFs%G_I_NS, dFa%G_I_NS)
    
    !print *, 'MAXVAL dRs Re', maxval(ABS(DBLE(dRs%chargedensity)))
    !print *, 'MAXVAL dRs Im', maxval(ABS(AIMAG(dRs%chargedensity)))
    !print *, 'MAXVAL dRa Re', maxval(ABS(DBLE(dRa%chargedensity)))
    !print *, 'MAXVAL dRa Im', maxval(ABS(AIMAG(dRa%chargedensity)))
    
    !print *, 'MAXVAL Fc Re', maxval(ABS(DBLE(dFs%CoulombPotential(:,:,:))))
    !print *, 'MAXVAL Fc Im', maxval(ABS(AIMAG(dFs%CoulombPotential(:,:,:))))
    !print *, 'MAXVAL Fa Re', maxval(ABS(DBLE(dFa%CoulombPotential(:,:,:))))
    !print *, 'MAXVAL Fa Im', maxval(ABS(AIMAG(dFa%CoulombPotential(:,:,:))))

  end subroutine calc_perturbed_potentials


  module function calc_perturbed_potentials_oneoff(R,R_pert) result (F)
    !---------------------------------------------------------------------------
    ! Calculate the perturbed potential vector (F) around a set of mean-field
    ! densities (R) affected by a perturbation (R_pert).
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   R    : density-vector containing the mean-field densities
    !   R_pert: density-vector containig the perturbation to the densities
    ! Output:
    !   F: potential-vector containing the linearised response of the mean-field
    !      potentials/
    !---------------------------------------------------------------------------
    use Coulombmod       , only : solve_coulomb
    use pairing_strengths, only : vmicro, vmicro_stored

    type (DensityVector), intent(in) :: R, R_pert
    type (PotentialVector)           :: F

    vmicro_stored = .false.

$CALCPOTENTIALS_PERTURBED
    
  end function calc_perturbed_potentials_oneoff
#endif 

  module subroutine combine_potentials(F)
    !----------------------------------------------------------------------------
    ! This routine makes sure that the potential vector F is ready to be used
    ! in the single-particle hamiltonian by summing various contributions that
    ! get calculated separately.
    !
    ! 1. add Coulomb potentials to F_I_I
    ! 2. add contribution from constraints to F_I_I (and possibly others)
    ! 3. add cranking contribution to F_I_S and G_I_N
    !
    ! Input:
    !   F : potentialvector, fresh from the calculation routines after possibly
    !       passing through preconditioning and mixing
    ! Output:
    !   F : potentialvector, ready to be used in the action of the single-particle
    !       hamiltonian
    !----------------------------------------------------------------------------
    type(PotentialVector), intent(inout) :: F

    ! Add the coulomb potential to F_I_I
    call add_coulomb_potential(F)
    ! Add the contribution of the multipole constraints
    call add_constraint_potential(F)
    ! cranking constraints add stuff to the spin and current potentials
    call add_cranking_potentials(F)

  end subroutine combine_potentials


  module subroutine add_constraint_potential(F)
    !----------------------------------------------------------------------------
    ! Add the contribution of constraints on the mass electric multipole moments
    ! to a potentialvector.
    !
    ! Input:
    !   F :  potentialvector with F_I_I
    ! Output:
    !   F : where the new F_I_I = old F_I_I + constraint contribution
    !-----------------------------------------------------------------------------
    use moments

    type(PotentialVector), intent(inout) :: F

    ! Add the contribution from the constraints on the multipole moments
    ! This line also sets Constraint_I_I globally!
    Constraint_I_I = constraints_sph_elmult(.false.)
    F%F_I_I(:,1:2) = F%F_I_I(:,1:2) + Constraint_I_I(:,1:2)

    ! We added stuff to the proton and neutron potentials, we should be
    ! consistent with the isospin 0 and 1 potentials
    F%F_I_I(:,3) = F%F_I_I(:,1) + F%F_I_I(:,2)
    F%F_I_I(:,4) = F%F_I_I(:,1) - F%F_I_I(:,2)

  end subroutine add_constraint_potential


  module subroutine add_coulomb_potential(F)
    !----------------------------------------------------------------------------
    ! Add the Coulomb contribution to the potential F_I_I.
    ! The index juggling is ugly, but necessary, because the Coulomb potential is
    ! defined on a larger mesh than F_I_I.
    !
    ! Input:
    !   F :  potentialvector with F_I_I (= associated with the Skyrme EDF)
    !        and calculated Coulomb potentials
    ! Output:
    !   F : where the new F_I_I = old F_I_I + coulomb potentials
    !-----------------------------------------------------------------------------
    type(PotentialVector), intent(inout) :: F

    F%F_I_I = F%F_I_I + transfer_coulomb_mesh(F, .true.) &
    &                 + transfer_coulomb_mesh(F, .false.)
  
  end subroutine add_coulomb_potential


  module subroutine add_cranking_potentials(F)
    !-----------------------------------------------------------------------------
    ! Add the cranking contributions to the
    !   F_I_S (spin potential)
    !   G_I_N (current potential)
    !
    ! Input:
    !    F : potentialvector
    ! Output:
    !    F : potentialvector, with modified F_I_S and G_I_N.
    !-----------------------------------------------------------------------------
    type(PotentialVector), intent(inout) :: F
$TR real(KIND=dp)                        :: trash
$TR trash = F%F_I_I(1,1) ! to stop compiler complaints when time-reversal is conserved
!
$TAUSCALAR $NTR F%F_I_S = F%F_I_S + crank_spin_potential()
$TAUSCALAR $NTR F%G_I_N = F%G_I_N + crank_current_potential()

$TAUTENSOR $NTR F%F_I_S = F%F_I_S + crank_spin_potential()
$TAUTENSOR $NTR F%G_I_N = F%G_I_N + crank_current_potential()

  end subroutine add_cranking_potentials

#if( $FAM == 0 )

  module function precondition_potentials(F_in, F_out) result(F)
    !---------------------------------------------------------------------------
    ! Precondition the change in potentials from one iteration to the next. 
    ! All code generated by Hephaestos. 
    !
    ! Input:
    !   F_in : input potentials for an SCF iteration
    !   F_out: output potentials for an SCF iteration
    !
    ! Output:
    !   F    : set of preconditioned potentials. For every 
    !            F_a = F_in,a + P_a [F_out,a - F_in,a]
    !          where P_a is the preconditioning operator for µ
    !          mean-field potential a.
    !------------------------------------------ ---------------------------------
    
    type(PotentialVector), intent(in) :: F_in, F_out
    type(PotentialVector)             :: F
    real(KIND=dp), allocatable        :: update(:,:)
    real(KIND=dp)                     :: coul_in(mv,4), coul_out(mv,4)

    call start_timer(T_potentials)
    call start_timer(T_pot_precon)

    !---------------------------------------------------------------------------
    ! Everything which is not specifically preconditioned below just gets 
    ! explicitly copied from the new values.
    F = F_out
    !---------------------------------------------------------------------------
    ! Preconditioning for the Skyrme potentials: we safeguard against 
    ! short wavelength modes with a high-pass filter.
$POTENTIALPRECON

    !---------------------------------------------------------------------------
    ! Preconditioning for the Coulomb potential
    ! TODO: also apply the preconditioners to the exchange potential!
    !---------------------------------------------------------------------------
    ! 1. get coulomb potentials on a typical mesh
    coul_out = transfer_coulomb_mesh(F_out,.false.) 
    coul_in  = transfer_coulomb_mesh(F_in,.false.)
    ! 2. calculate the difference
    update = coul_out - coul_in

    if(kerker_k0 .gt. 0.0d0) then    
      ! 3a. precondition with Kerker preconditioner 
      update = KerkerPreconditionPotential(update,mixstepsize,kerker_k0,sx_rho,sy_rho,sz_rho)
    endif

    ! 3b. precondition with the same preconditioner as F_I_I 
    ! OPTIONAL CHANGE FOR later
    ! update = PreconditionPotential(update,-preconfactor,1.0_dp, sx_rho,sy_rho,sz_rho)

    ! 4. save the result
    update = coul_in + update
    call set_coul(update(:,1:2), F)
    !---------------------------------------------------------------------------

    call stop_timer(T_pot_precon)
    call stop_timer(T_potentials)
 
  end function precondition_potentials
#endif


  module subroutine clean_potentials()
    !---------------------------------------------------------------------------
    ! Clean up the allocated potentials for multiple runs.
    !---------------------------------------------------------------------------

  end subroutine clean_potentials
  

module subroutine save_potential_history(F_in, F_out)
    !---------------------------------------------------------------------------
    ! Update the history of the mean-field densities, before calculating an
    ! updated version.
    !---------------------------------------------------------------------------
    integer :: i
    type(PotentialVector), intent(in) :: F_in, F_out

    if(memory .eq. 0) return
  
    if(.not.allocated(Potential_updates)) allocate(Potential_updates(memory))
    if(.not.allocated(Potential_iterates)) allocate(Potential_iterates(memory))

    do i=1,memory-1
      Potential_updates(memory-i+1)  = Potential_updates(memory-i)
      Potential_iterates(memory-i+1) = Potential_iterates(memory-i)
    enddo
    Potential_iterates(1)  = F_in
    Potential_updates(1)   = F_out + (-1.0d0) * F_in
!  print *, 'SAVING', maxval(abs(Potential_updates(1)%F_I_I))
!  print *,  maxval(abs(F_in%F_I_I))
!  print *,  maxval(abs(F_out%F_I_I))
end subroutine save_potential_history


end submodule functional_potentials
