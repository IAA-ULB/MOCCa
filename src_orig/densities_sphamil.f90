submodule (densities) densities_sphamil
!==============================================================================
!
! Implementation submodule of module densities: densities_sphamil
! The interfaces live in module densities (src_orig/densities.f90).
!
!==============================================================================


 implicit none

contains

module function calc_sphamil_me(denpsi, dendpsi, denddpsi, Fs, Fa, onthefly) result(sphamil_me)
    !------------------------------------------------------------------------------------
    ! This function calculates the single-particle matrix elements of the perturbed
    ! single-particle hamiltonian in total, i.e. including both symmetric and antisymmetric
    ! parts.
    !
    ! In general we have
    !
    ! \delta \langle a | h | b \rangle ~  \int d^3 r \delta F^{A,B} f^{A,B}_{ab}
    !
    ! where both the perturbed potentials and the spwf-factor get split into two different
    ! pieces with different symmetry properties
    !
    ! \delta F^{A,B} =  \delta F^{A,B}_{sym} + \delta F^{A,B}_{antisym}
    ! f^{A,B}_{ab}   =  \Re f^{A,B}_{ab} + i \Im f^{AB}_{ab}
    !
    ! In EV8-like calculations, only
    !   \int d^3r \delta F^{A,B}_{sym}     \Re f^{A,B}_{ab}
    !   \int d^3r \delta F^{A,B}_{antisym} \Im f^{A,B}_{ab}
    !
    ! will not vanish for symmetry reasons; this routine only implements these two
    ! contributions and will thus have to be rewritten to accomodate more general cases.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    ! Input:
    ! -------
    !   denpsi  : set of single-particle wavefunctions
    !   dendpsi : their first order derivatives
    !   denddpsi: their second order derivatives
    !   Fs      : potential vector containing the     SYMMETRIC linearised response
    !             of the mean-field potentials
    !   Fa      : potential vector containing the ANTISYMMETRIC linearised response
    !             of the mean-field potentials
    !   onthefly: [NOT IMPLEMENTED YET ]
    !
    ! Output:
    ! -------
    !  sphamil_me : matrix elements of the single-particle hamiltonian.
    !------------------------------------------------------------------------------------
    real(KIND=dp), intent(in)         :: denpsi(:,:,:), dendpsi(:,:,:,:), denddpsi(:,:,:,:)
    logical, intent(in)               :: onthefly
    type(PotentialVector), intent(in) :: Fs, Fa
    complex(KIND=dp), allocatable     :: sphamil_me(:,:), sp_sym(:,:), sp_asym(:,:)
    integer                           :: B, N, i, si

    call start_timer(T_spme_perturbed)
    sp_sym = calc_sphamil_me_sym    ( denpsi, dendpsi, denddpsi, Fs,  onthefly)
    sp_asym= calc_sphamil_me_antisym( denpsi, dendpsi, denddpsi, Fa,  onthefly)
    sphamil_me = sp_sym + sp_asym
    call stop_timer(T_spme_perturbed)

end function calc_sphamil_me

module function calc_sphamil_me_sym( denpsi, dendpsi, denddpsi, F,  onthefly) result(sphamil_me)
    !---------------------------------------------------------------------------------------
    ! This function calculates the single-particle matrix elements of the symmetric part
    ! of the single-particle hamiltonian as defined by a potentialvector dF, which should
    ! contain the symmetric part of a set of perturbed potentials.
    !
    ! Note:
    ! - this routine has no way of checking that you fed it a dF with correct properties.
    ! - this routine does not assume hermeticity of the matrix elements
    ! - the spwfs on input are named  "den[d/dd]psi" in order to have less changes
    !   in Hephaestos; these are not the pointers defined on top in this module.
    ! - this routine ASSUMES that F_I_I is the FULL potential, i.e. that
    !   combine_potentials has been called on F before using this!
    !
    ! TODO:
    !  - rename wavefunctions for clarity -> requires Hephaestos change
    !  - develop MPI parallelism
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    ! Input:
    ! -------
    !   denpsi  : set of single-particle wavefunctions
    !   dendpsi : their first order derivatives
    !   denddpsi: their second order derivatives
    !    F      : potential vector containing the SYMMETRIC linearised response
    !             of the mean-field potentials
    !   onthefly: [NOT IMPLEMENTED YET ]
    !
    ! Output:
    ! -------
    !  sphamil_me : matrix elements of the single-particle hamiltonian.
    !------------------------------------------------------------------------------------------
    real(KIND=dp), intent(in)         :: denpsi(:,:,:), dendpsi(:,:,:,:), denddpsi(:,:,:,:)
    logical, intent(in)               :: onthefly
    type(PotentialVector), intent(in) :: F

    complex(KIND=dp), allocatable     :: sphamil_me(:,:)
    integer                           :: it, B, si, N, wave_i, wave_j, i

$SPWF_DECLARATION

    call start_timer(T_spme_perturbed_sym)

    ! initialize
    allocate(sphamil_me(nwt,nwt)) ; sphamil_me = 0.0d0

!$OMP PARALLEL PRIVATE(si, B, N, it, i,             &
!$OMP                  wave_i, wave_j,              &
$OMP_DENSITY_VARS
!$OMP                 )                             &
!$OMP          SHARED ( HFBlocks, F, mv, dv,        &
!$OMP                   sphamil_me,                 &
!$OMP                   denpsi, dendpsi, denddpsi) DEFAULT(PRIVATE)
    do B=1,8
      N = HFBlocks(B); if(N.eq.0) cycle
      si = sum(HFBlocks(1:B-1)) ! offset has to be explicitly calculated for OpenMP
      !---------------------------------------------------------------------------
      ! Determine the isospin index
      if(B.ge.5) then
        it = 2
      else
        it = 1
      endif
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
      do wave_j=si+1,si+N  ! Note: no assumption of hermeticity here!
        do wave_i=si+1,si+N
          do i=1,mv
$EXPRESSION_SPH_SYM
          enddo
          sphamil_me(wave_j, wave_i) = sphamil_me(wave_j, wave_i) * dv
        enddo
      enddo
      !$OMP END DO
    enddo
!$OMP END PARALLEL
    call stop_timer(T_spme_perturbed_sym)

  end function calc_sphamil_me_sym

module function calc_sphamil_me_antisym( denpsi, dendpsi, denddpsi, F,  onthefly) result(sphamil_me)
    !---------------------------------------------------------------------------------------
    ! This function calculates the single-particle matrix elements of the ANTIsymmetric part
    ! of the single-particle hamiltonian as defined by a potentialvector dF, which should
    ! contain the ANTIsymmetric part of a set of perturbed potentials.
    !
    ! Note:
    ! - this routine has no way of checking that you fed it a dF with correct properties.
    ! - this routine does not assume hermeticity of the matrix elements
    ! - the spwfs on input are named  "den[d/dd]psi" in order to have less changes
    !   in Hephaestos; these are not the pointers defined on top in this module.
    ! - this routine ASSUMES that F_I_I is the FULL potential, i.e. that
    !   combine_potentials has been called on F before using this!
    !
    ! TODO:
    !  - rename wavefunctions for clarity -> requires Hephaestos change
    !  - develop MPI parallelism
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    ! Input:
    ! -------
    !   denpsi  : set of single-particle wavefunctions
    !   dendpsi : their first order derivatives
    !   denddpsi: their second order derivatives
    !    F      : potential vector containing the ANTISYMMETRIC linearised response
    !             of the mean-field potentials
    !   onthefly: [NOT IMPLEMENTED YET ]
    !
    ! Output:
    ! -------
    !  sphamil_me : matrix elements of the single-particle hamiltonian.
    !------------------------------------------------------------------------------------------
    real(KIND=dp), intent(in)         :: denpsi(:,:,:), dendpsi(:,:,:,:), denddpsi(:,:,:,:)
    logical, intent(in)               :: onthefly
    type(PotentialVector), intent(in) :: F

    complex(KIND=dp), allocatable     :: sphamil_me(:,:)
    integer                           :: it, B, si, N, wave_i, wave_j, i

$SPWF_DECLARATION

    call start_timer(T_spme_perturbed_asym)

    ! initialize
    allocate(sphamil_me(nwt,nwt)) ; sphamil_me = 0.0d0

!$OMP PARALLEL PRIVATE(si, B, N, it, i,             &
!$OMP                  wave_i, wave_j,              &
$OMP_DENSITY_VARS
!$OMP                 )                             &
!$OMP          SHARED ( HFBlocks, F, mv, dv,        &
!$OMP                   sphamil_me,                 &
!$OMP                   denpsi, dendpsi, denddpsi) DEFAULT(PRIVATE)
    do B=1,8
      N = HFBlocks(B)
      N = HFBlocks(B); if(N.eq.0) cycle
      si = sum(HFBlocks(1:B-1)) ! offset has to be explicitly calculated for OpenMP
      !---------------------------------------------------------------------------
      ! Determine the isospin index
      if(B.ge.5) then
        it = 2
      else
        it = 1
      endif
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
      do wave_j=si+1,si+N  ! Note: no assumption of hermeticity here!
        do wave_i=si+1,si+N
          do i=1,mv
$EXPRESSION_SPH_ANTISYM
          enddo
          sphamil_me(wave_j, wave_i) = sphamil_me(wave_j, wave_i) * dv
        enddo
      enddo
!$OMP END DO
    enddo
!$OMP END PARALLEL
    call stop_timer(T_spme_perturbed_asym)

  end function calc_sphamil_me_antisym

  module function calc_delta_me( denpsi, dendpsi, denddpsi, F,  onthefly) result(delta_me)
    !------------------------------------------------------------------------
    ! Calculate the matrix elements of \Delta in the Hartree-Fock basis.
    !     
    ! Attention: the pairing cutoffs figure in this function and those 
    !            are calculated in the Hartree-Fock basis. Even if this function 
    !            takes denpsi, dendpsi, denddpsi as input, these should be
    !            the Hartree-Fock wavefunctions!
    !
    ! Input:
    ! -------
    !   denpsi  : set of single-particle wavefunctions
    !   dendpsi : their first order derivatives
    !   denddpsi: their second order derivatives
    !    F      : potential vector containing the linearised pairing potentials
    !   onthefly: [NOT IMPLEMENTED YET ]
    !
    !
    ! Output:
    ! -------
    !  delta_me : matrix elements of the pairing tensor Delta.
    !
    !
    ! TODO: 
    ! - implement MPI parallelisation
    !
    !------------------------------------------------------------------------
    real(KIND=dp), intent(in)         :: denpsi(:,:,:), dendpsi(:,:,:,:), denddpsi(:,:,:,:)
    logical, intent(in)               :: onthefly
    type(PotentialVector), intent(in) :: F
    complex(KIND=dp), allocatable     :: delta_me(:,:)
    integer                           :: it, B, si, N, N2, wave_i, wave_j, i, T
$SPWF_DECLARATION
    call start_timer(T_spme_perturbed_pp)

    ! initialize
    allocate(delta_me(nwt,nwt)) ; delta_me = 0.0d0

!$OMP PARALLEL PRIVATE(si, B, N, it, i,             &
!$OMP                  wave_i, wave_j,              &
$OMP_DENSITY_VARS
!$OMP                 )                             &
!$OMP          SHARED ( HFBlocks, F, mv, dv,        &
!$OMP                   delta_me, Pcutoffs,         &
!$OMP                   denpsi, dendpsi, denddpsi) DEFAULT(PRIVATE)
    do B=1,8,2
      N = HFBlocks(B) ;  if (N.eq.0) cycle
      N2= HFBlocks(B+1)
      T = N+N2

      si = sum(HFBlocks(1:B-1)) ! OMP requires explicit offset calculation
      !---------------------------------------------------------------------------
      ! Determine the isospin index
      if(B.ge.5) then
        it = 2
      else
        it = 1
      endif
!$OMP DO SCHEDULE(STATIC)
!     Technical note: unfortunately, these two loops cannot be OpenMP
!     COLLAPSED as they form a non-rectangular iteration space. This is only
!     accepted by SOME compilers.
      do wave_i=si+1,si+N                       ! local index of the spwf
$TR          do wave_j=wave_i,si+N              ! symmetry-reduced
$NTR          do wave_j=si+N+1,si+N+N2      
          do i=1,mv
$EXPRESSION_DELTA_PP
          enddo
          ! The minus sign is because Hephaestos generates the expression for 
          ! 
          ! \tilde \rho_ji = \sum_{\sigma} \sigma psi_j(r',\sigma) \psi_i(r,-\sigma)
          !   
          ! whereas Delta is proportional to 
          !
          !  \Delta_ji \sim \tilde \rho_ij
          !
          ! This should be corrected in Hephaestos, but it is much harder than including this minus sign.
          delta_me(wave_j, wave_i) =  - delta_me(wave_j, wave_i) * dv * Pcutoffs(wave_i) * PCutoffs(wave_j)
          ! Sign to be clarified with better documentation of Hephaestos
$TR       delta_me(wave_j, wave_i) =  - delta_me(wave_j, wave_i) 
          ! Delta is globally antisymmetric in the case of time-reversal symmetry, but we 
          !  represent only half of the matrix explicitly!
$TR       delta_me(wave_i, wave_j) =  delta_me(wave_j, wave_i) 
$NTR      delta_me(wave_i, wave_j) = -delta_me(wave_j, wave_i) 
        enddo
      enddo
!$OMP END DO
   enddo
!$OMP END PARALLEL
    call stop_timer(T_spme_perturbed_pp)

  end function calc_delta_me

end submodule densities_sphamil
