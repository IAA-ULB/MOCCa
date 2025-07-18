module momentsofinertia
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
 ! TODO: 
 !    (1) improve on the efficiency of calcJ2andBelyaev routines through 
 !        switching the looping structure. The inner (j) loop currently has
 !        more communication/calculation to do and should be exchanged with 
 !        the outer (i) loop.
 !    (2) improve on the documentation of these routines; develop notes.
 !    (3) improve on the "transfer everything" MPI_ALLREDUCE calls, these can
 !        likely be exchanged with MPI_BCASTs per symmetry block. 
 !------------------------------------------------------------------------------

  use geninfo
  use densities
  use parameterization

  implicit none

  !-------------------------------------------------------------------------------
  ! Belyaev moments of inertia in units of hbar^2/MeV
  !
  ! The full definition of Belyaev moment of inertia around the Cartesian axes,
  !
  !           I_mm = - partial_{\omega} < psi' | J_m | psi' >
  !
  !        where
  !        a) the derivative should be evaluated at a reference frequency
  !            \omega
  !        b) psi' is a second order perturbation to the manybody state when
  !           the rotational frequency changes, i.e.
  !
  !               H - \omega J_m => H - \omega J_m  - do J_m
  !
  !
  !           hence
  !
  !          | psi'> - |psi_0 > =
  !             do               C
  !          -  --- sum_ab -------------  b^{\dagger}_a b^{\dagger}_a |psi_0 > (1)
  !              2           E_a + E_b
  !
  !           when T != 0 or in the presence of quasiparticle excitations
  !           there are also terms involving annihilation operators
  !
  !          with E_a and E_b the quasiparticle energies and
  !
  !           C =  < psi_0 | J_m  b^{\dagger}_a b^{\dagger}_a |psi_0 >
  !
  !  The array Belyaev corresponds to the calculation of Eq. (1), but this kind
  !  of perturbation theory is ill-defined when the reference many-body state
  !  is ane excitation, i.e. when one (or more) of the quasiparticle energies
  !  are negative.
  !
  !  The array Bely_coll "solves" this problem more or less by removing the
  !  blocked quasiparticles from Eq. (1); i.e. it tries to get the moment of inertia
  !  of the "even-even core".
  !
  real(KIND=dp) :: Belyaev(3,3), Bely_coll(3,3)
  !                        | |> Isospin index
  !                        |-> Cartesian index
  !-------------------------------------------------------------------------------
  ! Rigid rotor moment of inertia in units of hbar^2/MeV
  !
  !   I^rigid_x = m int d^3r rho * (y^2 + z^2)
  !   I^rigid_y = m int d^3r rho * (x^2 + z^2)
  !   I^rigid_z = m int d^3r rho * (x^2 + y^2)
  real(KIND=dp) :: Rigid_MOI(3,3)
  !                      | |> Isospin index
  !                      |-> Cartesian index
  !-------------------------------------------------------------------------------
  ! Dispersion of the angular momentum squared, i.e.
  !
  !   \Delta J_{\mu}  = < J^2_{\mu} > - < J_{\mu} >^2    (\mu =x,y,z)
  !
  ! Evaluating things in second quantisation, this can be written as
  !
  ! Equation (2):
  ! \Delta J_{\mu} =
  !   \sum_{ab} <a|J_{\mu}|b><b|J_{\mu}|a> [ \rho_aa ( 1  - \rho_bb)]
  ! - \sum_{ab} <a|J_{\mu}|b><\bar{a}|J_{\mu}|\bar{b} > \kappa^*_{a\bar{a}} \kappa_{b\bar{b}}
  !
  ! where the indices ab label single-particle states in the canonical basis and
  ! \bar{a}, \bar{b} label their canonical partners.
  !
  ! The first term on the rhs of (2) can be rewritten as
  !
  ! Equation (3):
  !  \sum_{ab} < a | J_{\mu}^2 | a > \rho_aa
  !      - \sum_{ab} < a | J_{\mu} | b > < b | J_{\mu} | a > [ \rho_aa \rho_bb]
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  ! This quantity gets calculated in several different ways
  !
  ! * J2
  !   = the ordinary, correct way: direct calculation of Equation (2), using
  !      Equation (3) for the first term.
  !
  ! * J2_pairing_cut
  !   = include the pairing cutoff in the definition all matrix
  !     elements of J_{\mu} and then use Eq.(2) WITHOUT using Eq.(3).
  !
  ! * J2_collective
  !   = a "collective" expectation value that is identical to J2_pairing_cut, but
  !     that removes all blocked states from the summations over a and b.
  !
  real(KIND=dp) ::  J2(3,3), J2_pairing_cut(3,3), J2_coll(3,3)
  !                    | |> Isospin index
  !                    |-> Cartesian index
  !-------------------------------------------------------------------------------
  ! Procedure pointer to perform the calculation of all quantities in this
  ! module; defined differently depending on whether we are doing HF,HF+BCS or HFB.
  procedure(calcJ2andBelyaev_HF), pointer :: calcJ2andBelyaev 

contains

  subroutine setBelyaevProcedure()
    !-----------------------------------------------------------------------------
    ! Set the procedure pointer calcJ2andBelyaev depending on the way
    ! we consider pairing correlations.
    !-----------------------------------------------------------------------------
    select case(Pairingtype)
    case(0)
      calcJ2andBelyaev => calcJ2andBelyaev_HF 
    case(1)
      calcJ2andBelyaev => calcJ2andBelyaev_BCS
    case(2)
      calcJ2andBelyaev => calcJ2andBelyaev_HFB
    end select
  end subroutine

  subroutine calcrigid(R)
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
    type(DensityVector), intent(in), target :: R
    real(KIND=dp), allocatable              :: rho(:)
    real(KIND=dp)                           :: xs(2), ys(2), zs(2)
    integer                                 :: it, i,j,k

    Rigid_MOI = 0

    xs = 0 ; ys = 0 ; zs = 0

    do it=1,2
      rho = DBLE(R%D_I_I(:,it)) ! taking the real part to ensure compatibility 
                                ! for complex densities
      do k=1,nz
        do j=1,ny
          do i=1,nx
            xs(it) = xs(it) + meshx_shifted(i)**2 * rho(meshindex(i, j, k))
            ys(it) = ys(it) + meshy_shifted(j)**2 * rho(meshindex(i, j, k))
            zs(it) = zs(it) + meshz_shifted(k)**2 * rho(meshindex(i, j, k))
          enddo
        enddo
      enddo
    enddo

    xs = xs * dv ; ys = ys *dv ; zs = zs * dv

    Rigid_MOI(1,1:2) = nucleonmass * (ys + zs)
    Rigid_MOI(2,1:2) = nucleonmass * (xs + zs)
    Rigid_MOI(3,1:2) = nucleonmass * (xs + ys)

    Rigid_MOI(:,3) = sum(Rigid_MOI(:,1:2),2)

    ! Converting to the correct units
    Rigid_MOI = Rigid_MOI/(hbarclum**2)

  end subroutine calcrigid

  pure real(KIND=dp) function rotcut(eps,it) result(cut)
      !-------------------------------------------------------------------------
      ! Cutoff for use in the calculation of < J^2 > and the Belyaev moment
      ! of inertia. Similar in shape to the pairing cutoff, but sharper.
      !
      ! eps : single-particle energy in the HF-basis
      ! it  : isospin index
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in) :: eps
      integer, intent(in)       :: it
      real(KIND=dp)             :: Up

      Up   =   2*  (eps - FermiEnergy(it) - RotCutWindow(it))/RotCutMu(it)
      cut  = sqrt(sqrt(1.0_dp/(1.0_dp + exp(Up))))
  end function rotcut

  subroutine calcJ2andBelyaev_HF
    !---------------------------------------------------------------------------
    ! Calculate the expectation value of J^2_mu in the many-body state and the
    ! Belyaev moment of inertia in the case of a Hartree-Fock calculation.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Note for the future: I suspect this routine can be sped up by switching 
    ! the loops over i <-> j, since the latter needs more communication.
    !
    ! In case of time-reversal breaking, I suspect the loops for the Z
    ! and (X,Y) matrix elements can be combined. 
    !---------------------------------------------------------------------------
    integer       :: i,j, b, it, ii, jj, si, N, N2, locali, localj
    integer       :: calc_rank, ranki, rankj, designated_rank(8)
    real(KIND=dp) :: ME(3), fi, fj, dfde
    real(KIND=dp) :: psi_i(mv,4), psi_j(mv,4), der_psi_j(mv,3,4)
#if(USE_MPI>0)
    integer :: mpi_err
#endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    designated_rank = -1 
    do B=1,8
      ! First, we designate a rank to do the calculations
      !    = the first rank storing spwfs in a symmetry block
      N = HFBlocks_global(B) ;  if(N.eq.0) cycle
      designated_rank(B) = rank_map(sum(HFBlocks_global(1:B-1))+1)
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    J2_pairing_cut = 0  ;  Belyaev = 0
 
    si = 0  
    do B = 1,8,2                     ! <----- loop over global indices of spwfs
      N = HFBlocks_global(B) ; if(N.eq.0) cycle 
      N2= HFBlocks_global(B+1)

      it = 1
      if(B.gt.4) it = 2

      calc_rank = designated_rank(B) ! <---- this rank will do the integrations
                                     !       of the spwf matrix elements
      do j=1, N
        jj    = si + j            ! global index of the spwf
        localj= spwf_inverse(jj)  ! local index of the spwf
        rankj = rank_map(jj)      ! MPI rank storing the spwf

#if(USE_MPI>0)
        call Transfer_psi(psi_j, localj, 'HF', rankj, calc_rank)
        call Transfer_derpsi_complete(der_psi_j,localj,'HF',rankj, calc_rank)
#else
        call Transfer_psi(psi_j, localj, 'HF')
        call Transfer_derpsi_complete(der_psi_j,localj,'HF')
#endif
        do i=1,N
          ii    = si + i            ! global index of the spwf
          locali= spwf_inverse(ii)  ! local index of the spwf
          ranki = rank_map(ii)      ! MPI rank storing the spwf

#if(USE_MPI>0)
          call Transfer_psi(psi_i, locali, 'HF', ranki, calc_rank)
#else
          call Transfer_psi(psi_i, locali, 'HF')
#endif

          if(MPI_RANK.eq.calc_rank) then
            ! Only one rank should do the integration of the matrix elements
$TR         fi = rho_can(ii)/2.0 ; fj = rho_can(jj)/2.
$NTR        fi = rho_can(ii)     ; fj = rho_can(jj)
            ! |< k | j_z |  l >|^2 
            ME(3)= angmom_z_real(psi_i,psi_j,der_psi_j)**2  
            J2_pairing_cut(3,it) = J2_pairing_cut(3,it) +  ME(3) * fi*(1-fj)

$TR         ! |< k | j_x | -l >|^2  
$TR         ME(1)= angmom_xt_real(psi_i,psi_j,der_psi_j)**2 
$TR         ! |< k | j_y | -l >|^2 
$TR         ME(2)= angmom_yt_imag(psi_i,psi_j,der_psi_j)**2
$TR         J2_pairing_cut(1:2,it) = J2_pairing_cut(1:2,it) +  ME(1:2) * fi*(1-fj)

            if(inversetemp.eq.-1) then
              dfdE = fj - fi
              if(abs(dfdE).gt.0) then
                  dfdE = dfdE/(spenergies(ii) - spenergies(jj))
              endif
            else
              if(abs(spenergies(ii) - spenergies(jj)).gt.1d-8) then
                dfdE = (fj - fi)/(spenergies(ii) - spenergies(jj))
              else                
                dfdE = fi**2 * inversetemp                                     &
                &            * exp(inversetemp*(spenergies(ii)-FermiEnergy(it)))
              endif
            endif
            Belyaev(3,it)   = Belyaev(3,it)   + ME(3)   * dfde
$TR         Belyaev(1:2,it) = Belyaev(1:2,it) + ME(1:2) * dfde
          endif
        enddo
      enddo

$NTR  do j=1,N2
$NTR    jj = si + N + j           ! global index of the spwf
$NTR    localj= spwf_inverse(jj)  ! local index of the spwf
$NTR    rankj = rank_map(jj)      ! MPI rank storing the spwf
$NTR
#if(USE_MPI>0)
$NTR    call Transfer_psi(psi_j, localj, 'HF', rankj, calc_rank)
$NTR    call Transfer_derpsi_complete(der_psi_j,localj,'HF',rankj, calc_rank)
#else
$NTR    call Transfer_psi(psi_j, localj, 'HF')
$NTR    call Transfer_derpsi_complete(der_psi_j,localj,'HF')
#endif
$NTR    do i=1,N2
$NTR      ii    = si + i            ! global index of the spwf
$NTR      locali= spwf_inverse(ii)  ! local index of the spwf
$NTR      ranki = rank_map(ii)      ! MPI rank storing the spwf
$NTR
$NTR      J2_pairing_cut(3,it) = J2_pairing_cut(3,it) +  ME(3) * fi*(1-fj)
$NTR
#if(USE_MPI>0)
$NTR      call Transfer_psi(psi_i, locali, 'HF', ranki, calc_rank)
#else
$NTR      call Transfer_psi(psi_i, locali, 'HF')
#endif
$NTR       if(MPI_RANK.eq.calc_rank) then
$NTR         ME(3)= angmom_z_real(psi_i,psi_j,der_psi_j)**2
$NTR         fi = rho_can(ii)     ; fj = rho_can(jj)
$NTR         J2_pairing_cut(3,it) = J2_pairing_cut(3,it) +  ME(3) * fi*(1-fj)
$NTR
$NTR         if(inversetemp.eq.-1) then
$NTR           dfdE = fj - fi
$NTR           if(abs(dfdE).gt.0) then
$NTR              dfdE = dfdE/(spenergies(ii) - spenergies(jj))
$NTR           endif
$NTR         else
$NTR           if(abs(spenergies(ii) - spenergies(jj)).gt.1d-8) then
$NTR             dfdE = (fj - fi)/(spenergies(ii) - spenergies(jj))
$NTR           else                
$NTR             dfdE = fi**2 * inversetemp                                   &
$NTR             &          * exp(inversetemp*(spenergies(ii)-FermiEnergy(it)))
$NTR           endif
$NTR         endif
$NTR         Belyaev(3,it) = Belyaev(3,it) + ME(3) * dfde
$NTR       endif
$NTR     enddo
$NTR   enddo

$NTR   do j=1,N2
$NTR     jj = si + N + j           ! global index of the spwf
$NTR     localj= spwf_inverse(jj)  ! local index of the spwf
$NTR     rankj = rank_map(jj)      ! MPI rank storing the spwf
$NTR
#if(USE_MPI>0)
$NTR     call Transfer_psi(psi_j, localj, 'HF', rankj, calc_rank)
$NTR     call Transfer_derpsi_complete(der_psi_j,localj,'HF',rankj, calc_rank)
#else
$NTR     call Transfer_psi(psi_j, localj, 'HF')
$NTR     call Transfer_derpsi_complete(der_psi_j,localj,'HF')
#endif
$NTR     do i=1,N
$NTR      ii    = si + i            ! global index of the spwf
$NTR      locali= spwf_inverse(ii)  ! local index of the spwf
$NTR      ranki = rank_map(ii)      ! MPI rank storing the spwf
#if(USE_MPI>0)
$NTR      call Transfer_psi(psi_i, locali, 'HF', ranki, calc_rank)
#else
$NTR      call Transfer_psi(psi_i, locali, 'HF')
#endif
$NTR       if(MPI_RANK.eq.calc_rank) then
$NTR         fi = rho_can(ii)     ; fj = rho_can(jj)
$NTR         ME(1)= angmom_x_real( psi_i,psi_j,der_psi_j)**2  
$NTR         ME(2)= angmom_y_imag( psi_i,psi_j,der_psi_j)**2  
$NTR         J2_pairing_cut(1:2,it) = J2_pairing_cut(1:2,it) + ME(1:2) * fi*(1-fj) +  ME(1:2)*fj*(1-fi)
$NTR
$NTR         if(inversetemp.eq.-1) then
$NTR           dfdE = fj - fi
$NTR           if(abs(dfdE).gt.0) then
$NTR              dfdE = dfdE/(spenergies(ii) - spenergies(jj))
$NTR           endif
$NTR         else
$NTR           if(abs(spenergies(ii) - spenergies(jj)).gt.1d-8) then
$NTR             dfdE = (fj - fi)/(spenergies(ii) - spenergies(jj))
$NTR           else                
$NTR             dfdE = fi**2 * inversetemp                                   &
$NTR             &          * exp(inversetemp*(spenergies(ii)-FermiEnergy(it)))
$NTR           endif
$NTR         endif
$NTR         Belyaev(1:2,it) = Belyaev(1:2,it) +  2*ME(1:2) * dfde
$NTR       endif
$NTR     enddo
$NTR   enddo
      si = si + N + N2
    enddo
    
#if(USE_MPI>0)
    call MPI_ALLREDUCE(MPI_IN_PLACE, Belyaev(:,1:2), 6, MPI_REAL8, MPI_SUM,    &
    &                  MPI_COMM_WORLD, MPI_ERR)
    call MPI_ALLREDUCE(MPI_IN_PLACE, J2_pairing_cut(:,1:2)     , 6, MPI_REAL8, MPI_SUM,    &
    &                  MPI_COMM_WORLD, MPI_ERR)
#endif

    ! Factor 2 for time-reversal
$TR J2_pairing_cut(:,1:2)      = 2 * J2_pairing_cut(:,1:2)
$TR Belyaev(:,1:2) = 2 * Belyaev(:,1:2)

    !  Sum for the total
    J2_pairing_cut(:,3)      = sum(J2_pairing_cut(:,1:2),2)
    Belyaev(:,3) = sum(Belyaev(:,1:2),2)
  end subroutine 

  subroutine calcJ2andBelyaev_BCS
    !---------------------------------------------------------------------------
    ! Calculate the expectation value of J^2_mu in the many-body state, as well
    ! as the Belyaev moment of inertia. 
    ! 
    ! Full expressions for both, in the BCS case at finite temperature 
    ! can be found in
    ! 
    !   Y. Alhassid et al., Phys. Rev. C 72, 064326 (2005).
    ! 
    !
    !  <J_i J_j> = sum_[ kl > 0 ]  ME_{kl} * weight_{kl}
    !    I_ij    = sum_[ kl > 0 ]  ME_{kl} * weight^B_{kl}
    !   
    ! where
    !      ME_{kl} = < k|j_i| l> < l|j_j| k>  + < l|j_i| k> < k|j_j| l>
    !              + < k|j_i|-l> <-l|j_j| k>  + <-k|j_i| l> < l|j_j|-k>
    !    
    ! and
    ! 
    !     weight_{kl}  = [u^2_k u_l^2 + u_k v_k u_l v_l] * f_k (1-f_l)     (a)
    !                  + [v^2_k v_l^2 + u_k v_k u_l v_l] * (1-f_k) f_l     (b)
    !                  + [u_k^2 v_l^2 - u_k v_k u_l v_l] * f_k f_l         (c)
    !                  + [v_k^2 u_l^2 - u_k v_k u_l v_l] * (1-f_k) (1-f_l) (d)
    !
    !
    !    weight^B_{kl} = (u_k u_l + v_k v_l)**2 * (f_l - f_k)/(E_k - E_l)  (Ba)
    !                  + (u_k v_l - v_k u_l)**2 * (1-f_k-f_l)/(E_k + E_l)  (Bb)
    !
    ! with the additional caveat that 
    !
    !   (f_l - f_k)/(E_k - E_l) => - df_k/dE_k as E_k => E_l
    !
    ! - df_k/dE_k = beta exp(beta * E_k)* f_k^2
    !
    ! Due to conserved signature and time-simplexes, the only matrix elements
    ! that enter the calculation of ME_{kl} are of the form
    !
    !  Re < i | J_x T  | j > 
    !  Im < i | J_y T  | j >
    !  Re < i | J_z    | j > 
    !
    ! In addition, we only calculate the diagonal moments of the Inertia tensor.
    !---------------------------------------------------------------------------
    ! Mental note to self.
    ! - - - - - - - - - - - -
    ! This routine has been checked
    !  *) Belyaev & J2 are close to values produced by EV8 for Zr84, 
    !     though not identical due to the difference in derivatives.
    !  *) BCS formulas and HF formulas agree
    !       a) when the pairing strength is set to zero, at zero temperature
    !       b) at high temperature after the pairing phase transition 
    !  *) Both HF and BCS formulas tend towards the zero temperature value as   
    !     beta gets increased.
    !  *) The formulas correctly produce 0 for spherical configurations at
    !     zero temperature.
    !---------------------------------------------------------------------------
    integer       :: i,j, b, it, ii, jj, si, N, locali, localj
    integer       :: calc_rank, ranki, rankj, designated_rank(8)
    real(KIND=dp) :: ME(3), uvi, uvj, ui, uj, vi, vj, fi, fj
    real(KIND=dp) :: wa, wb, wc, wd, Ba, Bb, dfde
    real(KIND=dp) :: psi_i(mv,4), psi_j(mv,4), der_psi_j(mv,3,4)
#if(USE_MPI>0)
    integer :: mpi_err
#endif

    J2_pairing_cut = 0  ;  Belyaev = 0

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    designated_rank = -1 
    do B=1,8
      ! First, we designate a rank to do the calculations
      !    = the first rank storing spwfs in a symmetry block
      N = HFBlocks_global(B) ;  if(N.eq.0) cycle
      designated_rank(B) = rank_map(sum(HFBlocks_global(1:B-1))+1)
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
    si = 0  
    do B = 1, Blocks   ! <----- loop over global indices of spwfs
      N = HFBlocks_global(B); if(N.eq.0) cycle
      calc_rank = designated_rank(B) ! <---- this rank will do the integrations
                                     !       of the spwf matrix elements

      do j=1, N
        jj    = si + j            ! global index of the spwf
        localj= spwf_inverse(jj)  ! local index of the spwf
        rankj = rank_map(jj)      ! MPI rank storing the spwf

#if(USE_MPI>0)
        call Transfer_psi(psi_j, localj, 'HF', rankj, calc_rank)
        call Transfer_derpsi_complete(der_psi_j,localj,'HF',rankj, calc_rank)
#else
        call Transfer_psi(psi_j, localj, 'HF')
        call Transfer_derpsi_complete(der_psi_j,localj,'HF')
#endif

        do i=1, N
          ii    = si + i            ! global index of the spwf
          locali= spwf_inverse(ii)  ! local index of the spwf
          ranki = rank_map(ii)      ! MPI rank storing the spwf
          it = 1

          if(B.gt.4) it = 2
#if(USE_MPI>0)
          call Transfer_psi(psi_i, locali, 'HF', ranki, calc_rank)
#else
          call Transfer_psi(psi_i, locali, 'HF')
#endif

          if(MPI_RANK.eq.calc_rank) then
            ! ^------ only one rank needs to to this calculation.

            ! BCS matrix elements ----------------------------------------------
            ! |< k | j_x | -l >|^2  
            ME(1)= angmom_xt_real(psi_i, psi_j, der_psi_j)**2
            ! |< k | j_y | -l >|^2 
            ME(2)= angmom_yt_imag(psi_i, psi_j, der_psi_j)**2
            ! |< k | j_z |  l >|^2 
            ME(3)= angmom_z_real (psi_i, psi_j, der_psi_j)**2
            !-------------------------------------------------------------------

            ! The factor two is because of the presence of the four terms in 
            ! ME_{kl}, which are pair-wise equal when J_i = J_j.
            ME = 2 * ME 

            vi  = BCSOccupations(ii)/2. ; vj  = BCSOccupations(jj)/2.
            ui  = 1 - vi                ; uj  = 1 - vj
            fi  = BCSf(ii)              ; fj  = BCSf(jj)
            uvi = ui*vi                 ; uvj = uj*vj         

            ! Take the square root, but take care for numerical errors producing 
            ! small negative numbers.
            if(uvi .gt. 0) then
              uvi = sqrt(uvi) 
            else
              uvi = 0
            endif
            if(uvj .gt. 0) then
               uvj = sqrt(uvj)
            else
               uvj = 0
            endif

            !-------------------------------------------------------------------
            ! <J^2>
            ! [u^2_k u_l^2 + u_k v_k u_l v_l] * f_k (1-f_l)                  (a)
            wa = (ui*uj + uvi*uvj) * fi * (1-fj)
            ! [v^2_k v_l^2 + u_k v_k u_l v_l] * (1-f_k) f_l                  (b)
            wb = (vi*vj + uvi*uvj) * (1-fi) * fj
            ! [u_k^2 v_l^2 - u_k v_k u_l v_l] * f_k f_l                      (c)
            wc = (ui*vj - uvi*uvj) * fi * fj
            ! [v_k^2 u_l^2 - u_k v_k u_l v_l] * (1-f_k) (1-f_l)              (d)
            wd = (vi*uj - uvi*uvj) * (1-fi) *(1-fj)

            J2_pairing_cut(:,it) = J2_pairing_cut(:,it) + ME * (wa+wb+wc+wd)
            !-------------------------------------------------------------------
            ! I_xx, I_yy and I_zz
            if(abs(BCSqps(ii) - BCSqps(jj)).gt.1d-5) then            
              dfde = (fj - fi)/(BCSqps(ii) - BCSqps(jj))
            else
              if(inversetemp.gt.0) then
                ! beta exp(beta * E_k)* f_k^2
                dfde = inversetemp * exp(inversetemp*BCSqps(ii)) * fi**2
              else
                dfde = 0 ! All QPS are unoccupied for T=0 calculation
              endif
            endif
            ! (u_k u_l + v_k v_l)**2 * (f_l - f_k)/(E_k - E_l)  (Ba)
            Ba = (ui*uj + vi*vj + 2*uvi*uvj) * dfdE 

            dfdE = (1 - fi - fj)/(BCSqps(ii) + BCSqps(jj))

            ! (u_k v_l - v_k u_l)**2 * (1-f_k-f_l)/(E_k + E_l)  (Bb)
            Bb = (ui*vj + uj*vi - 2*uvi*uvj) * dfde

            Belyaev(:,it) =  Belyaev(:,it) + ME * (Ba + Bb)
          endif
        enddo
      enddo
      si = si + N
    enddo
#if(USE_MPI>0)
    call MPI_ALLREDUCE(MPI_IN_PLACE, Belyaev(:,1:2), 6, MPI_REAL8, MPI_SUM,    &
    &                  MPI_COMM_WORLD, MPI_ERR)
    call MPI_ALLREDUCE(MPI_IN_PLACE, J2_pairing_cut(:,1:2)     , 6, MPI_REAL8, MPI_SUM,    &
    &                  MPI_COMM_WORLD, MPI_ERR)
#endif
    !  Sum for the total
    J2_pairing_cut(:,3) = sum(J2_pairing_cut(:,1:2),2)
    Belyaev(:,3)        = sum(Belyaev(:,1:2),2)
  end subroutine calcJ2andBelyaev_BCS

  subroutine calcJ2andBelyaev_HFB
    !---------------------------------------------------------------------------
    ! NOTES 
    !
    !  - The way this routine is parallelized is simple: parallel calculation
    !    of all relevant matrix elements of the angular momentum operators
    !    followed by a calculation of Belyaev by ALL MPI ranks.
    !---------------------------------------------------------------------------
    ! This routine calculates the 
    !
    !  (i) the dispersion of the diagonal elements of the inertia tensor 
    !
    !           < J^2_mm > - <J_m>^2
    ! 
    !   in the canonical basis
    !
    !              =  sum_{ab} |< a | j_m | b  >|^2 
    !                                             rho_aa (1 - rho_bb)   
    !              -  sum_{ab}  < a | j_m | b  > < \bar{b} | j_m | \bar{a}  > 
    !                                   kappa_{a\bar{a}} \kappa_{\bar{b} b}
    !            
    !            where \bar{a} is the canonically conjugate partner of a.
    !             (sum over all possible sps in the can. basis)
    ! 
    !   in the qp basis 
    !   
    !              =  1/2 sum_ab (1-f_a)(1-f_b) |J^{20}_ab|^2 
    !              +      sum_ab f_a (1-f_b) |J^{11}_ab|^2 
    !
    !     where the sums range over all possible combinations of qps. Note that
    !     for even-even nuclei, the f_a are zero and only the first sum 
    !     contributes. 
    !
    !      Some caveats apply:
    !
    !        (*) It is possible (and likely desireable) to include a cutoff
    !            in this calculation. This can (in my opinion) only be 
    !            meaningfully done in the HF basis, which is why the 
    !            the calculation starts with the calculation of the matrix
    !            elements of angular momentum in the HF basis. They are then
    !            transformed into the canonical basis. 
    !
    !        (*) For blocked calculations, we also include a collective 
    !            version of this dispersion, meaning we simply remove the 
    !            blocked quasi-particles from the summation. This is most 
    !            easily done in the canonical basis.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    !  and 
    !
    !  (ii) the Belyaev moment of inertia around the Cartesian axes, 
    !               
    !           I_mm = - partial_{\omega} < psi' | J_m | psi' > 
    !
    !        where 
    !        a) the derivative should be evaluated at a reference frequency 
    !            \omega 
    !        b) psi' is a second order perturbation to the HFB state when 
    !           the rotational frequency changes, i.e. 
    !       
    !               H - \omega J_m => H - \omega J_m  - do J_m
    !   
    !                   
    !           hence 
    !                                    
    !          | psi'> - |psi_0 > =  
    !             do               C
    !          -  --- sum_ab -------------  b^{\dagger}_a b^{\dagger}_a |psi_0 > 
    !              2           E_a + E_b
    !           
    !           when T != 0 or in the presence of quasiparticle excitations
    !           there are also terms involving annihilation operators
    !                
    !          with E_a and E_b the quasiparticle energies and 
    !
    !           C =  < psi_0 | J_m  b^{\dagger}_a b^{\dagger}_a |psi_0 > 
    !   
    ! This gives rise to 
    !
    !    I_mm = \sum_{ab} (1 - f_a  - f_b) (E_a + E_b)^{-1} |J^{20}|^2_{m,ab}
    !         + \sum_{ab} (f_b - f_a)      (E_a - E_b)^{-1} |J^{11}|^2_{m,ab}
    !
    ! where the sums are over all possible qp combinations.
    !
    ! --------------------------------------------------------------------------
    ! Notes to self about checks of this routine
    ! 
    !  (a) Collective and ordinary quantities are equal for non-blocked runs
    !  (b) Time-reversal conserving results equal the BCS ones
    !  (c) When pairing collapses, results are equal to the HF ones.
    !       (both with and without T conservation)
    !  (d) When including all quasiparticles in the sums, 
    !       < \Delta J^2> calculated in the HF basis and in the qp-basis
    !      are equal, also for blocked calculations.
    !
    !---------------------------------------------------------------------------
    integer       :: i,j, b, it, ii, iii, jjj, jj, si, N,k, sb, ibar, jbar, N2,T
    integer       :: locali, localj, ranki, rankj, calc_rank, designated_rank(8)
    integer       :: wave, der_index
    real(KIND=dp) :: ME(3),  fac, psi_i(mv,4), psi_j(mv,4), der_psi_j(mv,3,4)

    !----------------------------------------------------------------------------
    ! Single-particle matrix elements of Jx, Jy, Jz in the array labelled HFBasis
    real(KIND=dp) :: jx(nwt,nwt), jy(nwt,nwt), jz(nwt,nwt)
    ! Single-particle matrix elements of Jx, Jy, Jz in the canonical basis
    real(KIND=dp) :: jx_can(nwt,nwt), jy_can(nwt,nwt), jz_can(nwt,nwt)
    !  and WITH the pairing cutoff folded in
    real(KIND=dp) :: jx_can_cut(nwt,nwt), jy_can_cut(nwt,nwt), jz_can_cut(nwt,nwt)
    ! Same, but Jx^2, Jy^2 and Jz^2 and only the diagonal elements
    real(KIND=dp) :: jx2_can(nwt), jy2_can(nwt), jz2_can(nwt)
    !---------------------------------------------------------------------------
    real(KIND=dp) :: J20(nwt,nwt, 3), J11(nwt,nwt,3) , cut_cr
    logical       ::  blocked

#if(USE_MPI>0)
    integer       :: mpi_err
#endif

    Belyaev = 0         ; Bely_coll      = 0
    J2 = 0              ; J2_pairing_cut = 0  ; J2_coll = 0
    jx = 0              ; jy = 0              ; jz = 0
    jx2_can = 0.0d0     ; jy2_can = 0.0d0     ; jz2_can = 0.0d0

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    designated_rank = -1 
    do B=1,8
      ! First, we designate a rank to do the calculations
      !    = the first rank storing spwfs in a symmetry block
      N = HFBlocks_global(B) ;  if(N.eq.0) cycle
      designated_rank(B) = rank_map(sum(HFBlocks_global(1:B-1))+1)
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! ------------------- START of all things parallel -------------------------

    !---------------------------------------------------------------------------
    ! First, we calculate the full matrix elements of jx, jy and jz and their
    ! squares in the Hartree-Fock basis.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    ! If time-reversal is not conserved the matrices jx/jy/jz are 
    !
    !   jx_ij = Re < i | j_x | j >
    !   jy_ij = Im < i | j_y | j >
    !   jz_ij = Re < i | j_z | j >
    !
    ! jz_ij is a always a real number, and only non-zero if i and j are in the 
    ! same 'half' of the basis (equal signature in CR8-like symmetries).
    !   
    ! jx_ij is always real, jy_ij is always imaginary and both are only 
    ! non-zero if i and j belong to opposite "halves" of the basis 
    ! (different signatures in CR8-like symmetry options).
    !
    ! All matrix elements if the parity of i does not match the parity of j.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! If time-reversal is conserved the matrices jx/jy/jz are
    !       
    !  jx_ij = Re < i | j_x | \bar{j} >
    !  jy_ij = Im < i | j_y | \bar{j} >
    !  jz_ij = Re < i | j_z |      j  >
    !
    ! where \bar{j} is the partner under time-reversal of j. 
    !
    ! jz_ij and jx_ij are always real, and jy_ij is always imaginary.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    si = 0
    do B = 1,8,2
      N = HFBlocks_global(B)   ; if(N.eq.0) cycle   !<--- this is a loop over
      N2= HFBlocks_global(B+1)                      !     global indices
      T = N + N2
      it = 1 ; if(B.gt.4) it=2
      calc_rank = designated_rank(B) ! this is the rank doing the calculation
      !-------------------------------------------------------------------------
      ! If time-reversal is conserved, we can treat all directions equally:
      ! both the inner and outer summation are over the same symmetry block    
      ! for all directions.
      !-------------------------------------------------------------------------
      do j=1,N
        jj    = si + j            ! global index of the spwf
        localj= spwf_inverse(jj)  ! local index of the spwf
        rankj = rank_map(jj)      ! MPI rank storing the spwf

#if(USE_MPI>0)
        call Transfer_psi(psi_j, localj, 'HF', rankj, calc_rank)
        call Transfer_derpsi_complete(der_psi_j,localj,'HF',rankj, calc_rank)
#else
        call Transfer_psi(psi_j, localj, 'HF')
        call Transfer_derpsi_complete(der_psi_j,localj,'HF')
#endif
        do i=j,N
          ii    = si + i            ! global index of the spwf
          locali= spwf_inverse(ii)  ! local index of the spwf
          ranki = rank_map(ii)      ! MPI rank storing the spwf
#if(USE_MPI>0)
          call Transfer_psi(psi_i, locali, 'HF', ranki, calc_rank)
#else
          call Transfer_psi(psi_i, locali, 'HF')
#endif

          if(MPI_rank.eq.calc_rank) then
$TR         ! < k | j_x | -l >
$TR         jx(ii,jj)= angmom_xt_real(psi_i,psi_j,der_psi_j) 
$TR         jx(jj,ii)= jx(ii,jj)
$TR         ! < k | j_y | -l >
$TR         jy(ii,jj)= angmom_yt_imag(psi_i,psi_j,der_psi_j) 
$TR         jy(jj,ii)= jy(ii,jj)

            ! < k | j_z |  l >
            jz(ii  ,jj  ) = angmom_z_real(psi_i,psi_j,der_psi_j) 
            jz(jj  ,ii  ) = jz(ii,jj)
          endif
        enddo
      enddo
      ! If time-reversal is not conserved, we have only calculated half of the 
      ! necessary matrix elements of jz above
$NTR  calc_rank = designated_rank(B+1) ! this is the rank doing the calculation
$NTR  do j=1,N2
$NTR   jj = si + N + j           ! global index of the spwf
$NTR   localj= spwf_inverse(jj)  ! local index of the spwf
$NTR   rankj = rank_map(jj)      ! MPI rank storing the spwf
#if(USE_MPI>0)
$NTR   call Transfer_psi(psi_j, localj, 'HF', rankj, calc_rank)
$NTR   call Transfer_derpsi_complete(der_psi_j,localj,'HF',rankj, calc_rank)
#else
$NTR   call Transfer_psi(psi_j, localj, 'HF')
$NTR   call Transfer_derpsi_complete(der_psi_j,localj,'HF')
#endif
$NTR   do i=j,N2
$NTR      ii = si + N + i
$NTR      locali= spwf_inverse(ii)  ! local index of the spwf
$NTR      ranki = rank_map(ii)      ! MPI rank storing the spwf
#if(USE_MPI>0)
$NTR      call Transfer_psi(psi_i, locali, 'HF', ranki, calc_rank)
#else
$NTR      call Transfer_psi(psi_i, locali, 'HF')
#endif
$NTR      if(calc_rank.eq.MPI_RANK) then
$NTR        jz(ii,jj)= angmom_z_real(psi_i,psi_j,der_psi_j)
$NTR        jz(jj,ii)= jz(ii,jj)
$NTR      endif
$NTR     enddo
$NTR  enddo
      !-------------------------------------------------------------------------
      ! If time-reversal is not conserved, the inner and outer loops for the 
      ! jx and jy matrix elements are not the same.
      !-------------------------------------------------------------------------
$NTR  calc_rank = designated_rank(B) ! this is the rank doing the calculation
$NTR  do j=1, N2
$NTR    jj = si + N +  j          ! global index of the spwf
$NTR    localj= spwf_inverse(jj)  ! local index of the spwf
$NTR    rankj = rank_map(jj)      ! MPI rank storing the spwf
#if(USE_MPI>0)
$NTR    call Transfer_psi(psi_j, localj, 'HF', rankj, calc_rank)
$NTR    call Transfer_derpsi_complete(der_psi_j,localj,'HF',rankj, calc_rank)
#else
$NTR    call Transfer_psi(psi_j, localj, 'HF')
$NTR    call Transfer_derpsi_complete(der_psi_j,localj,'HF')
#endif
$NTR    do i=1,N
$NTR      ii = si + i               ! global index of the spwf
$NTR      locali= spwf_inverse(ii)  ! local index of the spwf
$NTR      ranki = rank_map(ii)      ! MPI rank storing the spwf
#if(USE_MPI>0)
$NTR      call Transfer_psi(psi_i, locali, 'HF', ranki, calc_rank)
#else
$NTR      call Transfer_psi(psi_i, locali, 'HF')
#endif
$NTR      if(calc_rank.eq.MPI_RANK) then
$NTR        !|< k | j_x | l >|^2
$NTR        jx(ii,jj)=angmom_x_real(psi_i, psi_j, der_psi_j)
$NTR        ! |< k | j_y | l >|^2
$NTR        jy(ii,jj)=angmom_y_imag(psi_i, psi_j, der_psi_j)
$NTR        jx(jj,ii) =  jx(ii,jj)
$NTR        jy(jj,ii) = -jy(ii,jj)  ! attention to the sign(s)!
$NTR      endif
$NTR    enddo
$NTR  enddo

      si = si + N + N2
    enddo

    !----------------------------------------------------------------------
    ! The calculation of the diagonal matrix elements <J^2> is much easier
    ! since every rank can fend for its own.
    do i=1,nwt_local     ! loop over local indices
      wave = spwf_map(i) ! global index

      if(store_derivatives) then
        ! In this case we have the derivatives in storage
        der_index = i
      else
        ! In this case, we recalculate the derivatives on the fly
        call Derive_X_spwf(canpsi(:,:,i), sx(:,i), candpsi(:,1,:,1))
        call Derive_Y_spwf(canpsi(:,:,i), sy(:,i), candpsi(:,2,:,1))
        call Derive_Z_spwf(canpsi(:,:,i), sz(:,i), candpsi(:,3,:,1))
        der_index = 1
      endif

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Angular momenta
      jx2_can(wave)     = angmom_x_quad(canpsi(:,:,i), candpsi(:,:,:,der_index), &
      &                                 canpsi(:,:,i), candpsi(:,:,:,der_index))
      jy2_can(wave)     = angmom_y_quad(canpsi(:,:,i), candpsi(:,:,:,der_index), &
      &                                 canpsi(:,:,i), candpsi(:,:,:,der_index))
      jz2_can(wave)     = angmom_z_quad(canpsi(:,:,i), candpsi(:,:,:,der_index), &
      &                                 canpsi(:,:,i), candpsi(:,:,:,der_index))
    enddo

#if(USE_MPI>0)
    ! These MPI ALLREDUCE calls can be improved upon:
    !   by using BCAST from each calc_rank and restricting to symmetry blocks
    call MPI_ALLREDUCE(MPI_IN_PLACE, jx, nwt**2, MPI_REAL8, MPI_SUM,    &
    &                  MPI_COMM_WORLD, MPI_ERR)
    call MPI_ALLREDUCE(MPI_IN_PLACE, jy, nwt**2, MPI_REAL8, MPI_SUM,    &
    &                  MPI_COMM_WORLD, MPI_ERR)
    call MPI_ALLREDUCE(MPI_IN_PLACE, jz, nwt**2, MPI_REAL8, MPI_SUM,    &
    &                  MPI_COMM_WORLD, MPI_ERR)

    call MPI_ALLREDUCE(MPI_IN_PLACE, jx2_can, nwt, MPI_REAL8, MPI_SUM,    &
    &                  MPI_COMM_WORLD, MPI_ERR)
    call MPI_ALLREDUCE(MPI_IN_PLACE, jy2_can, nwt, MPI_REAL8, MPI_SUM,    &
    &                  MPI_COMM_WORLD, MPI_ERR)
    call MPI_ALLREDUCE(MPI_IN_PLACE, jz2_can, nwt, MPI_REAL8, MPI_SUM,    &
    &                  MPI_COMM_WORLD, MPI_ERR)
#endif

    ! ------------------- END of all things parallel ---------------------------

    si = 0
    do B = 1,8,2
      N = HFBlocks_global(B)   ; if(N.eq.0) cycle   !<--- this is a loop over
      N2= HFBlocks_global(B+1)                      !     global indices
      T = N + N2
      it = 1 ; if(B.gt.4) it=2
      !-------------------------------------------------------------------------
      ! First calculate matrix elements of Jx,Jy,Jz in the canonical basis
      jx_can(si+1:si+T, si+1:si+T) = &
      &  matmul(transpose(cantransfo(si+1:si+T, si+1:si+T)), &
      &                              jx(si+1:si+T, si+1:si+T))
      jy_can(si+1:si+T, si+1:si+T) = &
      &  matmul(transpose(cantransfo(si+1:si+T, si+1:si+T)), &
      &                              jy(si+1:si+T, si+1:si+T))
      jz_can(si+1:si+T, si+1:si+T) = &
      & matmul(transpose(cantransfo(si+1:si+T, si+1:si+T)), &
      &                              jz(si+1:si+T, si+1:si+T))

      jx_can(si+1:si+T, si+1:si+T) = &
      & matmul(jx_can(si+1:si+T, si+1:si+T), &
      &        cantransfo(si+1:si+T, si+1:si+T))

      jy_can(si+1:si+T, si+1:si+T) = &
      & matmul(jy_can(si+1:si+T, si+1:si+T), &
      &        cantransfo(si+1:si+T, si+1:si+T))

      jz_can(si+1:si+T, si+1:si+T) = &
      & matmul(jz_can(si+1:si+T, si+1:si+T), &
      &        cantransfo(si+1:si+T, si+1:si+T))

      !-------------------------------------------------------------------------
      ! ... and then proceed to calculate the matrix elements with cutoff
      if(rotcorr_cut) then
        !-----------------------------------------------------------------------
        ! The matrix elements computed above are the matrix elements in the
        ! sp-basis in storage. This is not necessarily the HF-basis.
        if(.not. diagsphamil) then
             jx(si+1:si+T,si+1:si+T) = &
             &    matmul(        jx(si+1:si+T,si+1:si+T),          &
             &                   HFtransfo(si+1:si+T,si+1:si+T))
             jx(si+1:si+T,si+1:si+T) = &
             &  matmul(transpose(HFtransfo(si+1:si+T,si+1:si+T)), &
             &                   jx(si+1:si+T,si+1:si+T))

             jy(si+1:si+T,si+1:si+T) = &
             &    matmul(        jy(si+1:si+T,si+1:si+T),          &
             &                   HFtransfo(si+1:si+T,si+1:si+T))
             jy(si+1:si+T,si+1:si+T) = &
             &  matmul(transpose(HFtransfo(si+1:si+T,si+1:si+T)), &
             &                   jy(si+1:si+T,si+1:si+T))

             jz(si+1:si+T,si+1:si+T) = &
             &    matmul(        jz(si+1:si+T,si+1:si+T),          &
             &                   HFtransfo(si+1:si+T,si+1:si+T))
             jz(si+1:si+T,si+1:si+T) = &
             &  matmul(transpose(HFtransfo(si+1:si+T,si+1:si+T)), &
             &                   jz(si+1:si+T,si+1:si+T))
        endif

        ! Apply the cutoff in the HF basis
        do i=1,T
          do j=1,T
            cut_cr = rotcut(spenergies(si+i), it) * rotcut(spenergies(si+j), it)
            jx(si+i,si+j) = jx(si+i,si+j) * cut_cr
            jy(si+i,si+j) = jy(si+i,si+j) * cut_cr
            jz(si+i,si+j) = jz(si+i,si+j) * cut_cr
          enddo
        enddo
        ! Transform back (if necessary)
        if(.not. diagsphamil) then
             jx(si+1:si+T,si+1:si+T) = &
             &    matmul(          jx(si+1:si+T,si+1:si+T),          &
             &           transpose(HFtransfo(si+1:si+T,si+1:si+T)))
             jx(si+1:si+T,si+1:si+T) = &
             &    matmul(          HFtransfo(si+1:si+T,si+1:si+T), &
             &                     jx(si+1:si+T,si+1:si+T))

             jy(si+1:si+T,si+1:si+T) = &
             &    matmul(          jy(si+1:si+T,si+1:si+T),          &
             &           transpose(HFtransfo(si+1:si+T,si+1:si+T)))
             jy(si+1:si+T,si+1:si+T) = &
             &    matmul(          HFtransfo(si+1:si+T,si+1:si+T), &
             &                     jy(si+1:si+T,si+1:si+T))

             jz(si+1:si+T,si+1:si+T) = &
             &    matmul(          jz(si+1:si+T,si+1:si+T),          &
             &           transpose(HFtransfo(si+1:si+T,si+1:si+T)))
             jz(si+1:si+T,si+1:si+T) = &
             &  matmul(            HFtransfo(si+1:si+T,si+1:si+T), &
             &                     jz(si+1:si+T,si+1:si+T))
        endif
      endif
      !-------------------------------------------------------------------------
      ! Transform the sp. matrix elements into the canonical basis.
      ! Note that we could have directly calculated these matrix elements in
      ! the canonical basis, but for one thing: the pairing cutoff which cannot
      ! be meaningully defined in the canonical basis. 
      !-------------------------------------------------------------------------
      jx_can_cut(si+1:si+T, si+1:si+T) = &
      &  matmul(transpose(cantransfo(si+1:si+T, si+1:si+T)), &
      &                              jx(si+1:si+T, si+1:si+T))
      jy_can_cut(si+1:si+T, si+1:si+T) = &
      &  matmul(transpose(cantransfo(si+1:si+T, si+1:si+T)), &
      &                              jy(si+1:si+T, si+1:si+T))

      jz_can_cut(si+1:si+T, si+1:si+T) = &
      & matmul(transpose(cantransfo(si+1:si+T, si+1:si+T)), &
      &                              jz(si+1:si+T, si+1:si+T))

      jx_can_cut(si+1:si+T, si+1:si+T) = &
      & matmul(jx_can_cut(si+1:si+T, si+1:si+T), &
      &        cantransfo(si+1:si+T, si+1:si+T))

      jy_can_cut(si+1:si+T, si+1:si+T) = &
      & matmul(jy_can_cut(si+1:si+T, si+1:si+T), &
      &        cantransfo(si+1:si+T, si+1:si+T))

      jz_can_cut(si+1:si+T, si+1:si+T) = &
      & matmul(jz_can_cut(si+1:si+T, si+1:si+T), &
      &        cantransfo(si+1:si+T, si+1:si+T))

      si = si + N + N2
    enddo
    !---------------------------------------------------------------------------
    ! Then we calculate the dispersion of J^2, in the canonical basis.
    !           < \Delta J^2 >  = < J^2 > - < J >^2
    ! 
    !---------------------------------------------------------------------------
    si = 0 
    do b = 1, Blocks,2
      N  = HFBlocks_global(b)   ; if(N.eq.0) cycle
      N2 = HFBlocks_global(b+1)

      it = 1 ; if(B.gt.4) it = 2

      do i=1, N+N2
        ii   = si + i
$TR     ibar = ii
$NTR    ibar = conjugp(ii)

        !--------------------------------------------------------------------
        ! Calculation of J2: terms with one index
        !
        !
        !  J2 ~  \sum_{ab} < a | J_{\mu}^2 | a > \rho_aa
        !--------------------------------------------------------------------
        fac=  rho_can(ii)      ! There is an additional factor 2  hidden in here
                               ! when TR is conserved, but it's easier to let
                               ! it be here and not have to cancel it in the
                               ! J2 summation below.
        ME(1) = jx2_can(ii)
        ME(2) = jy2_can(ii)
        ME(3) = jz2_can(ii)

        J2(:,it) = J2(:,it) + ME * fac

        do j=1,N+N2
          jj = si + j
$TR       jbar = jj
$NTR      jbar = conjugp(jj)
          !--------------------------------------------------------------------
          ! Calculation of J2: terms with two indices
          !
          !   J2 ~ - \sum_{ab} < a | J_{\mu} | b > < b | J_{\mu} | a > [ \rho_aa \rho_bb]
          !
          ! If time-reversal is conserved,
          !  we can treat the rho-rho term and kappa-kappa term equally
          !
          ! Factors 1./2 due to time-reversal
          !--------------------------------------------------------------------
$TR       fac=  - rho_can(ii)*rho_can(jj)/4.0d0 - kappa_can(ii)*kappa_can(jj)
$TR
$TR       ME(1) = 2*jx_can(ii,jj)**2 ! Factor two for time-reversal
$TR       ME(2) = 2*jy_can(ii,jj)**2 ! Factor two for time-reversal
$TR       ME(3) = 2*jz_can(ii,jj)**2 ! Factor two for time-reversal
$TR
$TR       J2(:,it) = J2(:,it) + ME * fac

          ! If time-reversal is broken; we can not do things quite that easily.
          !  => the rho-rho term and kappa-kappa term use matrix elements of
          !     of different states
$NTR      fac= -rho_can(ii) * rho_can(jj)
$NTR
$NTR      ME(1) = jx_can(ii,jj)**2
$NTR      ME(2) = jy_can(ii,jj)**2
$NTR      ME(3) = jz_can(ii,jj)**2
$NTR
$NTR      J2(:,it) = J2(:,it) + ME * fac
          ! We haven't necessarily found canonical partners for all states
          ! If a partner is absent, this means that the relevant matrix
          ! elements of kappa are too small anyway; we can safely forget about
          ! this term.
$NTR      if(ibar .ne. 0 .and. jbar .ne. 0) then
$NTR        fac= -kappa_can(ii)*kappa_can(jbar)
$NTR        ME(1) = jx_can(ii,jj)*jx_can(jbar,ibar)
$NTR        ME(2) = jy_can(ii,jj)*jy_can(jbar,ibar)
$NTR        ME(3) = jz_can(ii,jj)*jz_can(jbar,ibar)
$NTR        J2(:,it) = J2(:,it) + ME(:) * fac
$NTR      endif
          !--------------------------------------------------------------------
          ! b) calculation of J2_pairing_cut, i.e. with pairing cutoffs
          !--------------------------------------------------------------------
          ! If time-reversal is conserved,
          !  we can treat the rho-rho term and kappa-kappa term equally
          ! 
$TR       ! Factors 1./2 due to time-reversal
$TR       fac=  rho_can(ii)/2.*(1-rho_can(jj)/2.)-kappa_can(ii)*kappa_can(jj)
$TR       ME(1) = 2*jx_can_cut(ii,jj)**2 ! Factor two for time-reversal
$TR       ME(2) = 2*jy_can_cut(ii,jj)**2 ! Factor two for time-reversal
$TR       ME(3) = 2*jz_can_cut(ii,jj)**2 ! Factor two for time-reversal
$TR
$TR       J2_pairing_cut(:,it) = J2_pairing_cut(:,it) + ME * fac
          ! If time-reversal is broken; we can not do things quite that easily.
          !  => the rho-rho term and kappa-kappa term use matrix elements of
          !     of different states
$NTR      fac=  rho_can(ii)   *(1-rho_can(jj))
$NTR      ME(1) = jx_can_cut(ii,jj)**2
$NTR      ME(2) = jy_can_cut(ii,jj)**2
$NTR      ME(3) = jz_can_cut(ii,jj)**2
$NTR      J2_pairing_cut(:,it) = J2_pairing_cut(:,it) + ME * fac
          ! We haven't necessarily found canonical partners for all states
          ! If a partner is absent, this means that the relevant matrix
          ! elements of kappa are too small anyway; we can safely forget about
          ! this term.
$NTR      if(ibar .ne. 0 .and. jbar .ne. 0) then
$NTR        fac= -kappa_can(ii)*kappa_can(jbar)
$NTR        ME(1) = jx_can_cut(ii,jj)*jx_can_cut(jbar,ibar)
$NTR        ME(2) = jy_can_cut(ii,jj)*jy_can_cut(jbar,ibar)
$NTR        ME(3) = jz_can_cut(ii,jj)*jz_can_cut(jbar,ibar)
$NTR        J2_pairing_cut(:,it) = J2_pairing_cut(:,it) + ME(:) * fac
$NTR      endif
        enddo
      enddo
      si = si +  N + N2
    enddo
    J2(:,3)              = sum(J2(:,1:2),2)
    J2_pairing_cut(:,3)  = sum(J2_pairing_cut(:,1:2),2)
    !---------------------------------------------------------------------------
    ! I also calculate some approximation for the collective angular momentum, 
    ! which I define as <J^2> without the contribution from the blocked qps. 
    ! 
    !  To safely remove contributions of individual qps, this calculation is 
    !  done in the qp basis. 
    !  
    !  <J^2> = 1/2 sum_{ab} (1-f_a)(1-f_b) |J^20_ab|^2
    !        +     sum_{ab}  f_a (1-f_b)   |J^11_ab|^2     
    ! 
    !  where the sum simply does not include the blocked qps, but in general
    ! ranges over all possible other combinations.
    !   
    ! This in addition serves as an additional sanity check: without blocking 
    ! this collective value should equal the ordinary one calculated above.
    !---------------------------------------------------------------------------

    ! First, construct J20
    call calcJ20(jx,bogoliubov, j20(:,:,1))
    call calcJ20(jy,bogoliubov, j20(:,:,2)) 
    call calcJ20(jz,bogoliubov, j20(:,:,3)) 

    !Then , construct J11
    call calcJ11(jx,bogoliubov, j11(:,:,1)) 
    call calcJ11(jy,bogoliubov, j11(:,:,2)) 
    call calcJ11(jz,bogoliubov, j11(:,:,3)) 

    if(inversetemp .lt. 0) then
      si = 0 ; sb = 0
      do b=1,Blocks,2
        N  = hfblocks_global(b)    ; if(N .eq. 0) cycle
        N2 = hfblocks_global(b+1)
        
        it = 1 ; if(B.gt.4) it = 2
        do i=1, N + N2
          ii = si + i
          iii= sb + N + N2 + i

          do j=1,N+N2
            jj = si + j
            jjj= sb + N + N2 + j

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Don't include the contribution from the blocked qps and their
            ! partner qps.
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            blocked= .false.
            if(allocated(blocked_qps)) then
              do k=1,size(blocked_qps)
                if(ii.eq.blocked_qps(k) .or.  jj.eq.blocked_qps(k)) then
                  blocked = .true.
                endif
                if(ii.eq.partner_qps(k) .or.  jj.eq.partner_qps(k)) then
                  blocked = .true.
                endif
              enddo              
            endif
            if(blocked) cycle
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !J^20 contribution
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            fac =  configmatrix(iii) * configmatrix(jjj)
            ME = 0.5 * J20(ii,jj,:)**2  * fac
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! J^11 contribution
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ! Note: this should NEVER contribute to the collective J^2. 
            !       I've coded it only as a sanity check, to check whether
            !       the full qp-basis expression is the same as the 
            !       sp-basis summation when the blocked QPS are not ommitted.
            fac=  configmatrix(jjj)*(1 - configmatrix(iii))
            ME = ME + J11(ii,jj,:)**2  * fac
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            J2_coll(:,it) = J2_coll(:,it) + ME      
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          enddo
        enddo
        si = si +   N +   N2
        sb = sb + 2*N + 2*N2
      enddo
$TR   J2_coll(:,1:2) = 2 * J2_coll(:,1:2)   ! Time-reversal factor two 
      J2_coll(:,3) = sum(J2_coll(:,1:2), 2)
    endif
    !---------------------------------------------------------------------------
    ! Then the Belyaev moment of inertia in the ordinary sp. basis.
    ! 
    !  From expanding the many-body wave-function around the HFB minimum for 
    !  small rotational frequency omega, we get the following expression
    ! 
    !     I_mm = \sum_{ab} (1 - f_a  - f_b) (E_a + E_b)^{-1} |J^{20}|^2_{m,ab}
    !          + \sum_{ab} (f_b - f_a)      (E_a - E_b)^{-1} |J^{11}|^2_{m,ab}
    ! 
    ! based on (the generalisation of) equation 3.92 on pg 131 in Ring and 
    ! Schuck, . In the case of time-reversal conservation (and no blocking), 
    ! the sum over (ab) gets restricted and the f_a vanish. We have 
    ! 
    !   I_{mm} = 2 \sum_{ab>0}  (|J^{20}|^2_{m,ab} + |J^{20}|^2_{m,a\bar{b}})
    !                           ---------------------------------------------
    !                                               E_a + E_b
    !---------------------------------------------------------------------------

    si = 0 ; sb = 0
    do b = 1, Blocks, 2
      N = HFBlocks_global(b)   ; if(N.eq.0) cycle
      N2= HFBlocks_global(b+1)
      it = 1 ; if(B.gt.4) it = 2

      do i=1, N + N2
        ii = si + i
        iii= sb + N + N2 + i

        do j=1, N + N2
          jj = si + j
          jjj= sb + N + N2 + j
          !---------------------------------------------------------------------
          ! Ordinary  (full) Belyaev calculation
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !        (1 - f_i - f_j) J^{20}^2 contribution
          fac = 1 - configmatrix(sb+i) - configmatrix(sb+j)
          Belyaev(:,it) = Belyaev(:,it) + &
          &             fac*J20(ii,jj,:)**2 /(Qpenergies(iii) + Qpenergies(jjj))    
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !     (f_j - f_i) J^{11} contribution (with numerical safety)
          if(abs(Qpenergies(iii) - Qpenergies(jjj)) .gt. 1d-8) then
            fac =  configmatrix(sb+j) - configmatrix(sb+i)
            Belyaev(:,it) = Belyaev(:,it) + &
            &       fac*J11(ii,jj,:)**2 /(Qpenergies(iii)-Qpenergies(jjj))  
          elseif(inversetemp .gt. 0) then
           ! 28/12/2020, WR: I'm unsure whether there should be a factor 2
           ! here or not.... To be doublechecked.
           ! degen = inversetemp * configmatrix(sb+i)**2 *                     &
           ! &                                 exp(inversetemp * Qpenergies(ii))
           ! Belyaev(:,it) = Belyaev(:,it) +   J11(ii,jj,:)**2 * degen   	
           call stp('The calculation of moments of the inertia at finite T is &
           &         not supported yet.')
          endif 
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! Collective Belyaev calculation  
          if(inversetemp.lt.0) then
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Don't include the blocked_qps and their partner qps
            blocked= .false.

            if(allocated(blocked_qps)) then
              do k=1,size(blocked_qps)
                if(ii.eq.blocked_qps(k) .or.  jj.eq.blocked_qps(k)) then
                  blocked = .true.
                endif
                if(ii.eq.partner_qps(k) .or.  jj.eq.partner_qps(k)) then
                  blocked = .true.
                endif
              enddo
            endif
            if(blocked) cycle

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! (1 - f_i - f_j) J^{20}^2 contribution
            fac = 1 - configmatrix(sb+i) - configmatrix(sb+j) 
            Bely_coll(:,it) = Bely_coll(:,it) + &
            &           fac*J20(ii,jj,:)**2 /(Qpenergies(iii) + Qpenergies(jjj))  
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  (f_j - f_i) J^{11} contribution (with numerical safety)
            fac =  configmatrix(sb+j) - configmatrix(sb+i)
            if(abs(Qpenergies(iii) - Qpenergies(jjj)) .gt. 1d-8) then
              Bely_coll(:,it) = Bely_coll(:,it) + &
              &       fac*J11(ii,jj,:)**2 /(Qpenergies(iii)-Qpenergies(jjj))  
            endif  
          endif        
          !---------------------------------------------------------------------
        enddo
      enddo
      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
    enddo
    !---------------------------------------------------------------------------
    ! Final considerations    
$TR Belyaev(:,1:2)   = 2*Belyaev(:,1:2)   ! Time-reversal factor 2's
$TR Bely_coll(:,1:2) = 2*Bely_coll(:,1:2) ! Time-reversal factor 2's

    Belyaev(:,3)   = sum(Belyaev(:,1:2),2)
    Bely_coll(:,3) = sum(Bely_coll(:,1:2), 2)
  end subroutine calcJ2andBelyaev_HFB  

  subroutine calcJ20(j,bogo, j20)
    !---------------------------------------------------------------------------
    ! Subroutine to calculate the matrix J20 in the Ring and Schuck 
    ! notation.
    !
    !  J20 = U^\dagger j V^* - V^\dagger j^t U^*
    !
    ! The matrix j on input contains the single-particle matrix elements of 
    ! the angular momentum operator, bogo is the bogoliubov transformation.
    !
    !---------------------------------------------------------------------------
    ! When time-reversal is conserved, one should exercice caution. 
    ! The full U and V matrices are 
    !
    !   ( U^+ 0  )    and  ( 0    V^-)
    !   ( 0   U^-)         ( V^+  0  )
    !
    ! and only the U^+ and V^+ are actually stored by the code, and 
    !
    !      V^- = - V^+
    !      U^- = + U^+
    !
    ! So, if the single-particle matrix elements passed in are just 
    !
    !                < a | j_mu | b > 
    !
    ! this routine actually returns - J^{20}_{a\bar{b}} where \bar{b} is the 
    ! time-reversed partner of b.  
    ! 
    ! If instead, the single-particle matrix elements that are passed in are
    ! 
    !                < a | j_mu T | b >, 
    ! 
    ! the result of this routine is indeed J^{20}_{ab}.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! When time-reversal is not conserved, the full U and V matrices are 
    ! stored, but not pair-wise equal anymore. In that case, the routine 
    ! behaves as one would expect naively.
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in)  :: j(nwt,nwt)
    real(KIND=dp), intent(in)  :: bogo(2*nwt, 2*nwt)
    real(KIND=dp), intent(out) :: j20(nwt,nwt)
    real(KIND=dp)              :: tmp(nwt,nwt)
    integer :: N, N2, si, sb, b, T

    si = 0 ; sb = 0
    do B = 1, 8, 2
      N = HFBlocks_global(B)   ; if(N.eq.0) cycle
      N2= HFblocks_global(B+1)

      T =  N + N2
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! U^\dagger j 
      j20(si+1:si+T, si+1:si+T) = &
      & matmul(transpose(bogo(sb+1:sb+T, sb+T+1:sb+2*T)),j(si+1:si+T,si+1:si+T))

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! U^\dagger j V^*
      j20(si+1:si+T, si+1:si+T) = &
      & matmul( j20(si+1:si+T, si+1:si+T), bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T))

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! V^\dagger j^t 
      tmp(si+1:si+T, si+1:si+T) = transpose(j(si+1:si+T, si+1:si+T))
      tmp(si+1:si+T, si+1:si+T) = &
      & matmul(transpose(bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)),                  &
      &                             tmp(si+1:si+T,si+1:si+T))
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! V^\dagger j^T U
      j20(si+1:si+T, si+1:si+T) = j20(si+1:si+T, si+1:si+T) - &
      &  matmul( tmp(si+1:si+T, si+1:si+T), bogo(sb+1:sb+T, sb+T+1:sb+2*T))

      si = si +   T 
      sb = sb + 2*T 
    enddo
  end subroutine calcJ20

  subroutine calcJ11(j,bogo, j11)
    !---------------------------------------------------------------------------
    ! Small subroutine to calculate the matrix J11 in the Ring and Schuck 
    ! notation.
    !
    !  J11 = U^\dagger j U - V^\dagger j^t V^*
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    ! When time-reversal is conserved, one should exercice caution. 
    ! The full U and V matrices are 
    !
    !   ( U^+ 0  )    and  ( 0    V^-)
    !   ( 0   U^-)         ( V^+  0  )
    !
    ! and only the U^+ and V^+ are actually stored by the code, and 
    !
    !      V^- = - V^+
    !      U^- = + U^+
    !
    ! So the full calculation of the second term in that case would be
    !
    ! ( 0   V^+) ( j+      j+-  ) ( 0   V^- )
    ! ( V^- 0  ) ( j+-^*  -j-   ) ( V^+ 0   )
    !
    !  =  ( - V^+ j+  V^+  -V^+ j+-^* V^+ )
    !     ( - V^+ j+- V^+   V^+ j+    V^+ )
    !
    ! So, at least when dealing with real sp matrix elements, if we naively
    ! calculate with things that are stored we get
    !
    !               V_stored^\dagger j+  V_stored^* = V^+ j+   V^+
    !         or    V_stored^\dagger j+- V_stored^* = V^+ j+-  V^+
    !
    ! which is negative of what we should get!
    !
    ! (The first term involving U^\dagger and U do not have such pitfalls, they
    !  are simply block diagonal.) 
    ! 
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in)  :: j(nwt,nwt)
    real(KIND=dp), intent(in)  :: bogo(2*nwt, 2*nwt)
    real(KIND=dp), intent(out) :: j11(nwt,nwt)
    real(KIND=dp)              :: tmp(nwt,nwt)
    integer                    :: N, si, sb, b,  T, N2, s

    si = 0 ; sb = 0
    do b = 1, Blocks,2 
      N = HFBlocks_global(b)    ; if(N.eq.0) cycle
      N2 = HFBlocks_global(b+1)

      T = N + N2
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! U^\dagger j 
      j11(si+1:si+T, si+1:si+T) = &
      & matmul(transpose(bogo(sb+1:sb+T, sb+T+1:sb+2*T)),j(si+1:si+T,si+1:si+T))
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! U^\dagger j U
      j11(si+1:si+T, si+1:si+T) = &
      & matmul( j11(si+1:si+T, si+1:si+T), bogo(sb+1:sb+T, sb+T+1:sb+2*T))
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! V^\dagger j^t 
      tmp(si+1:si+T, si+1:si+T) = transpose(j(si+1:si+T, si+1:si+T))
      tmp(si+1:si+T, si+1:si+T) = &
      & matmul(transpose(bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T)),                  &
      &                             tmp(si+1:si+T,si+1:si+T))

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! V^\dagger j^T V
      ! Here we take care of the hidden minus sign as commented above.
$TR   s = -1
$NTR  s = +1
      j11(si+1:si+T, si+1:si+T) = j11(si+1:si+T, si+1:si+T) - s * &
      &  matmul( tmp(si+1:si+T, si+1:si+T), bogo(sb+T+1:sb+2*T, sb+T+1:sb+2*T))

      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
    enddo
  end subroutine calcJ11

  subroutine PrintMomentsofInertia()
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    1 format (21('-'), 'Rotational properties', 20('-'))
    2 format ('    Rigid rotor                       (hbar^2/MeV)')
    3 format ('              n', 14x ,'p',14x,'t ')
    4 format (' I_R X ', 3f15.7)
    5 format (' I_R Y ', 3f15.7)
    6 format (' I_R Z ', 3f15.7)
    7 format ('    Belyaev                          (hbar^2/MeV) ')
   71 format ('    Belyaev (collective)             (hbar^2/MeV) ')
    8 format (' I_B X ', 3f15.7)
    9 format (' I_B Y ', 3f15.7)
   10 format (' I_B Z ', 3f15.7)
   11 format ('     <J^2> - <J>^2                        (hbar^2)')
   12 format ('     <J^2> - <J>^2 (with pairing cutoff)  (hbar^2)')
   13 format (' J2_X  ', 3f15.7)
   14 format (' J2_Y  ', 3f15.7)
   15 format (' J2_Z  ', 3f15.7)
   16 format (' J2_t  ', 3f15.7)
   17 format ('     <J^2> - <J>^2 (collective)           (hbar^2)')

    print 1
    print 2
    print 3
    print 4, Rigid_MOI(1,:)
    print 5, Rigid_MOI(2,:)
    print 6, Rigid_MOI(3,:)
    print *
    print 7
    print 3
    print  8, Belyaev(1,:)
    print  9, Belyaev(2,:)
    print 10, Belyaev(3,:)
    print *
    print 71
    print  3
    print  8, Bely_coll(1,:)
    print  9, Bely_coll(2,:)
    print 10, Bely_coll(3,:)
    print *
    print 11
    print 3
    print 13, J2(1,:)
    print 14, J2(2,:)
    print 15, J2(3,:)
    print 16, sum(J2,1)
    print *
    print 12
    print 3
    print 13, J2_pairing_cut(1,:)
    print 14, J2_pairing_cut(2,:)
    print 15, J2_pairing_cut(3,:)
    print 16, sum(J2_pairing_cut,1)
    print *
    print 17
    print 3
    print 13, J2_coll(1,:)
    print 14, J2_coll(2,:)
    print 15, J2_coll(3,:)
    print 16, sum(J2_coll, 1)

  end subroutine PrintMomentsofInertia

  subroutine select_J2_and_MOI(J2, MOI)
    !-------------------------------------------------------------------------
    ! Selecting the 'right' moment of inertia and <J^2> to be used in the
    ! calculation of the corrections for spurious motion.
    !
    ! This depends on the type of pairing we consider.
    !
    ! Input :
    !    NONE
    ! Output:
    !    J2 : the expectation value of J^2 in all three Cartesian directions
    !    MOI: the appropriate Belyaev moment of inertia along all axes
    !-------------------------------------------------------------------------

    real(KIND=dp), intent(out) :: J2(3), MOI(3)

    select case(pairingtype)
    case(0,1)
      ! HF or BCS
      MOI  = Belyaev(:,3)
      J2   = J2_pairing_cut(:,3)
      ! Sanity check: no collective sense of rotational correction implemented
      !               yet for HF/BCStype calculations
      if(blocktype.ne.0) then
        call stp('Rotational correction for odd nuclei not incorporated into HF/BCS.')
      endif
    case (2)
      ! HFB
      MOI  = Bely_coll(:,3)
      J2   = J2_coll(:,3)   ! Take the collective value.
    end select
  end subroutine select_J2_and_MOI
end module momentsofinertia
