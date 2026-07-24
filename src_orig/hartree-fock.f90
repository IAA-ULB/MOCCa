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
module hartreefock
 !==============================================================================
 !
 ! Module containing all of the necessary routines to correctly decide on the 
 ! occupation of single-particle wave-functions in the case of a Hartree-Fock
 ! calculation.
 !
 !
 !------------------------------------------------------------------------------
 ! Hephaestos keywords:
 ! 
 ! Switching between conservation of an antilinear, antihermitian symmetry
 ! TR : $TR
 ! NTR: $NTR
 !==============================================================================
 
 use vectors, only: DensityVector, PotentialVector
 use wavefunctions
 
 implicit none

 real(KIND=dp) :: HFdispersion(2) = 0.0
 real(KIND=dp) :: FermiEnergyHF(2) = 0.0
 
contains
 
 subroutine NaiveFill(occupations)
    !---------------------------------------------------------------------------
    ! This subroutine finds the orbitals with the lowest single particle
    ! energy and fills them, after sorting all the levels.
    !
    ! At zero temperature, the occupations are simply set to be one or zero.
    ! Note that this does not necessarily play well near spherical symmetry for
    ! non-magic numbers of nucleons.
    !---------------------------------------------------------------------------

    integer :: i,j,jp1,n,p, ProtonUpperBound, NeutronUpperBound
    integer :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp), intent(out) :: occupations(nwt)

    n=0; p=0
    !---------------------------------------------------------------------------
    !Setting all occupation numbers to Zero
    Occupations = 0.0_dp
    !---------------------------------------------------------------------------
    ProtonUpperBound  = floor(Protons)
    NeutronUpperBound = floor(Neutrons)
    !---------------------------------------------------------------------------
    !Finding the order of the spwfs, in terms of energy
    ProtonOrder = OrderSpwfsISO(+1)
    NeutronOrder= OrderSpwfsISO(-1)
    !---------------------------------------------------------------------------
    !Filling in the lowest Proton orbitals. This is easy, since we know
    !the order of Spwfs.
    i=1 
    do while(p.lt.ProtonUpperBound .and. i.le.nwp)
      j              = ProtonOrder(i)
      $TR  Occupations(j) = 2.0_dp
      $NTR Occupations(j) = 1.0_dp
      p = p + int(occupations(j))
      i = i + 1
    enddo
    !NS: defining FermiEnergy even in HF case as right between last occupied
    ! and the first unoccupied level.
    if(ProtonUpperBound .gt. 0) then
        jp1=ProtonOrder(i)
        FermiEnergyHF(2)=spenergies(j)+(spenergies(jp1)-spenergies(j))/2.d0
    endif
    
    i=1
    do while(n.lt.NeutronUpperBound .and. i.le.nwn)
      j              = NeutronOrder(i)
      $TR  Occupations(j) = 2.0_dp
      $NTR Occupations(j) = 1.0_dp
      n = n + int(occupations(j))
      i = i + 1
    enddo
    if(NeutronUpperBound .gt. 0) then
        jp1=NeutronOrder(i)
        FermiEnergyHF(1)=spenergies(j)+(spenergies(jp1)-spenergies(j))/2.d0
    endif
        
    return
  end subroutine NaiveFill

  subroutine FiniteTemperatureHF(occupations, Fermi, particles_in_gas)
    !---------------------------------------------------------------------------
    ! Simple bisection routine to find the correct Fermi energy for a finite
    ! temperature HF calculation. A Newton type method might be more efficient, 
    ! but I fear for rounding errors and instability with exp() of either large
    ! or small numbers.
    !---------------------------------------------------------------------------

    use parameterization, only : hbm

    real(KIND=dp), intent(out) :: occupations(nwt)
    integer, intent(in)        :: particles_in_gas
    real(KIND=dp)              :: energies(nwt,2), N
    real(KIND=dp)              :: Fermi(2), Fmin, Fmax, betaE, Nmin, Nmax
    integer                    :: Order(nwt,2), nw, nwn_alt, nwp_alt, nw_alt
    integer                    :: it, i

    nwn_alt = 0 ; nwp_alt = 0
    do i=1,nwn
        if(spenergies(i).lt.0) nwn_alt = nwn_alt+1  
    enddo
    do i=nwn+1,nwt
        if(spenergies(i).lt.0) nwp_alt = nwp_alt+1  
    enddo

    if( fixfermi ) then
        !-----------------------------------------------------------------------
        ! We perform a calculation at fixed chemical potential
        Fermi(1) = mun ; Fermi(2) = mup
    else
        !-----------------------------------------------------------------------
        ! We fix the particle number to <N> = neutrons, <Z> = protons.

        ! Finding the order of the spwfs, in terms of energy
        Order = 0    
        Order(1:nwn,1) = OrderSpwfsISO(-1)
        Order(1:nwp,2) = OrderSpwfsISO(+1)

        do i=1,nwn
            energies(i,1) = spenergies(Order(i,1))
        enddo
        do i=1,nwp
            energies(i,2) = spenergies(Order(i,2)) 
        enddo      

        do it=1,2
            if(it .eq. 1) then
                N = neutrons ; nw = nwn ; nw_alt = nwn_alt
            else 
                N = protons  ; nw = nwp ; nw_alt = nwp_alt
            endif
            !-------------------------------------------------------------------
            ! Establish a search interval
            select case(particles_in_gas)
            case(0,1)
              ! Analytically, these are guaranteed to be lower and upper bounds
              ! if we are looking for a traditional solution.
              Fmin = minval(spenergies) - log(2*nw/N - 1)/inversetemp
              Fmax = maxval(spenergies) - log(2*nw/N - 1)/inversetemp
            case(2)
              ! But the upper bound is no longer good if we are only occupying
              ! bound states.
              Fmin = minval(spenergies) - log(2*nw_alt/N - 1)/inversetemp
              Fmax =                    - log(2*nw_alt/N - 1)/inversetemp
            end select 

            Nmin = FToccupations(Fmin, energies(1:nw,it), particles_in_gas) - N
            Nmax = FToccupations(Fmax, energies(1:nw,it), particles_in_gas) - N

            if(Nmin .gt. 0 .or. Nmax .lt. 0) then
              call stp('The Fermi energy was not correctly bracketed')
            endif

            select case (particles_in_gas)
            case(0,2)
              ! Ordinary case & only bound states case
              Fermi(it) = FermiBisection(Fmin,Nmin,Fmax,Nmax, energies(1:nw,it)& 
              &                                    ,N,particles_in_gas, hbm(it))
            case(1)
              ! We do this in two steps, as
              !  N_nucleus = N_total - N_gas is not a monotonic function at all.
              ! Note that this particular implementation is a bit stupid: we 
              ! assume the corrected Fermi energy to be not too far (2 MeV) from
              ! the uncorrected one.
              ! First solve the ordinary problem
              Fermi(it)=FermiBisection(Fmin,Nmin,Fmax,Nmax, energies(1:nw,it), &
              &                        N, 0, hbm(it))

              Fmin = Fermi(it) - 2
              Fmax = Fermi(it) + 2
              Fermi(it) = FermiBisection(Fmin,Nmin,Fmax,Nmax, energies(1:nw,it)& 
              &                                    ,N,particles_in_gas, hbm(it))
            end select
        enddo
    endif
    !---------------------------------------------------------------------------
    ! Actually calculate the occupations for the fixed chemical potential
    do i = 1,nwt
        if (i .gt. nwn) then
            it    = 2 
        else
            it    = 1
        endif        

        select case(particles_in_gas)
        case(0,1)
          betaE = inversetemp * (spenergies(i) - Fermi(it))
          occupations(i) = 2.0d0/(1 + exp(betaE))
        case(2)
          betaE = inversetemp * (spenergies(i) - Fermi(it))
          if(spenergies(i) .gt. 0) then
            occupations(i) = 0  
          else
            occupations(i) = 2.0d0/(1 + exp(betaE))
          endif
        end select
    enddo

    !---------------------------------------------------------------------------
    ! Calculate the HF dispersion
    !  DN = 2  Tr (rho (1-rho))
    !---------------------------------------------------------------------------
    HFdispersion = 0.0
    do i=1,nwt
        if (i .gt. nwn) then
            it    = 2 
        else
            it    = 1
        endif   
        ! Note the time-reversal factors of two
        HFdispersion(it) = HFdispersion(it) +                                  &
        &                        occupations(i)/2.0d0 * (1-occupations(i)/2.0d0)
    enddo
    ! Factor of two from the formula
    HFdispersion = 2 * HFDispersion
    ! Factor two for time-reversal
    HFdispersion = 2 * HFDispersion
  end subroutine FiniteTemperatureHF

  function FToccupations(mu, energies, gas) result (N)
    !---------------------------------------------------------------------------
    ! Simple function that sums the occupations for given mu and sp-energies.
    !---------------------------------------------------------------------------

    integer, intent(in)       :: gas
    real(KIND=dp), intent(in) :: energies(:), mu
    real(KIND=dp)             :: N, betaE
    integer                   :: i

    N = 0

    select case(gas)
    case(0,1)
      do i=1, size(energies)
          betaE = inversetemp * (energies(i) - mu)
          N     = N + 1.0d0/(1 + exp(betaE))
      enddo
    case(2) 
      ! Only take into account bound states
      do i=1, size(energies)
          betaE = inversetemp * (energies(i) - mu)
          if(energies(i) .lt. 0) then 
            N     = N + 1.0d0/(1 + exp(betaE))
          endif
      enddo
    end select
    ! Time-reversal gives the factor 2
    N = 2 * N

  end function FToccupations
  
  recursive function FermiBisection(xmin,Nmin,xmax,Nmax, energies, N,gas,hbm)  &
  &                  result(x)
    !---------------------------------------------------------------------------
    ! Bisection search for an appropriate chemical potential for the FT HF.
    !---------------------------------------------------------------------------

    real(KIND=dp)             :: xmin, xmax, x, Nmin, Nmax, xnew, Nnew, Ngas
    real(KIND=dp), intent(in) :: energies(:), N, hbm  
    integer, intent(in)       :: gas

    ! We allow a larger error here; in my experience there are calculations 
    ! where this routine keeps on calling itself because differences get
    ! smaller than machine precision...
    if(abs(Nmin) .lt. 10 * pairing_prec) then 
        x = xmin        
        return 
    endif
    if(abs(Nmax) .lt. 10 * pairing_prec) then
        x = xmax
        return
    endif   
    
    xnew = 0.5d0 * (xmin + xmax)
    Nnew = FToccupations(xnew, energies, gas)

    !---------------------------------------------------------------------------
    ! We correct for the number of particles suspended in the gas around the
    ! (compound) nucleus.
    select case(gas)
    case(0,2)
      !  No correction
    case(1)
      !  Subtraction method
      Ngas = gasoccupations(xnew, hbm)
      Nnew = Nnew - Ngas
    end select
    Nnew =  Nnew - N

    if(Nnew .gt. 0) then
        x = FermiBisection(xmin, Nmin, xnew, Nnew, energies,N, gas, hbm)
    else
        x = FermiBisection(xnew, Nnew, xmax, Nmax, energies,N, gas, hbm)            
    endif

  end function FermiBisection

  function gasoccupations(fermi, hbm) result(s)
    !---------------------------------------------------------------------------
    ! We count the number of particles in the free gas.
    !
    ! For Lagrange derivatives in a 1D box, the single-particle energies are
    !
    !    epsilon_kx = hbar^2/(2*m) (2*pi * kx/(2*nx * dx))**2
    !
    ! where nx is the size of the EV8-box and 
    !
    !    kx = +/- 1/2, +/-3/2, ....
    !
    ! and we should not forget the spin degree of freedom, resulting in a
    ! overall degeneracy of two!
    !
    ! The levels in the full (3D) box can thus be indexed as
    ! 
    !    psi_kxkykz with a single-particle energy of 
    !     eps_kxkykz  = (eps_kx + eps_ky + eps_kz) 
    !
    ! Since the box-size is not necessarily the same in all dimensions, we 
    ! cannot profit from the degeneracy in the levels for exactly cubic 
    ! boxes. 
    !-----------------------------------------------------------------------
    real(KIND=dp), intent(in) :: hbm, fermi
    real(KIND=dp)             :: s, maxe,ez,ey,ex,checkz, checky, checkx, etot
    real(KIND=dp)             :: trash
    integer                   :: kx, ky, kz

    s = 0
    maxe = maxval(spenergies)

    ! trash statement to stop the compiler complaining
    trash = fermi

    kz=1
    do while(.true.)
      ez = hbm*(pi*kz/(2*nz*dx))**2

      checkz = 1.d0/(1.0d0 + exp(inversetemp*(ez - fermi))) 
      if(checkz.lt. 1d-8) exit
 
      ky = 1
      do while(.true.)
        ey     = hbm*(pi*ky/(2*ny*dx))**2 

        checky = 1.d0/(1.0d0 + exp(inversetemp*(ey - fermi))) 
        if(checky.lt. 1d-8) exit

        kx = 1
        do while(.true.) 
          ex = hbm*(pi*kx/(2*nx*dx))**2 

          checkx = 1.d0/(1.0d0 + exp(inversetemp*(ex - fermi))) 
          if(checkx .lt. 1d-8) exit

          etot = ex + ey + ez
          if(etot .gt. maxe) then
            ! There is an implicit energy cutoff in our calculation:
            ! the highest single-particle energy of the spwfs we consider.
            exit
          endif          
          s  = s + 1.d0/(1.0d0 + exp(inversetemp*(etot - fermi))) 
          kx = kx + 1
        enddo
        ky = ky + 1
      enddo
      kz = kz + 1
    enddo   
    ! Spin degree of freedom
    s= 2*s
  end function gasoccupations

  subroutine CalcHFgaps(Fermi, stabfactor, F)
    !---------------------------------------------------------------------------
    ! Dummy routine.
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)         :: Fermi(2), stabfactor(2)
    type(PotentialVector), intent(in) :: F
    real(KIND=dp)                     :: trash(2) 
    ! trash statements to stop the compiler complaining about unused arguments
    trash = fermi ; trash=stabfactor ; trash = F%F_I_I(1,1)

  end subroutine calcHFgaps

end module hartreefock
