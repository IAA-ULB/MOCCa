module functional
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
 !
 ! Module containing the means to calculate (and print) the mean-field energy.
 ! Note that the actual coupling constants are contained in the constants.f90
 ! file. 
 !==============================================================================
 
 use compilation
 use geninfo
 use densities
 

 implicit none
 
  
    real(KIND=dp) :: t0=-2488.913 
    real(KIND=dp) :: x0=0.834
    real(KIND=dp) :: t1=486.818
    real(KIND=dp) :: x1=-0.344
    real(KIND=dp) :: t2=-546.395
    real(KIND=dp) :: x2=-1.0
    real(KIND=dp) :: t3=13777.0
    real(KIND=dp) :: x3=1.354
    real(KIND=dp) :: yt3a= 0.166666666666666666667 
    real(KIND=dp) :: te=0.0
    real(KIND=dp) :: to=0.0
    real(KIND=dp) :: wso=123.0
    real(KIND=dp) :: wsoq=123.0
    real(KIND=dp) :: t1n2=24.3409
    real(KIND=dp) :: t2n2=-27.31975
    real(KIND=dp) :: x1n2=-0.344
    real(KIND=dp) :: x2n2=-1.0   
    real(KIND=dp) :: hbm(2)         = 20.73551910_dp 
         
    integer       :: COM1body = 2
    real(KIND=dp) :: nucleonmass(2) = (/939.565379_dp , 938.272046_dp /)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Energies are calculated in the BFH representation using the B coupling 
    ! coefficients, and then recombined into the isospin representation using 
    ! the C coefficients. By default only the latter is printed, but the BFH
    ! representation can be asked for for debugging purposes. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    real(KIND=dp) :: Kinetic(2), Skyrme(2), TotalE, SpwfEnergy, COMCorrection(2,2)

    ! Declaration of the energy terms and the coupling coefficients 
$DECLARATION
 

 contains
 
 subroutine calcedfcoefs()
     
$CALCCOEF
    
    nucleonmass = 0.5*(nucleonmass(1) + nucleonmass(2))
 
 end subroutine calcedfcoefs
 
 subroutine printedfcoefs
    !---------------------------------------------------------------------------
    ! Print the values of the EFD coefs used.
   
    1 format (90('-'))
    2 format (' Skyrme coupling constants ')
    3 format (26x, 'Isospin representation BFH representation')
    4 format (26x, 'scalar      vector      total       (n-p)')
    
    print 1
    print 2
    print 3
    print 4
    print 1
$PRINTCOEF    
    print 1
 end subroutine printedfcoefs

 subroutine PrintEnergy()
    1 format (90('-'))
    5 format (30x, '      neutron        proton         total')
    6 format (15x, 'Kinetic Energy:', 3f15.6)
   61 format (15x, '    COM 1-body:', 3f15.6)
    7 format (15x, '  Total energy:', 30x, f15.6)
    8 format (15x, '    from spwfs:', 30x, f15.6)

    call printSkyrme

    print 1
    print 5
    print 6, Kinetic, sum(Kinetic)
    print 61, COMcorrection(:,1), sum(COMcorrection(:,1))
    print 7, TotalE
    print 8, spwfenergy
    print 1
 end subroutine PrintEnergy
 
 subroutine CalcEnergy
    !---------------------------------------------------------------------------
    ! Calculate all of the relevant energies.
    !---------------------------------------------------------------------------
    
    ! Kinetic energy
    Kinetic = CompKinetic()
    ! COM correction
    call CompCOMCorrection()
    ! Skyrme functional
    call compSkyrme()
    ! Total energy from single-particle energies
    SpwfEnergy = calcspwfenergy()
    
 end subroutine CalcEnergy
 
 subroutine CompSkyrme()
    !---------------------------------------------------------------------------
    !
    ! 
    !
    !
    !
    integer       :: it,m,n,k
    real(KIND=dp) :: Edensity(mv,3)
    
$CALCULATION
    
    Skyrme = &
$TOTAL
    
    TotalE = sum(Skyrme + Kinetic) + sum(COMCorrection)
    
 end subroutine CompSkyrme
 
 subroutine PrintSkyrme()
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Print all contributions to the energy.
    !
    
    1 format (90('-'))
    2 format (' Skyrme Energy ')
    3 format (35x, 'isoscalar      isovector      total')
    4 format (17x, 'Total Skyrme:', 3f15.6)
    
    print 1
    print 2
    print 3
    print 1
    
$PRINT 

    print 1
    print 4, Skyrme, sum(Skyrme)
    print 1
 end subroutine PrintSkyrme

 function CompKinetic() result(kinetic)
    !---------------------------------------------------------------------------
    ! This subroutine computes the total kinetic energy,
    ! according to the following formula:
    !    E_k = -\hbar/2m \int d^3x \sum_{k} v_{k} \Psi_k^* \Delta \Psi_k
    !---------------------------------------------------------------------------
    ! Note that the 1-body c.o.m. correction is not taken into account here!
    !---------------------------------------------------------------------------
    use Constants

    integer          :: wave, it,k,i
    real(KIND=dp)    :: Inproduct
    real(KIND=dp)    :: Kinetic(2)

    ! Kinetic Energy
    Kinetic = 0.0_dp
    do wave=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave.le.sum(HFBlocks(1:Blocks/2))) it = 1

        Inproduct = 0.0_dp
        do k=1,4          
                do i=1,mv
                       Inproduct = Inproduct + HFPsi(i,1,1,k,wave) *  & 
                       &  ( HFddPsi(i,1,1,1,1,k,wave) + &
                       &    HFddPsi(i,1,1,2,2,k,wave) + &
                       &    HFddPsi(i,1,1,3,3,k,wave))
                enddo
        enddo
        Kinetic(it)= Kinetic(it) + Occupations(wave)*Inproduct
    enddo
    Kinetic=-Kinetic * hbm * dv
    return
  end function CompKinetic
  
  subroutine CompCOMCorrection()
    !---------------------------------------------------------------------------
    ! M. Bender et al., Eur. Phys. J. A 7, 467-478 (2000)
    !---------------------------------------------------------------------------
   
    COMCorrection = 0.0_dp
    if(COM1Body .gt. 0) then
      !Deduce 1-body COM correction from the Kinetic Energy
      COMCorrection(1,:) = - Kinetic(:) * nucleonmass/ &
      &                 (neutrons * nucleonmass(1) + protons * nucleonmass(2))
    endif

  end subroutine CompCOMCorrection

  subroutine calcFields()
        
        integer :: m, n, k, it

$CALCFIELDS


  end subroutine calcFields 
  
  function sphamil(psi, dpsi, ddpsi, sx,sy,sz,iso) result(hpsi)
    !---------------------------------------------------------------------------
    ! Apply the action of the single-particle hamiltonian to the 
    ! single-particle wave-functions.
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in) :: psi(nx,ny,nz,4),     dpsi(nx,ny,nz,3,4)
    real(KIND=dp), intent(in) :: ddpsi(nx,ny,nz,3,3,4)
    integer, intent(in)       :: sx(4),sy(4),sz(4),   iso
    real(KIND=dp)             :: hpsi(nx,ny,nz,4)
    real(KIND=dp)             :: temp(nx,ny,nz,4), dtemp(nx,ny,nz,3,4)
    real(KIND=dp)             :: ddtemp(nx,ny,nz,3,3,4), laptemp(nx,ny,nz,4)
    
    real(KIND=dp)             :: ReducedMass
    
    integer :: it, i,k
    
    !---------------------------------------------------------------------------
    ! Determine the isospin index
    it = (iso + 3)/2
    
    ! Reduced mass in case of self-consistent 1-body COM correction
    if(COM1body .eq. 2) then
        Reducedmass = (1.0_dp-nucleonmass(it)/                                 &
        &                      (neutrons*nucleonmass(1)+protons*nucleonmass(2)))
    else
        Reducedmass = 1.0_dp
    endif
    
    !---------------------------------------------------------------------------
    ! Action of the kinetic energy
    do k=1,4
        do i=1,mv
            hpsi(i,1,1,k) = - hbm(it)* reducedmass *       (ddpsi(i,1,1,1,1,k) &
            &                                             + ddpsi(i,1,1,2,2,k) &
            &                                             + ddpsi(i,1,1,3,3,k)) 
        enddo
    enddo
    !---------------------------------------------------------------------------
    ! Action of the Skyrme fields
    !
    ! Note that 
    ! a) Coulomb is included in the F_I_I field
    ! b) Every density contains the contributions from constraints on that 
    !    density. 
    ! c) The kinetic energy is NOT included in the F_N_N field, because 
    !    a constant is not in the Lagrange basis; so the current way of
    !    deriving stuff is not correct for a term
    !          hbar^2_2m
    !---------------------------------------------------------------------------
    
$SKYRMEACTION
    
  end function sphamil
  
  function calcspwfenergy() result(spwfenergy)
    !---------------------------------------------------------------------------
    ! Calculates the total energy from the single-particle energies. 
    !
    !---------------------------------------------------------------------------
    
    use wavefunctions
    
    integer       :: wave
    real(KIND=dp) :: spwfenergy, e_rear
    
    ! Start by summing the single-particle energies
    spwfenergy = 0
    do wave=1,nwt
        spwfenergy = spwfenergy + occupations(wave) * spenergies(wave)
    enddo
    
    ! Calculation of rearrangement energy
$EREAR   
    spwfenergy = spwfenergy - e_rear
    
    spwfenergy = 0.5 * spwfenergy + 0.5 * sum(kinetic)
    
    ! Always add the 1-body COMcorrection. In case it is used iteratively, it
    ! is double counted along with the kinetic energy!
    if(COM1body.gt.0) then
    	SpwfEnergy = SpwfEnergy  + sum(COMCorrection(1,:))/2.0_dp
    endif
    
  end function calcspwfenergy
end module functional
