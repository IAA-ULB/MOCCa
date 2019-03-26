module temperature_projection
!-------------------------------------------------------------------------------
!
!
!
!
!
!
!-------------------------------------------------------------------------------

use geninfo
use HartreeFock
use BCS
use pairing
use functional
use wavefunctions
use moments

implicit none

    !---------------------------------------------------------------------------
    ! Contains the natural logarithms of the various partition functions.
    !---------------------------------------------------------------------------
    real*16 :: partition, partition_tilde
    real*16 :: projectedpartition, projectedpartition_nov(2)

    real(KIND=dp) :: projection_cutoff = 10

contains

    subroutine ProjectThermal()
        !-----------------------------------------------------------------------
        !
        !
        !
        !
        !-----------------------------------------------------------------------

        select case (PairingType)

        case(0)
            call ProjectThermalHartreeFock
        case(1)
            call ProjectThermalBCS
        case(2)
            call ProjectThermalHFB
        end select

        call PrintThermalProjection
        
    end subroutine ProjectThermal
    
    subroutine ProjectThermalBCS()
        !-----------------------------------------------------------------------
        !
        !
        !-----------------------------------------------------------------------

        complex*32, allocatable :: xi(:)
        complex*32              :: fac, Iimag, temp

        integer :: iphi, wave, i, it
        real*16 :: u, v, eqp, Z(2), phi, Venergy, avn, avp, f, constraints

        type(moment), pointer   :: current

        !-----------------------------------------------------------------------
        ! We estimate <V>
        !
        ! <V>  = <H_HFB> - E
        !      = sum_k>0 E^qp_k f_k - 1/2 * sum_k E^qp_k  
        !                           + 1/2 * Tr (h - mu) - E
        !                           + mu * N 
        !                           
        ! Missing: contribution of constraints
        !-----------------------------------------------------------------------
        Venergy = - totalE
        do i=1,nwt
            it = 1
            if(i.gt.nwn) it = 2

            !  All of the factors 2 due to time-reversal
            f  = 1.0/(1+exp(inversetemp* BCSqps(i)))
            Venergy = Venergy + 2 * BCSqps(i) * (f - 0.5)                      &
            &                 + 2.0 * 0.5 * (spenergies(i) - FermiEnergy(it))                  
        enddo
        Venergy = Venergy + FermiEnergy(1) * neutrons + FermiEnergy(2) * protons
  
        constraints = 0
        Current =>Root
        do while(associated(Current%next)) 
            Current => Current%next
            if(Current%constrainttype .ne.0) then
              constraints = constraints + Current%multiplier*sum(Current%value)
            endif
        enddo

        Venergy = Venergy + constraints

        Iimag = dcmplx(0,1.0)  
        !-----------------------------------------------------------------------
        ! We calculate first the unprojected partition function
        ! lnZ = -Beta * E_HF - beta * mu * <N> 
        partition = - inversetemp    *  totalE                                 & 
        &           + FermiEnergy(1) * neutrons                                &
        &           + FermiEnergy(2) * protons                                 & 
        &           + sum(entropy)   + constraints

        partition_tilde = - inversetemp * totalE   + sum(entropy) + constraints
        !-----------------------------------------------------------------------
        ! Neutrons
        allocate(xi(2*nwn)) ; xi =0.0

        do iphi = 1, 2*nwn
          phi = pi * (iphi-1)/nwn
          
          ! Calculate xi^BCS for every n
          do wave=1,nwn
              
              u   = 1 - BCSoccupations(wave)/2.0
              v   =     BCSoccupations(wave)/2.0
              eqp = BCSqps(wave)

              fac = u     
              fac = fac +     exp(  2 * Iimag * phi)    * v 
              fac = fac + 2 * exp( -  inversetemp*eqp + Iimag * phi)
              fac = fac +     exp( -2*inversetemp*eqp) * (v+exp(2*Iimag*phi)* u)

              temp = inversetemp * Eqp
              temp = temp - inversetemp * ( spenergies(wave) - FermiEnergy(1) )
              temp = temp + log(fac)
 
              xi(iphi) = xi(iphi) + temp 
          enddo
          xi(iphi) = xi(iphi) - Iimag * phi * neutrons
        enddo
     
        avn = sum(abs(xi))/(2*nwn)
        xi  = xi - avn

        Z(1) = avn + DBLE(log(sum(exp(xi))))
        Z(1) = Z(1) - inversetemp * FermiEnergy(1) * neutrons - log(2*nwn * 1.0)    

        deallocate(xi)
        !-----------------------------------------------------------------------
        ! Protons
        allocate(xi(2*nwp)) ; xi = 0.0

        do iphi = 1, 2*nwp
          phi = pi * (iphi-1)/nwp
          
          do wave=nwn+1,nwt
              u   = 1 - BCSoccupations(wave)/2
              v   =     BCSoccupations(wave)/2
              eqp = BCSqps(wave)

              fac = u     
              fac = fac +     exp(  2 * Iimag * phi)   * v 
              fac = fac + 2 * exp( -  inversetemp*eqp  + Iimag * phi)
              fac = fac +     exp( -2*inversetemp*eqp) *(v + exp(2*Iimag*phi)*u)

              temp = inversetemp * Eqp
              temp = temp - inversetemp * ( spenergies(wave) - FermiEnergy(2) )
              temp = temp + log(fac)
 
              xi(iphi) = xi(iphi) + temp 
          enddo
          xi(iphi) = xi(iphi) - Iimag * phi * protons 
        enddo

        avp = sum(abs(xi))/(2*nwp)
        xi  = xi - avp

        Z(2) = avp + DBLE(log(sum(exp(xi))))
        Z(2) = Z(2) - inversetemp * FermiEnergy(2) * protons - log(2*nwp * 1.0)
        !-----------------------------------------------------------------------
  
        projectedpartition_nov = Z
        projectedpartition     = sum(Z) + inversetemp * Venergy

    end subroutine ProjectThermalBCS    

    subroutine ProjectThermalHartreeFock()
        !-----------------------------------------------------------------------
        !
        !-----------------------------------------------------------------------
    
        integer          :: i, iphi, wave
        complex*32       :: terms_n(2*nwn), terms_p(2*nwp)
        real*16          :: Venergy, fac, Z(2), avn, avp, constraints
        complex*32       :: temp, Iimag
        complex*32       ::  phi
        type(moment), pointer  :: Current

        Iimag = dcmplx(0,1.0)

        constraints = 0
        Current =>Root
        do while(associated(Current%next)) 
            Current => Current%next
            if(Current%constrainttype .ne.0) then
              constraints = constraints + Current%multiplier*sum(Current%value)
            endif
        enddo
        !-----------------------------------------------------------------------
        ! We calculate first the unprojected partition function
        ! lnZ = -Beta * E_HF - beta * mu * <N> 
        partition = - inversetemp    *  totalE                                 & 
        &           + FermiEnergy(1) * neutrons                                &
        &           + FermiEnergy(2) * protons                                 & 
        &           + sum(entropy) + constraints

        partition_tilde = - inversetemp    *  totalE  + sum(entropy)           &
        &                                             + constraints
        !-----------------------------------------------------------------------
        ! We estimate <V>
        Venergy = - totalE
        do i=1,nwt
            Venergy = Venergy +  spenergies(i) * rho_can(i)
        enddo
        Venergy = Venergy + constraints

        Z = 0
        !-----------------------------------------------------------------------       
        ! We sum the logarithms of all the factors for numerical stability.  
        do iphi=1,2*nwn
            phi = pi*(iphi-1)/nwn
            temp = 0                        
            do wave=1,nwn
               fac  = inversetemp * (spenergies(wave) - fermienergy(1))
               temp = temp + log(1+exp(-fac + Iimag*phi))  
            enddo
            !  Factor 2 is time-reversal
            terms_n(iphi) =  2 * temp + (-Iimag * phi * neutrons) 
        enddo
        
        do iphi=1,2*nwp
            phi = pi*(iphi-1)/nwp
            temp = 0            
            do wave=nwn+1,nwt
               fac  = inversetemp * (spenergies(wave) - fermienergy(2))
               temp = temp + log(1+exp(-fac + Iimag*phi))  
            enddo
            terms_p(iphi) =  2 * temp + (-Iimag * phi * protons) 
        enddo

        ! We subtract the average of all logarithms from the terms....
        avn = sum(abs(terms_n))/(2*nwn)
        avp = sum(abs(terms_p))/(2*nwp)
      
        terms_n = terms_n - avn
        terms_p = terms_p - avp        
          
        ! ... and add it again at the end
        Z(1) = DBLE(log(sum(exp(terms_n)))) + avn
        Z(2) = DBLE(log(sum(exp(terms_p)))) + avp
    
        !-----------------------------------------------------------------------
        !  Further, easy corrections
        Z(1) = Z(1) - inversetemp * FermiEnergy(1)*neutrons - log(2*nwn * 1.0)
        Z(2) = Z(2) - inversetemp * FermiEnergy(2)*protons  - log(2*nwp * 1.0)

        projectedpartition_nov = Z
        projectedpartition     = sum(Z) + inversetemp * Venergy

    end subroutine ProjectThermalHartreeFock
    
    subroutine ProjectThermalHFB()
      !-------------------------------------------------------------------------
      ! Note that this routine relies heavily on timereversal.
      !-------------------------------------------------------------------------
      integer                 :: i,j, it
      real(KIND=dp)           :: Venergy, f, Eqp
      real*16, allocatable    :: neutron_terms(:), proton_terms(:)
      real*16                 :: avn, avp, Z(2)
      !-------------------------------------------------------------------------
      ! We calculate first the unprojected partition function
      ! lnZ = -Beta * E_HF - beta * mu * <N> 
      partition = - inversetemp    *  totalE                                   & 
      &           + FermiEnergy(1) * neutrons                                  &
      &           + FermiEnergy(2) * protons                                   & 
      &           + sum(entropy)
      partition_tilde = - inversetemp    *  totalE  + sum(entropy)
      !-------------------------------------------------------------------------
      ! We estimate <V>
      !
      ! <V>  = <H_HFB> - E
      !      = sum_k>0 E^qp_k f_k - 1/2 * sum_k E^qp_k  
      !                           + 1/2 * Tr (h - mu) - E
      !                           + mu * N 
      Venergy = - totalE
      do i=1,nwt
          it = 1
          if(i.gt.nwn) it = 2

          !  All of the factors 2 are due to time-reversal
          Eqp = qpenergies(i)
          f  = 1.0/(1+exp(inversetemp*Eqp))
          Venergy=Venergy+ 2*(Eqp*(f-0.5) + 0.5*(spenergies(i)-FermiEnergy(it)))                  
      enddo
      Venergy = Venergy + FermiEnergy(1) * neutrons + FermiEnergy(2) * protons

      !-------------------------------------------------------------------------
      neutron_terms = &
      &        HFBdeterminant(Bogoliubov(      1:2*nwn,      1:2*nwn),neutrons,&
      &                                        HFBsizes(1:4), Qpenergies(1:nwn))
      proton_terms = &
      &        HFBdeterminant(Bogoliubov(2*nwn+1:2*nwt,2*nwn+1:2*nwt),protons, &
      &                                   HFBsizes(5:8),  Qpenergies(nwn+1:nwt))

      !-------------------------------------------------------------------------
      avn = sum(neutron_terms)/(2*nwn)
      avp = sum(proton_terms)/(2*nwp)

      neutron_terms = neutron_terms - avn
      proton_terms  = proton_terms  - avp
    
      Z = 0
      do i=1,2*nwn
        Z(1) = Z(1) + (-1)**(i-1) * exp(neutron_terms(i))
      enddo
      do i=1,2*nwp
        Z(2) = Z(2) + (-1)**(i-1) * exp(proton_terms(i))
      enddo
      Z(1) = log(Z(1)) + avn 
      Z(2) = log(Z(2)) + avp
      !-------------------------------------------------------------------------
      do i=1,nwt
        it = 1
        if(i.gt. nwn) it = 2
        Z(it) = Z(it) - inversetemp * (spenergies(i) - FermiEnergy(it)) 
      enddo
      !-------------------------------------------------------------------------      
      Z(1) = Z(1) - inversetemp * FermiEnergy(1)*neutrons - log(2*nwn * 1.0)
      Z(2) = Z(2) - inversetemp * FermiEnergy(2)*protons  - log(2*nwp * 1.0)

      projectedpartition_nov = Z
      projectedpartition     = sum(Z) + inversetemp * Venergy

    end subroutine ProjectThermalHFB

    function HFBdeterminant(Bogo, particles, blocks,Eqp) result(detMR)
        !-----------------------------------------------------------------------
        !
        !-----------------------------------------------------------------------
        real(KIND=dp), intent(in) :: Bogo(:,:), Eqp(:)
        real(KIND=dp), intent(in) :: particles
        real(KIND=dp)             :: phi
        integer, intent(in)       :: blocks(4)

        !-----------------------------------------------------------------------
        ! Attention: this routine does not operate at the full possible 
        ! precision for complex numbers (complex*32). The reason is that 
        ! Lapack does not provide a subroutine for the QR decomposition. 
        ! Feeding larger precision complex numbers into LAPACK routines gives
        ! nonsense results.
        !-----------------------------------------------------------------------
        complex*16, allocatable :: matrix(:,:), tau(:,:), work(:),detM(:)  
        complex*16              :: Iimag
        real*16, allocatable    :: detMR(:)
        integer                 :: iphi, N, Ntotal, sb, si, B, info, i, it

        Iimag = cmplx(0, 1.0)

        Ntotal = sum(blocks)
        allocate(detM(2*Ntotal)) ; detM = 0

        do iphi=1,2*Ntotal
          sb = 0 ; si = 0
          phi = pi * (iphi - 1)/Ntotal

          detM(iphi) = 0
          ! We calculate the logarithm of the derminant of every subblock
          do B = 1,4
            N = blocks(B) ; if (N .eq. 0) cycle

            !-------------------------------------------------------------------
            ! Building the matrix in this symmetry block
            allocate(matrix(2*N, 2*N)) ; matrix = 0        
            do i=1,N
              !M  = e^{- i phi N}
              matrix(i  ,i  ) = exp(-Iimag * phi)
              matrix(i+N,i+N) = exp( Iimag * phi)
            enddo

            matrix = matmul(           matrix,(Bogo(sb+1:sb+2*N, sb+1:sb+2*N))) 
            matrix = matmul(transpose(Bogo(sb+1:sb+2*N, sb+1:sb+2*N)), matrix)

            do i=1,N
              ! M = W^dagger e^{-i phi N} W + e^{-\beta E}
              matrix(i  ,i  ) = matrix(i  ,i  ) + exp( -inversetemp * Eqp(si+i))
              matrix(i+N,i+N) = matrix(i+N,i+N) + exp(  inversetemp * Eqp(si+i))
            enddo    

            !-------------------------------------------------------------------
            !  Performing a QR decomposition
            allocate(tau(2*N,2*N), work(4*N)) 
            call ZGEQRF(2*N, 2*N, matrix, 2*N, tau, work, 4*N,info)

            if(info.ne.0) then
              print *, 'ZGERQF failed. Error=' , info
              stop
            endif

            ! Logarithm of the determinant
            do i=1,2*N
              detM(iphi) = detM(iphi) + log(matrix(i,i))
            enddo
            deallocate(matrix, work, tau)
            sb = sb + 2*N
            si = si +   N
          enddo
          detM(iphi) = detM(iphi) - Iimag * phi * particles
        enddo
        !-----------------------------------------------------------------------
        allocate(detMR(2*Ntotal))
        detMR = DBLE(detM)

    end function HFBdeterminant

    subroutine PrintThermalProjection
        !-----------------------------------------------------------------------
        !
        !
        !
        !-----------------------------------------------------------------------

        1 format (80('-'))
        2 format (' Thermal particle number projection')
        3 format (' Partition function')
        4 format (' lnZ (constant mu)  : ', f20.12)
        5 format (' lnZ (constant  N)  : ', f20.12)
!        5 format (' Projected  lnZ (No V) (n): ', f20.12,/, &
!        &         '                       (p): ', f20.12,/, &
!        &         '                       (t): ', f20.12)
        6 format (' lnZ (canonical)    : ', f20.12)
       61 format (' Full precision     : ', f20.12) 



        print 1
        print 2
        print 3
        print 4, partition
        print 5, partition_tilde
!        print 5, projectedpartition_nov, sum(projectedpartition_nov)
        print 6, projectedpartition
        print 61, TotalE
        print 1
        
    end subroutine PrintThermalProjection

    subroutine QRDeterminant
        !-----------------------------------------------------------------------
        ! Use a QR-decomposition to calculate
        !
        !   log(det(A))
        ! 
        ! for A a complex matrix, in a stable fashion. 
        !
        ! We calculate the
        !-----------------------------------------------------------------------

    end subroutine 
end module temperature_projection
