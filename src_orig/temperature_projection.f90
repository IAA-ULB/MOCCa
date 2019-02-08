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

implicit none

    !---------------------------------------------------------------------------
    ! Contains the natural logarithms of the various partition functions.
    !---------------------------------------------------------------------------
    real(KIND=dp) :: partition
    real(KIND=dp) :: projectedpartition, projectedpartition_nov

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
            return
        case(2)
            return
        end select

        call PrintThermalProjection
        
    end subroutine ProjectThermal
    
    subroutine ProjectThermalHartreeFock()
        !-----------------------------------------------------------------------
        !
        !
        !
        !
        !-----------------------------------------------------------------------
    
        integer          :: it, i, gaussn, gaussp, iphi, wave
        real*16          :: V, terms_n(2*nwn), terms_p(2*nwp), Z(2)
        complex*32       :: prod, Iimag
        complex*32       :: fac, N, phi
        
        !-----------------------------------------------------------------------
        ! We calculate first the unprojected partition function
        ! lnZ = -Beta * E_HF - beta * mu * <N> + S_HF
        partition = -inversetemp* (totalE          + FermiEnergy(1) * neutrons &
        &                                          + FermiEnergy(2) * protons )&
        &           + sum(entropy)

        !-----------------------------------------------------------------------
        ! We estimate <V>  as the difference between the total energy and 
        ! the single-particle energies
        V = totalE
        do i=1,nwt
            V = V - spenergies(i) * rho_can(i)
        enddo

        !-----------------------------------------------------------------------
        ! Collocation points of the Fourier sum are equal to the number states
        ! in our modelspace. Notice the factor two for time-reversal 
        gaussn = 2 * nwn
        gaussp = 2 * nwp

        Iimag = dcmplx(0,1.0)
            
        Z = 0
            
        do iphi=1,gaussn
            phi = 2*pi*(iphi-1)/gaussn
    
            prod = dcmplx(1.0,0.0)
            do wave=1,nwn
                fac      = inversetemp*(spenergies(wave) - fermienergy(1))
                prod     = prod    *  (1.0 + exp(-fac)*exp(Iimag * phi))
            enddo
            terms_n(iphi) = (exp(-Iimag * phi * neutrons) * prod**2)/(2*nwn)      
        enddo

        do iphi=1,gaussp
            phi = 2*pi*(iphi-1)/gaussp
    
            prod = dcmplx(1.0,0.0)
            do wave=nwn+1,nwt
                fac      = inversetemp*(spenergies(wave) - fermienergy(2))
                prod     = prod    *  (1.0 + exp(-fac)*exp(Iimag * phi))
            enddo
            terms_n(iphi) = (exp(-Iimag * phi * protons) * prod**2)/(2*nwp)     
        enddo

        Z(1) = log(sum(terms_n))
        Z(2) = log(sum(terms_p))

        Z(1) = Z(1)- inversetemp * FermiEnergy(1)*neutrons 
        Z(2) = Z(2)- inversetemp * FermiEnergy(2)*protons  

        projectedpartition_nov = sum(Z) 
        projectedpartition     = sum(Z) + inversetemp * V
  
    end subroutine ProjectThermalHartreeFock

    subroutine PrintThermalProjection
        !-----------------------------------------------------------------------
        !
        !
        !
        !-----------------------------------------------------------------------

        1 format (80('-'))
        2 format (' Thermal particle number projection')
        3 format (' Partition function')
        4 format (' Ordinary   lnZ       : ', f15.7)
        5 format (' Projected  lnZ (No V): ', f15.7)
        6 format (' Projected  lnZ (   V): ', f15.7)

        print 1
        print 2
        print 3
        print 4, partition
        print 6, projectedpartition
        print 5, projectedpartition_nov
        print 1
        
    end subroutine PrintThermalProjection

end module temperature_projection
