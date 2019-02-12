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
    
    subroutine ProjectThermalHFB()
      !-------------------------------------------------------------------------
      !
      !
      !
      ! Note that this routine relies heavily on timereversal.
      !-------------------------------------------------------------------------
      
      complex*16, allocatable :: matrix(:,:), A(:,:), Tau(:,:), work(:)
      real*16, allocatable    :: scales(:)
      complex*16              :: Iimag, phi, detM, detA
      integer                 :: i,j,si,sb, B, N, info

      allocate(matrix(2*nwt, 2*nwt))

      Iimag = cmplx(0, 1.0)

      si = 0 ; sb = 0
      !-------------------------------------------------------------------------
      ! a) Building the matrix of which we need the determinant
      do B=1,8
          N = HFBsizes(B) ;  if (N.eq.0) cycle

          do i=1,N
              !M  = e^{- i phi N}
              matrix(sb+i  ,sb+i  ) = exp(-Iimag * phi)
              matrix(sb+i+N,sb+i+N) = exp( Iimag * phi)
          enddo
          ! M = W^\dagger e^{- i phi N} W
          matrix(sb+1:sb+2*N, sb+1:sb+2*N) = matmul(                           &
          &                                   matrix(sb+1:sb+2*N, sb+1:sb+2*N),&
          &                                Bogoliubov(sb+1:sb+2*N, sb+1:sb+2*N)) 

          matrix(sb+1:sb+2*N, sb+1:sb+2*N) = matmul(                           &
          &                    transpose(Bogoliubov(sb+1:sb+2*N, sb+1:sb+2*N)),&
          &                                    matrix(sb+1:sb+2*N, sb+1:sb+2*N))


          do i=1,N
              ! M = W^dagger e^{-i phi N} W + e^{-\beta E}
              matrix(sb+i  ,sb+i  ) = matrix(sb+i  ,sb+i  )                    &
              &                          + exp( -inversetemp * qpenergies(si+i))
              matrix(sb+i+N,sb+i+N) = matrix(sb+i  ,sb+i  )                    &
              &                          + exp(  inversetemp * qpenergies(si+i))
          enddo
          si = si +  N
          sb = sb +2*N
      enddo
      !------------------------------------------------------------------------- 
      ! b) Piece-wise QR decomposition with the zgerqf lapack routine
      si = 0 ; sb = 0

      detM = 0.0
      do B=1,8
         N = HFBsizes(B) ;  if (N.eq.0) cycle
         
         allocate(A(2*N,2*N), tau(2*N,2*N), work(2*N), scales(2*N)) 
         A = matrix(sb+1:sb+2*N,sb+1:sb+2*N)
         call ZGERQF(2*N, 2*N,A, 2*N, tau, work, 2*N,info)

         if(info.ne.0) then
            print *, 'ZGERQF failed. Error code = ', info
            stop
         endif

         !  Extract the scales of the problem
         do i=1,2*N
            scales(i)  = maxval(abs(A(i,1:2*N)))
            A(i,1:2*N) = A(i,1:2*N)/scales(i)
         enddo

         ! Compute the logarithm of the determinant
         do i=1,2*N
            detM = detM + log(scales(i)) + log(A(i,i))
         enddo
         deallocate(A, tau, work, scales)
      enddo
      
      deallocate(matrix)      
    end subroutine ProjectThermalHFB

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

        !-----------------------------------------------------------------------
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
