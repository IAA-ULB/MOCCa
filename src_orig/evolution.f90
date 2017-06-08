module evolution
!===============================================================================
!  #######   ##   #    # #####   ##   #      #    #  ####
!     #     #  #  ##   #   #    #  #  #      #    # #
!     #    #    # # #  #   #   #    # #      #    #  ####
!     #    ###### #  # #   #   ###### #      #    #      #
!     #    #    # #   ##   #   #    # #      #    # #    #
!     #    #    # #    #   #   #    # ######  ####   ####
!
!  Copyright W. Ryssens & M. Bender
!
!===============================================================================
!
! Module that governs the evolution of the single-particle wavefunctions from 
! one iteration to the next. 
!
! Currently possible:
!   a) Gradient descent aka imaginary time step.
!
!===============================================================================


    use wavefunctions
    use functional

    implicit none
    
    !---------------------------------------------------------------------------
    ! Parameters of the iteration scheme
    real(KIND=dp):: dt   = 0.06
    real(KIND=dp):: hbar =  6.58211928_dp
    
    !---------------------------------------------------------------------------
    !Procedure that determines the evolution of a Spwf under imaginary time.
    abstract interface
      subroutine Evolve_interface(Iteration)
        integer, intent(in)       :: iteration
      end subroutine
    end interface
    procedure(Evolve_Interface),pointer :: Evolve    
    !---------------------------------------------------------------------------

contains
    
    subroutine Evolve_graddesc(iteration)
        !-----------------------------------------------------------------------
        ! 
        ! a) For every wave-function do a gradient step
        !    
        !    psi => ( 1 - dt/hbar h ) psi
        ! b) Calculate values:
        !    < psi | h   | psi >
        !    < psi | h^2 | psi >
        ! b) Orthonormalize within symmetry blocks
        !-----------------------------------------------------------------------
        
        use wavefunctions
        
        integer, intent(in) :: iteration
        integer             :: wave, iso
        real(KIND = dp)     :: hpsi(nx,ny,nz,4)
        
        do wave=1,nwt
            if(wave .lt. nwn) then
                iso = -1
            else
                iso = +1
            endif

            hpsi = sphamil( hfpsi(:,:,:,:,wave)     ,                          &
            &              hfdpsi(:,:,:,:,:,wave)   ,                          &
            &              hfddpsi(:,:,:,:,:,:,wave),                          &
            &              sx(:,wave), sy(:,wave), sz(:,wave),iso)

            spenergies(wave) = sum(hfpsi(:,:,:,:,wave) * hpsi(:,:,:,:)) * dv
            
            hfpsi(:,:,:,:,wave) = hfpsi(:,:,:,:,wave) - (dt/hbar) * hpsi
        enddo

        call GramSchmidt

    end subroutine 
end module evolution
