subroutine calc_DP_I_I(coupl_constant, sigma, sigma_b, DENLIST, ngrid, DP_I_I)
    !--------------------------------------------------------------------------------------------
    ! Calculate the mean-field potential DP_I_I from the densities
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !    coupl_constant: an real*8 array of all the coupling constants
    !    ngrid         : total number of grid points
    !
    ! Output:
    !    DP_I_I    : the calculated mean-field potential DP_I_I at every grid point
    !               dimension ( ngrid, 2 )
    !               with the final index indicating the isospin (1 = neutrons, 2 = protons)
    !--------------------------------------------------------------------------------------------

    implicit none

DECLARATION
    real*8, intent(in)  :: coupl_constant(:)
    integer, intent(in) :: ngrid
    real*8, intent(out) :: DP_I_I(ngrid,2)

    real*8 :: temp


    !---------------------------------------------------------------------------
    ! Calculation of F_I_I 
        F_I_I = 0.0 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Isospin = 0 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant( 1)  * D_I_I(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant( 1)  * D_I_I(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant(3) * (sigma)  * pow(D_I_I(:,3), (sigma - 1)) * D_I_I(:,3) * D_I_I(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant( 3)  * pow(D_I_I(:,3), sigma) * D_I_I(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant( 3)  * pow(D_I_I(:,3), sigma) * D_I_I(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant(4) * (sigma)  * pow(D_I_I(:,3), (sigma - 1)) * D_I_I(:,4) * D_I_I(:,4) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant(5) * (sigma_b)  * pow(D_I_I(:,3), (sigma_b - 1)) * D_I_I(:,3) * D_I_I(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant( 5)  * pow(D_I_I(:,3), sigma_b) * D_I_I(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant( 5)  * pow(D_I_I(:,3), sigma_b) * D_I_I(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant(6) * (sigma_b)  * pow(D_I_I(:,3), (sigma_b - 1)) * D_I_I(:,4) * D_I_I(:,4) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant( 7)  * D_Nm_Nm(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant( 9)  * Lap_D_I_I(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant( 9)  * Lap_D_I_I(:,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant(13)  * Der_C_I_NS(:,1,2,3,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & - coupl_constant(13)  * Der_C_I_NS(:,1,3,2,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant(13)  * Der_C_I_NS(:,2,3,1,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & - coupl_constant(13)  * Der_C_I_NS(:,2,1,3,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & + coupl_constant(13)  * Der_C_I_NS(:,3,1,2,3) 
 
        F_I_I(:,3) = F_I_I(:,3)  & 
        & - coupl_constant(13)  * Der_C_I_NS(:,3,2,1,3) 
 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Isospin = 1 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant( 2)  * D_I_I(:,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant( 2)  * D_I_I(:,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant( 4)  * pow(D_I_I(:,3), sigma) * D_I_I(:,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant( 4)  * pow(D_I_I(:,3), sigma) * D_I_I(:,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant( 6)  * pow(D_I_I(:,3), sigma_b) * D_I_I(:,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant( 6)  * pow(D_I_I(:,3), sigma_b) * D_I_I(:,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant( 8)  * D_Nm_Nm(:,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant(10)  * Lap_D_I_I(:,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant(10)  * Lap_D_I_I(:,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant(14)  * Der_C_I_NS(:,1,2,3,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & - coupl_constant(14)  * Der_C_I_NS(:,1,3,2,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant(14)  * Der_C_I_NS(:,2,3,1,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & - coupl_constant(14)  * Der_C_I_NS(:,2,1,3,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & + coupl_constant(14)  * Der_C_I_NS(:,3,1,2,4) 
 
        F_I_I(:,4) = F_I_I(:,4)  & 
        & - coupl_constant(14)  * Der_C_I_NS(:,3,2,1,4) 
 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Recombination 
        F_I_I(:,1) = F_I_I(:,3) +F_I_I(:,4)
        F_I_I(:,2) = F_I_I(:,3) -F_I_I(:,4) 
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Calculation of F_Nm_Nm 
        F_Nm_Nm = 0.0 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Isospin = 0 
        F_Nm_Nm(:,3) = F_Nm_Nm(:,3)  & 
        & + coupl_constant( 7)  * D_I_I(:,3) 
 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Isospin = 1 
        F_Nm_Nm(:,4) = F_Nm_Nm(:,4)  & 
        & + coupl_constant( 8)  * D_I_I(:,4) 
 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Recombination 
        F_Nm_Nm(:,1) = F_Nm_Nm(:,3) +F_Nm_Nm(:,4)
        F_Nm_Nm(:,2) = F_Nm_Nm(:,3) -F_Nm_Nm(:,4) 
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Calculation of G_I_NS 
        G_I_NS = 0.0 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Isospin = 0 
        G_I_NS(:,1,1,3) = G_I_NS(:,1,1,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,1,1,3) 
 
        G_I_NS(:,1,2,3) = G_I_NS(:,1,2,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,1,2,3) 
 
        G_I_NS(:,1,3,3) = G_I_NS(:,1,3,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,1,3,3) 
 
        G_I_NS(:,2,1,3) = G_I_NS(:,2,1,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,2,1,3) 
 
        G_I_NS(:,2,2,3) = G_I_NS(:,2,2,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,2,2,3) 
 
        G_I_NS(:,2,3,3) = G_I_NS(:,2,3,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,2,3,3) 
 
        G_I_NS(:,3,1,3) = G_I_NS(:,3,1,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,3,1,3) 
 
        G_I_NS(:,3,2,3) = G_I_NS(:,3,2,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,3,2,3) 
 
        G_I_NS(:,3,3,3) = G_I_NS(:,3,3,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,3,3,3) 
 
        G_I_NS(:,1,1,3) = G_I_NS(:,1,1,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,1,1,3) 
 
        G_I_NS(:,1,2,3) = G_I_NS(:,1,2,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,1,2,3) 
 
        G_I_NS(:,1,3,3) = G_I_NS(:,1,3,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,1,3,3) 
 
        G_I_NS(:,2,1,3) = G_I_NS(:,2,1,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,2,1,3) 
 
        G_I_NS(:,2,2,3) = G_I_NS(:,2,2,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,2,2,3) 
 
        G_I_NS(:,2,3,3) = G_I_NS(:,2,3,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,2,3,3) 
 
        G_I_NS(:,3,1,3) = G_I_NS(:,3,1,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,3,1,3) 
 
        G_I_NS(:,3,2,3) = G_I_NS(:,3,2,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,3,2,3) 
 
        G_I_NS(:,3,3,3) = G_I_NS(:,3,3,3)  & 
        & + coupl_constant(11)  * C_I_NS(:,3,3,3) 
 
        G_I_NS(:,2,3,3) = G_I_NS(:,2,3,3)  & 
        & - coupl_constant(13)  * Der_D_I_I(:,1,3) 
 
        G_I_NS(:,3,2,3) = G_I_NS(:,3,2,3)  & 
        & + coupl_constant(13)  * Der_D_I_I(:,1,3) 
 
        G_I_NS(:,3,1,3) = G_I_NS(:,3,1,3)  & 
        & - coupl_constant(13)  * Der_D_I_I(:,2,3) 
 
        G_I_NS(:,1,3,3) = G_I_NS(:,1,3,3)  & 
        & + coupl_constant(13)  * Der_D_I_I(:,2,3) 
 
        G_I_NS(:,1,2,3) = G_I_NS(:,1,2,3)  & 
        & - coupl_constant(13)  * Der_D_I_I(:,3,3) 
 
        G_I_NS(:,2,1,3) = G_I_NS(:,2,1,3)  & 
        & + coupl_constant(13)  * Der_D_I_I(:,3,3) 
 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Isospin = 1 
        G_I_NS(:,1,1,4) = G_I_NS(:,1,1,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,1,1,4) 
 
        G_I_NS(:,1,2,4) = G_I_NS(:,1,2,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,1,2,4) 
 
        G_I_NS(:,1,3,4) = G_I_NS(:,1,3,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,1,3,4) 
 
        G_I_NS(:,2,1,4) = G_I_NS(:,2,1,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,2,1,4) 
 
        G_I_NS(:,2,2,4) = G_I_NS(:,2,2,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,2,2,4) 
 
        G_I_NS(:,2,3,4) = G_I_NS(:,2,3,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,2,3,4) 
 
        G_I_NS(:,3,1,4) = G_I_NS(:,3,1,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,3,1,4) 
 
        G_I_NS(:,3,2,4) = G_I_NS(:,3,2,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,3,2,4) 
 
        G_I_NS(:,3,3,4) = G_I_NS(:,3,3,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,3,3,4) 
 
        G_I_NS(:,1,1,4) = G_I_NS(:,1,1,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,1,1,4) 
 
        G_I_NS(:,1,2,4) = G_I_NS(:,1,2,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,1,2,4) 
 
        G_I_NS(:,1,3,4) = G_I_NS(:,1,3,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,1,3,4) 
 
        G_I_NS(:,2,1,4) = G_I_NS(:,2,1,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,2,1,4) 
 
        G_I_NS(:,2,2,4) = G_I_NS(:,2,2,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,2,2,4) 
 
        G_I_NS(:,2,3,4) = G_I_NS(:,2,3,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,2,3,4) 
 
        G_I_NS(:,3,1,4) = G_I_NS(:,3,1,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,3,1,4) 
 
        G_I_NS(:,3,2,4) = G_I_NS(:,3,2,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,3,2,4) 
 
        G_I_NS(:,3,3,4) = G_I_NS(:,3,3,4)  & 
        & + coupl_constant(12)  * C_I_NS(:,3,3,4) 
 
        G_I_NS(:,2,3,4) = G_I_NS(:,2,3,4)  & 
        & - coupl_constant(14)  * Der_D_I_I(:,1,4) 
 
        G_I_NS(:,3,2,4) = G_I_NS(:,3,2,4)  & 
        & + coupl_constant(14)  * Der_D_I_I(:,1,4) 
 
        G_I_NS(:,3,1,4) = G_I_NS(:,3,1,4)  & 
        & - coupl_constant(14)  * Der_D_I_I(:,2,4) 
 
        G_I_NS(:,1,3,4) = G_I_NS(:,1,3,4)  & 
        & + coupl_constant(14)  * Der_D_I_I(:,2,4) 
 
        G_I_NS(:,1,2,4) = G_I_NS(:,1,2,4)  & 
        & - coupl_constant(14)  * Der_D_I_I(:,3,4) 
 
        G_I_NS(:,2,1,4) = G_I_NS(:,2,1,4)  & 
        & + coupl_constant(14)  * Der_D_I_I(:,3,4) 
 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Recombination 
        G_I_NS(:,:,:,1) = G_I_NS(:,:,:,3) +G_I_NS(:,:,:,4)
        G_I_NS(:,:,:,2) = G_I_NS(:,:,:,3) -G_I_NS(:,:,:,4) 
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Calculation of FP_I_I 
        FP_I_I = 0.0 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Isospin = n 
        FP_I_I(:,1) = FP_I_I(:,1)  & 
        & + coupl_constant(15)  * DP_I_I(:,1) 
 
        FP_I_I(:,1) = FP_I_I(:,1)  & 
        & + coupl_constant(15)  * DP_I_I(:,1) 
 
        FP_I_I(:,1) = FP_I_I(:,1)  & 
        & + coupl_constant(17)  * pow(D_I_I(:,3), sigmap) * DP_I_I(:,1) 
 
        FP_I_I(:,1) = FP_I_I(:,1)  & 
        & + coupl_constant(17)  * pow(D_I_I(:,3), sigmap) * DP_I_I(:,1) 
 
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Isospin = p 
        FP_I_I(:,2) = FP_I_I(:,2)  & 
        & + coupl_constant(16)  * DP_I_I(:,2) 
 
        FP_I_I(:,2) = FP_I_I(:,2)  & 
        & + coupl_constant(16)  * DP_I_I(:,2) 
 
        FP_I_I(:,2) = FP_I_I(:,2)  & 
        & + coupl_constant(18)  * pow(D_I_I(:,3), sigmap) * DP_I_I(:,2) 
 
        FP_I_I(:,2) = FP_I_I(:,2)  & 
        & + coupl_constant(18)  * pow(D_I_I(:,3), sigmap) * DP_I_I(:,2) 
 
    !---------------------------------------------------------------------------



end subroutine calc_DP_I_I

pure function pow(f,alpha) result(pf)
    !---------------------------------------------------------------------------
    ! Safely take powers of a REAL density f, avoiding the raising of negative
    ! numbers to powers that are 0 or negative. This is achieved by adding a
    ! small (positive value) to the density.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: f(mv), alpha
    real(KIND=dp)             :: pf(mv)

    if(alpha .lt. 0) then
      pf = (f + eps)**(alpha)
    else
      pf = (f)**(alpha)
    endif
end function pow
