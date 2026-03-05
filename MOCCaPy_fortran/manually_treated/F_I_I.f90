module F_I_I

   implicit none

   contains

function pow(f,alpha,eps) result(pf)
  !---------------------------------------------------------------------------
  ! Safely take powers of a REAL density f, avoiding the raising of negative
  ! numbers to powers that are 0 or negative. This is achieved by adding a
  ! small (positive value) to the density.
  !---------------------------------------------------------------------------
  implicit none

  real*8, intent(in)  :: f(:), alpha, eps
  real*8              :: pf(size(f))

  if(alpha .lt. 0) then
     pf = (f + eps)**(alpha)
  else
     pf = (f)**(alpha)
  endif
end function pow

subroutine calc_F_I_I(D_I_I, Lap_D_I_I, D_Nm_Nm, Der_C_I_NS, coupl_constant, sigma, sigma_b, eps, ngrid, F_I_I)
  !--------------------------------------------------------------------------------------------
  ! Calculate the mean-field potential F_I_I from the densities
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Input:
  !    D_I_I         : real*8 array, dimension (ngrid, 2)
  !                    the mean-field density D_I_I at every mesh point
  !                    the final index indicating ISOSPIN, i.e. 1-> t=0 & 2-> t=1.
  !
  !    Lap_D_I_I     : real*8 array, dimension (ngrid,2)
  !                    the laplacian of the mean-field density D_I_I at every mesh point
  !                    the final index indicating ISOSPIN, i.e. 1-> t=0 & 2-> t=1.
  !
  !    D_Nm_Nm       : real*8 array, dimension(ngrid, 2)
  !                    the mean-field density D_Nm_Nm at every mesh point
  !                    the final index indicating ISOSPIN, i.e. 1-> t=0 & 2-> t=1.
  !
  !    Der_C_I_NS    : real*8 array, dimension (ngrid,3,3,3,2)
  !                    the derivative of the mean-field density Der_C_I_NS at every mesh point
  !                    1st index -> space
  !                    2nd index -> the external derivative
  !                    3rd index -> the \nabla internal to C_I_NS
  !                    4th index -> the sigma matrix
  !                    5th index -> isospin (t=0/1)
  !
  !    coupl_constant: real*8 array of length 14
  !                    listing all the coupling constants in the ordering of the BXL.func file
  !
  !    sigma, sigma_b: real*8
  !                    powers of the density dependence
  !
  !    eps           : real*8
  !                    numerical safeguard value to use when taking powers of the density
  !
  !    ngrid         : integer
  !                    total number of grid points
  !
  ! Output:
  !    F_I_I        : real*8 array dimension (ngrid,2)
  !                   the calculated mean-field potential F_I_I at every grid point
  !                   with the final index indicating the isospin
  !--------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: ngrid
  real*8, intent(in)  :: D_I_I(ngrid,2), Lap_D_I_I(ngrid,2), D_Nm_Nm(ngrid,2), Der_C_I_NS(ngrid,3,3,3,2)
  real*8, intent(in)  :: coupl_constant(14)
  real*8, intent(in)  :: sigma, sigma_b, eps
  real*8, intent(out) :: F_I_I(ngrid,2)

  F_I_I = 0.0
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Isospin = 0
  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant( 1)  * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant( 1)  * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant(3) * (sigma)  * pow(D_I_I(:,1), (sigma - 1), eps) * D_I_I(:,1) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant( 3)  * pow(D_I_I(:,1), sigma, eps) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant( 3)  * pow(D_I_I(:,1), sigma, eps) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant(4) * (sigma)  * pow(D_I_I(:,1), (sigma - 1), eps) * D_I_I(:,2) * D_I_I(:,2)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant(5) * (sigma_b)  * pow(D_I_I(:,1), (sigma_b - 1), eps) * D_I_I(:,1) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant( 5)  * pow(D_I_I(:,1), sigma_b, eps) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant( 5)  * pow(D_I_I(:,1), sigma_b, eps) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant(6) * (sigma_b)  * pow(D_I_I(:,1), (sigma_b - 1), eps) * D_I_I(:,2) * D_I_I(:,2)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant( 7)  * D_Nm_Nm(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant( 9)  * Lap_D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant( 9)  * Lap_D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant(13)  * Der_C_I_NS(:,1,2,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & - coupl_constant(13)  * Der_C_I_NS(:,1,3,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant(13)  * Der_C_I_NS(:,2,3,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & - coupl_constant(13)  * Der_C_I_NS(:,2,1,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & + coupl_constant(13)  * Der_C_I_NS(:,3,1,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
       & - coupl_constant(13)  * Der_C_I_NS(:,3,2,1,1)

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Isospin = 1
  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant( 2)  * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant( 2)  * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant( 4)  * pow(D_I_I(:,1), sigma, eps) * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant( 4)  * pow(D_I_I(:,1), sigma, eps) * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant( 6)  * pow(D_I_I(:,1), sigma_b, eps) * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant( 6)  * pow(D_I_I(:,1), sigma_b, eps) * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant( 8)  * D_Nm_Nm(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant(10)  * Lap_D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant(10)  * Lap_D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant(14)  * Der_C_I_NS(:,1,2,3,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & - coupl_constant(14)  * Der_C_I_NS(:,1,3,2,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant(14)  * Der_C_I_NS(:,2,3,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & - coupl_constant(14)  * Der_C_I_NS(:,2,1,3,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & + coupl_constant(14)  * Der_C_I_NS(:,3,1,2,2)

  F_I_I(:,2) = F_I_I(:,2)  &
       & - coupl_constant(14)  * Der_C_I_NS(:,3,2,1,2)

end subroutine calc_F_I_I

end module F_I_I
