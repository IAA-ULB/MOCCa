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

subroutine calc_F_I_I(D_I_I, Der_D_I_I, Der_Der_D_I_I, &
    &                        Lap_D_I_I, D_Nm_Nm, Der_C_I_NS, &
    &                        coupl_constant, sigma, beta, gamma, eps, ngrid, F_I_I)
  !--------------------------------------------------------------------------------------------
  ! Calculate the mean-field potential F_I_I from the densities
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Input:
  !    D_I_I         : real*8 array, dimension (ngrid, 2)
  !                    the mean-field density D_I_I at every mesh point
  !                    the final index indicating ISOSPIN, i.e. 1-> t=0 & 2-> t=1.
  !
  !    Der_D_I_I     : real*8 array, dimension (ngrid, 3, 2)
  !                    the gradient of the density D_I_I at every mesh point
  !                    the final index indicating ISOSPIN, i.e. 1-> t=0 & 2-> t=1.
  !
  !    Der_Der_D_I_I  : real*8 array, dimension (ngrid, 6, 2)
  !                    the complete Hessian of the density D_I_I at every mesh point
  !                    in packed storage, i.e. the second index is arranged in lexicographical ordering
  !                     1 : xx
  !                     2 : xy
  !                     3 : xz
  !                     4 : yy
  !                     5 : yz
  !                     6 : zz
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
  !    coupl_constant: real*8 array of length 21
  !                    listing all the coupling constants in the ordering of the BXL.func file
  !
  !    sigma, beta : real*8
  !    gamma         powers of the density dependence
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
  real*8, intent(in)  :: D_I_I(ngrid,2), der_D_I_I(ngrid,3,2), der_der_D_I_I(ngrid,6,2)
  real*8, intent(in)  :: Lap_D_I_I(ngrid,2), D_Nm_Nm(ngrid,2), Der_C_I_NS(ngrid,3,3,3,2)
  real*8, intent(in)  :: coupl_constant(22)
  real*8, intent(in)  :: sigma, beta, gamma, eps
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
  & + coupl_constant( 5)  * D_Nm_Nm(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant( 7)  * Lap_D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant( 7)  * Lap_D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(11) * (beta)  * pow(D_I_I(:,1), (beta - 1), eps) * D_I_I(:,1) * D_Nm_Nm(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(11)  * pow(D_I_I(:,1), beta, eps) * D_Nm_Nm(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(12) * (beta)  * pow(D_I_I(:,1), (beta - 1), eps) * D_I_I(:,2) * D_Nm_Nm(:,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(13) * (gamma)  * pow(D_I_I(:,1), (gamma - 1), eps) * D_I_I(:,1) * D_Nm_Nm(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(13)  * pow(D_I_I(:,1), gamma, eps) * D_Nm_Nm(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(14) * (gamma)  * pow(D_I_I(:,1), (gamma - 1), eps) * D_I_I(:,2) * D_Nm_Nm(:,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(15) * (beta)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(15) * (beta)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(15) * (beta)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta) * coupl_constant(15)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta) * coupl_constant(15)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta) * coupl_constant(15)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(15)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(15)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,4,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(15)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,6,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta) * coupl_constant(15)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta) * coupl_constant(15)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta) * coupl_constant(15)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(15)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(15)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,4,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(15)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,6,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(16) * (beta)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,1,2) * Der_D_I_I(:,1,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(16) * (beta)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,2,2) * Der_D_I_I(:,2,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(16) * (beta)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,3,2) * Der_D_I_I(:,3,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(17) * (gamma)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(17) * (gamma)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(17) * (gamma)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (gamma) * coupl_constant(17)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (gamma) * coupl_constant(17)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (gamma) * coupl_constant(17)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(17)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(17)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,4,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(17)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,6,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (gamma) * coupl_constant(17)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (gamma) * coupl_constant(17)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (gamma) * coupl_constant(17)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(17)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(17)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,4,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(17)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,6,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(18) * (gamma)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,1,2) * Der_D_I_I(:,1,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(18) * (gamma)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,2,2) * Der_D_I_I(:,2,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(18) * (gamma)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,3,2) * Der_D_I_I(:,3,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(19) * (beta-1)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,1,1) * D_I_I(:,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(19) * (beta-1)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,2,1) * D_I_I(:,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(19) * (beta-1)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,3,1) * D_I_I(:,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta-1) * coupl_constant(19)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,1,1) * D_I_I(:,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta-1) * coupl_constant(19)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,2,1) * D_I_I(:,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta-1) * coupl_constant(19)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,3,1) * D_I_I(:,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * D_I_I(:,1) * Der_Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * D_I_I(:,1) * Der_Der_D_I_I(:,4,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * D_I_I(:,1) * Der_Der_D_I_I(:,6,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta-1) * coupl_constant(19)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta-1) * coupl_constant(19)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta-1) * coupl_constant(19)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_Der_D_I_I(:,1,1) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_Der_D_I_I(:,4,1) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_Der_D_I_I(:,6,1) * D_I_I(:,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(19)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(20) * (beta-1)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,1,1) * D_I_I(:,2) * Der_D_I_I(:,1,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(20) * (beta-1)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,2,1) * D_I_I(:,2) * Der_D_I_I(:,2,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(20) * (beta-1)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,3,1) * D_I_I(:,2) * Der_D_I_I(:,3,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta-1) * coupl_constant(20)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,1,1) * D_I_I(:,2) * Der_D_I_I(:,1,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta-1) * coupl_constant(20)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,2,1) * D_I_I(:,2) * Der_D_I_I(:,2,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - (beta-1) * coupl_constant(20)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,3,1) * D_I_I(:,2) * Der_D_I_I(:,3,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,1,2) * Der_D_I_I(:,1,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,2,2) * Der_D_I_I(:,2,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,3,2) * Der_D_I_I(:,3,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * D_I_I(:,2) * Der_Der_D_I_I(:,1,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * D_I_I(:,2) * Der_Der_D_I_I(:,4,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * D_I_I(:,2) * Der_Der_D_I_I(:,6,2)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(21)  * Der_C_I_NS(:,1,2,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(21)  * Der_C_I_NS(:,1,3,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(21)  * Der_C_I_NS(:,2,3,1,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(21)  * Der_C_I_NS(:,2,1,3,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & + coupl_constant(21)  * Der_C_I_NS(:,3,1,2,1)

  F_I_I(:,1) = F_I_I(:,1)  &
  & - coupl_constant(21)  * Der_C_I_NS(:,3,2,1,1)

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
  & + coupl_constant( 6)  * D_Nm_Nm(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & + coupl_constant( 8)  * Lap_D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & + coupl_constant( 8)  * Lap_D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & + coupl_constant(12)  * pow(D_I_I(:,1), beta, eps) * D_Nm_Nm(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & + coupl_constant(14)  * pow(D_I_I(:,1), gamma, eps) * D_Nm_Nm(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (beta) * coupl_constant(16)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (beta) * coupl_constant(16)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (beta) * coupl_constant(16)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(16)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(16)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,4,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(16)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,6,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (beta) * coupl_constant(16)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (beta) * coupl_constant(16)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (beta) * coupl_constant(16)  * pow(D_I_I(:,1), (beta - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(16)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(16)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,4,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(16)  * pow(D_I_I(:,1), beta, eps) * Der_Der_D_I_I(:,6,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (gamma) * coupl_constant(18)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (gamma) * coupl_constant(18)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (gamma) * coupl_constant(18)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(18)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(18)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,4,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(18)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,6,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (gamma) * coupl_constant(18)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (gamma) * coupl_constant(18)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (gamma) * coupl_constant(18)  * pow(D_I_I(:,1), (gamma - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(18)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(18)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,4,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(18)  * pow(D_I_I(:,1), gamma, eps) * Der_Der_D_I_I(:,6,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & + coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & + coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & + coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (beta-1) * coupl_constant(20)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,1) * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (beta-1) * coupl_constant(20)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,1) * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - (beta-1) * coupl_constant(20)  * pow(D_I_I(:,1), (beta-1 - 1), eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,1) * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_Der_D_I_I(:,1,1) * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_Der_D_I_I(:,4,1) * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_Der_D_I_I(:,6,1) * D_I_I(:,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,1,1) * Der_D_I_I(:,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,2,1) * Der_D_I_I(:,2,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(20)  * pow(D_I_I(:,1), beta-1, eps) * Der_D_I_I(:,3,1) * Der_D_I_I(:,3,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & + coupl_constant(22)  * Der_C_I_NS(:,1,2,3,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(22)  * Der_C_I_NS(:,1,3,2,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & + coupl_constant(22)  * Der_C_I_NS(:,2,3,1,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(22)  * Der_C_I_NS(:,2,1,3,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & + coupl_constant(22)  * Der_C_I_NS(:,3,1,2,2)

  F_I_I(:,2) = F_I_I(:,2)  &
  & - coupl_constant(22)  * Der_C_I_NS(:,3,2,1,2)

end subroutine calc_F_I_I

end module F_I_I
