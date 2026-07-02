module F_Nm_Nm

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

subroutine calc_F_Nm_Nm(D_I_I, coupl_constant, beta, gamma, eps, ngrid, F_Nm_Nm)
  !--------------------------------------------------------------------------------------------
  ! Calculate the mean-field potential F_Nm_Nm from the densities
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Input:
  !    D_I_I         : real*8 array, dimension (ngrid, 2)
  !                    the mean-field density D_I_I at every mesh point
  !                    the final index indicating ISOSPIN, i.e. 1-> t=0 & 2-> t=1.
  !
  !
  !    coupl_constant: real*8 array of length 14
  !                    listing all the coupling constants in the ordering of the BXL.func file
  !
  !    beta, gamma   : real*8
  !                    two of the density dependencies
  !
  !    eps           : real*8
  !                    numerical safeguard value to use when taking powers of the density
  !
  !    ngrid         : integer
  !                    total number of grid points
  !
  ! Output:
  !    F_Nm_Nm       : real*8 array dimension (ngrid,2)
  !                    the calculated mean-field potential F_Nm_Nm at every point grid
  !                    with the final index indicating the isospin
  !--------------------------------------------------------------------------------------------

  integer, intent(in) :: ngrid
  real*8, intent(in)  :: D_I_I(ngrid,2), coupl_constant(22)
  real*8, intent(in)  :: eps, beta, gamma
  real*8, intent(out) :: F_Nm_Nm(ngrid,2)

  !------------------------------------------------------------------------
  ! Calculation of F_Nm_Nm
  F_Nm_Nm = 0.0
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Isospin = 0
  F_Nm_Nm(:,1) = F_Nm_Nm(:,1)  &
  & + coupl_constant( 5)  * D_I_I(:,1)

  F_Nm_Nm(:,1) = F_Nm_Nm(:,1)  &
  & + coupl_constant(11)  * pow(D_I_I(:,1), beta, eps) * D_I_I(:,1)

  F_Nm_Nm(:,1) = F_Nm_Nm(:,1)  &
  & + coupl_constant(13)  * pow(D_I_I(:,1), gamma, eps) * D_I_I(:,1)

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Isospin = 1
  F_Nm_Nm(:,2) = F_Nm_Nm(:,2)  &
  & + coupl_constant( 6)  * D_I_I(:,2)

  F_Nm_Nm(:,2) = F_Nm_Nm(:,2)  &
  & + coupl_constant(12)  * pow(D_I_I(:,1), beta, eps) * D_I_I(:,2)

  F_Nm_Nm(:,2) = F_Nm_Nm(:,2)  &
  & + coupl_constant(14)  * pow(D_I_I(:,1), gamma, eps) * D_I_I(:,2)

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Isospin = 0
  F_Nm_Nm(:,1) = F_Nm_Nm(:,1)  &
       & + coupl_constant( 7)  * D_I_I(:,1)

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Isospin = 1
  F_Nm_Nm(:,2) = F_Nm_Nm(:,2)  &
       & + coupl_constant( 8)  * D_I_I(:,1)

end subroutine calc_F_Nm_Nm
end module
