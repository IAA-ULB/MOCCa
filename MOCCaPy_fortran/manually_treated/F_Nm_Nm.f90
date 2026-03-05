subroutine calc_F_Nm_Nm(D_I_I, coupl_constant, ngrid, F_Nm_Nm)
  !--------------------------------------------------------------------------------------------
  ! Calculate the mean-field potential F_I_I from the densities
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Input:
  !    D_I_I         : real*8 array, dimension (ngrid, 2)
  !                    the mean-field density D_I_I at every mesh point
  !                    the final index indicating ISOSPIN, i.e. 1-> t=0 & 2-> t=1.

  !
  !    coupl_constant: real*8 array of length 14
  !                    listing all the coupling constants in the ordering of the BXL.func file
  !    ngrid         : integer
  !                    total number of grid points
  !
  ! Output:
  !    F_Nm_Nm       : real*8 array dimension (ngrid,2)
  !                    the calculated mean-field potential F_Nm_Nm at every grid point
  !                    with the final index indicating the isospin
  !--------------------------------------------------------------------------------------------

  integer, intent(in) :: ngrid
  real*8, intent(in)  :: D_I_I(ngrid,2), coupl_constant(14)
  real*8, intent(out) :: F_Nm_Nm(ngrid,2)

  !------------------------------------------------------------------------
  ! Calculation of F_Nm_Nm
  F_Nm_Nm = 0.0
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Isospin = 0
  F_Nm_Nm(:,1) = F_Nm_Nm(:,1)  &
       & + coupl_constant( 7)  * D_I_I(:,1)

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Isospin = 1
  F_Nm_Nm(:,2) = F_Nm_Nm(:,2)  &
       & + coupl_constant( 8)  * D_I_I(:,1)

end subroutine calc_F_Nm_Nm
