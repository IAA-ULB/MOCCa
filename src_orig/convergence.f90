!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
!    Copyright (C) 2026 W. Ryssens and M. Bender
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Affero General Public License as published
!    by the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!===============================================================================
module convergence

 use geninfo

 implicit none

 real(KIND=dp), allocatable :: con_rates(:)   

contains

  subroutine Converged(C, iter)
    !---------------------------------------------------------------------------
    ! Checks if the code has converged using the following convergence criteria.
    !
    ! Input :
    !     iter :  Iteration count
    !
    ! Output : 
    !       C  :  True if the convergence criteria are all satisfied.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !   Keyword         Default    Quantity
    !  ------------    ---------  -------------
    !   energy_prec     1d-9     abs((E^(i) - E^(i-1))/E^(i))     < energy_prec 
    !                            
    !   moment_prec     1d-3     abs((Q2m^(i) - Q2m^(i))/Q2m^(i)) < moment_prec
    !                                         if abs(beta_2m^(i)) > 0.01 
    !                            
    !
    !   disp_prec       1d-5     abs(sum_i v^2_i <psi|h^2|psi> - epsilon^2)
    !                                     < disp_prec
    !
    !   gradient_prec   1d+0     |s.p. gradient|  <    gradient_prec  
    !
    !   fermi_prec      1d-3     abs(lambda^(i) - lambda^(i-1)) < fermi_prec
    !                                    for both nucleon species
    !
    !   angmom_prec     1d-3     abs(<J_mu>^(i) - <J_mu>^(i-1)) < angmom_prec
    !                                    for all cartesian directions
    !
    !   min_iter_conv   -1       A minimum number of iterations to perform
    !
    ! Additional notes:
    !  *   The convergence criterion on the energy is checked for the past 
    !      five (5) iterations, not only for the last one.
    !  *   The angular momentum convergence criterion is trivially satisfied
    !      when the angular momentum values are restricted by symmetry.
    !  *   The multipole moment convergence used to include all multipole
    !      moments, but now only looks at the quadrupole moments.
    !  *   The dispersion condition is trivially satisfied when using the 
    !      HFB gradient solver, as it is simply set to zero in that case.
    !  *   Similarly, the gradient condition is trivially satisfied when using
    !      the HFB solver, as it is simply set to zero in that case.
    !      Note that this criterion is not particularly tight by default, as
    !      I suspect the definition of |s.p. gradient| might evolve.
    !
    ! Things that could be thought about
    ! - - - - - - - - - - - - - - - - - - 
    !  -> Convergence criteria for constraints, both multipole and cranking, 
    !     on values as well as multipliers.
    !---------------------------------------------------------------------------
    use Moments
    use functional
    use evolution
    use cranking

    logical, intent(out)  :: C

    integer             :: i
    integer, intent(in) :: iter
    real(KIND=dp)       :: dE(5), dQ

    type(Moment), pointer  :: Current 

    C = .true.

    ! Checking the evolution of the energy
    do i =1,4
            dE(i) = abs(Ehistory(i) - Ehistory(i+1))/abs(totalE)
    enddo
    dE(5) = abs(TotalE - Ehistory(1))/abs(totalE)
    
    if(.not. all(dE .lt. energy_prec)) then
      C = .false.
    endif

    ! Check all of the quadrupole moments that are large enough
    Current => Root

    do while(associated(Current%next)) 
        Current => Current%next
        if(abs(Current%Beta(3)).gt.0.01 .and. Current%l .eq. 2) then
            dQ = abs(sum(Current%history)-sum(Current%value))
            dQ = dQ/abs(sum(Current%value))
            if(dQ > moment_prec) C = .false.
        endif
    enddo   
    
    ! Checking the weighted dispersion
    if(d2H .gt. disp_prec) C = .false.

    ! Checking the norm of the s.p. gradient
    if(gradientnorm .gt. gradient_prec) C = .false.
    
    ! Check the Fermi energy
    if(any(abs(Fermienergy - FermiHistory).gt.fermi_prec)) C = .false.
        
    ! Check the angular momentum
    if(.not. crank_smooth) then
      ! angular momenta calculated by integration of spin and current densities 
      do i=1,3
        if(abs(TotalAngMom_dens(i)-AngMomOld_dens(i)).gt.angmom_prec ) C = .false.
      enddo
    else
      ! angular momenta calculated by summation of spwf properties
      do i=1,3
        if(abs(TotalAngMom(i) - AngMomOld(i)).gt.angmom_prec ) C = .false.
      enddo    
    endif  
    ! Check if we have performed at least a minimum of iterations
    if(iter.le. min_iter_conv) C = .false.
        
  end subroutine Converged

end module convergence
