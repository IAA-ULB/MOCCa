module convergence

 !==============================================================================
 !_________ _______  _       _________ _______  _                 _______ 
 !\__   __/(  ___  )( (    /|\__   __/(  ___  )( \      |\     /|(  ____ \
 !   ) (   | (   ) ||  \  ( |   ) (   | (   ) || (      | )   ( || (    \/
 !   | |   | (___) ||   \ | |   | |   | (___) || |      | |   | || (_____ 
 !   | |   |  ___  || (\ \) |   | |   |  ___  || |      | |   | |(_____  )
 !   | |   | (   ) || | \   |   | |   | (   ) || |      | |   | |      ) |
 !   | |   | )   ( || )  \  |   | |   | )   ( || (____/\| (___) |/\____) |
 !   )_(   |/     \||/    )_)   )_(   |/     \|(_______/(_______)\_______)
 !                                                                       
 !  Copyright W. Ryssens & M. Bender
 !
 !==============================================================================
 use Geninfo

 implicit none


 real(KIND=dp), allocatable :: con_rates(:)   

contains

!  subroutine monitor_convergence(iter)
!    !---------------------------------------------------------------------------
!    ! We calculate estimates for various convergence rates. 
!    !
!    !
!    !
!    !---------------------------------------------------------------------------

!    use functional
!    use moments

!    integer, intent(in)   :: iter
!    integer               :: N, i
!    real(KIND=dp)         :: rate, ds, dsold, D
!    type(Moment), pointer :: mom 
!    
!    !if(allocated(con_rates)) deallocate(con_rates)
!    
!    !---------------------------------------------------------------------------
!    ! Count how many rates we need
!    !N = 3     ! Energy + complete Routhian + (E_func - E_spwf)

!    !mom => Root
!    !do while(associated(mom%next))
!    !    mom => mom%next
!    !    if(mom%constrainttype.ne.0) N = N+1
!    !enddo

!    !---------------------------------------------------------------------------
!    if(.not. allocated(con_rates)) then
!        allocate(con_rates(4)) ; con_rates = 0.0d0
!    endif

!    ! Rate of total energy
!    D = log(abs(totalE-EHistory(1))/abs(EHistory(1)-EHistory(2)))     
!    D = D/log(abs(Ehistory(1)-EHistory(2))/abs(EHistory(2)-EHistory(3)))     
!    D = D - con_rates(1)
!    con_rates(1) =  con_rates(1)  + 0.1 * D 
!    ! Rate of the Routhian
!    D = log(abs(Routhian-RHistory(1))/abs(RHistory(1) - RHistory(2)))
!    D = D/log(abs(RHistory(1)-RHistory(2))/abs(RHistory(2) - RHistory(3)))
!    D = D  - con_rates(2)
!    con_rates(2) = con_rates(2)   + 0.1 * D
!    ! Rate of E_spwf - E_fun
!    ds    =   Spwfenergy     - Spwfhistory(1) - totalE      + Ehistory(1) 
!    dsold =   Spwfhistory(1) - Spwfhistory(2) - Ehistory(1) + Ehistory(2)
!    D = log(abs(ds)/abs(dsold)) 

!    ds    =   Spwfhistory(1) - Spwfhistory(2) - Ehistory(1) + Ehistory(2)
!    dsold =   Spwfhistory(2) - Spwfhistory(3) - Ehistory(2) + Ehistory(3)
!    D = D / log(abs(ds)/abs(dsold)) 
!    D = D - con_rates(3)
!    con_rates(3) = con_rates(3) + 0.1 * D
!  
!    ! Rate of the change in density
!    D = sqrt(sum((D_I_I - D_I_I_hist(:,:,1))**2)*dv) / sqrt(sum((D_I_I_hist(:,:,1) - D_I_I_hist(:,:,2))**2)*dv)
!    D = log(D) / &
!    &   log(sqrt(sum((D_I_I_hist(:,:,1) - D_I_I_hist(:,:,2))**2)*dv) / &
!    &       sqrt(sum((D_I_I_hist(:,:,3) - D_I_I_hist(:,:,2))**2)*dv))
!    D = D - con_rates(4)
!    con_rates(4) = con_rates(4) + 0.1 * D



!    do i=1, size( con_rates)
!        if(abs(con_rates(i)) .lt. 1d-16) con_rates(i) = 0.0d0
!    enddo
!  end subroutine monitor_convergence

  subroutine Converged(C) 
    !---------------------------------------------------------------------------
    ! Checks if the code has converged using the following convergence 
    ! criteria, all of which need to be verified across 5 iterations.
    !
    !   Keyword         Default    Quantity
    !  ------------    ---------  -------------
    !   energy_prec     1d-9     abs((E^(i) - E^(i-1))/E^(i))     < energy_prec 
    !
    !   moment_prec     1d-3     abs((Qlm^(i) - Qlm^(i))/Qlm^(i)) < moment_prec
    !                                if Qlm^(i) is large enough
    !
    !   disp_prec       1d-5     abs(sum_i v^2_i <psi|h^2|psi> - epsilon^2)
    !                                     < disp_prec
    !
    ! 
    !---------------------------------------------------------------------------
    use Moments
    use functional
    use evolution

    logical       :: C
    integer       :: i
    real(KIND=dp) :: dE(5), dQ

    type(Moment), pointer  :: Current 

    C = .true.

    !---------------------------------------------------------------------------
    ! Checking the evolution of the energy
    do i =1,4
            dE(i) = abs(Ehistory(i) - Ehistory(i+1))/abs(totalE)
    enddo
    dE(5) = abs(TotalE - Ehistory(1))/abs(totalE)
    
    if(.not. all(dE .lt. energy_prec)) then
     C = .false.
    endif

    !---------------------------------------------------------------------------
    ! Checking the weighted dispersion
    if(d2H .gt. disp_prec) C = .false.
    !---------------------------------------------------------------------------
    ! Check all of the multipole moments that are large enough
    Current => Root

    do while(associated(Current%next)) 
        Current => Current%next
        if(Current%Beta(3).gt.0.05) then
            dQ = abs(sum(Current%history)-sum(Current%value))
            dQ = dQ/abs(sum(Current%value))
            if(dQ > moment_prec) C = .false.
        endif
    enddo   
  end subroutine Converged

end module convergence


