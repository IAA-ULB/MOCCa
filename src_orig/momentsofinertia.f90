module momentsofinertia

  use geninfo
  use densities
  use parameterization

  implicit none

  real(KIND=dp) :: Belyaev(3,3), Rigid(3,3), J2(3,3)
  real(KIND=dp) :: J2_coll(3,3)

  procedure(calcJ2andBelyaev_HF), pointer :: calcJ2andBelyaev 

contains

  subroutine setBelyaevProcedure()
    select case(Pairingtype)
    case(0)
      calcJ2andBelyaev => calcJ2andBelyaev_HF 
    case(1)
      calcJ2andBelyaev => calcJ2andBelyaev_BCS
    case(2)
      calcJ2andBelyaev => calcJ2andBelyaev_HFB
    end select
  end subroutine

  subroutine calcrigid()
    !---------------------------------------------------------------------------
    ! Calculate the rgid-rotor moment of inertia of the density.
    ! 
    !   I^rigid_x = m int d^3r rho * (y^2 + z^2)
    !   I^rigid_y = m int d^3r rho * (x^2 + z^2)
    !   I^rigid_z = m int d^3r rho * (x^2 + y^2)
    !
    ! where m is the mass of the nucleon species under consideration. 
    !
    ! When naively calculated by the code, these numbers have units
    !       MeV * fm^2/c^2
    ! We divide by (hbar c)**2 to obtain the final units of
    !       hbar^2 / MeV
    !---------------------------------------------------------------------------
    real(KIND=dp) :: xs(2), ys(2), zs(2)
    real(KIND=dp), pointer :: rho(:,:,:)
    integer       :: it, i,j,k

    Rigid = 0

    xs = 0 ; ys = 0 ; zs = 0

    do it=1,2
      ! Assigning the storage structure in D_I_I a more readable form
      rho(1:nx, 1:ny, 1:nz) => D_I_I(1:nx*ny*nz,it)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            xs(it) = xs(it) + meshx(i)**2 * rho(i,j,k)
            ys(it) = ys(it) + meshy(j)**2 * rho(i,j,k)
            zs(it) = zs(it) + meshz(k)**2 * rho(i,j,k)
          enddo
        enddo
      enddo
    enddo

    xs = xs * dv ; ys = ys *dv ; zs = zs * dv
      
    Rigid(1,1:2) = nucleonmass * (ys + zs)
    Rigid(2,1:2) = nucleonmass * (xs + zs)
    Rigid(3,1:2) = nucleonmass * (xs + ys)

    Rigid(:,3) = sum(Rigid(:,1:2),2)

    ! Converting to the correct units
    Rigid = Rigid/(hbarclum**2)

  end subroutine calcrigid

  subroutine calcJ2andBelyaev_HF
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------

    integer       :: i,j, b, it, ii, jj, si
    real(KIND=dp) :: ME(3), fi, fj, dfde

    J2 = 0  ;  Belyaev = 0
 
    si = 0  
    do b = 1, Blocks
      do i=1, HFBlocks(b)
        ii = si + i
        it = 1
        if(ii.gt.nwn) it = 2
        do j=1,HFblocks(b)
          jj = si + j
        
          ! |< k | j_x | -l >|^2  
          ME(1)= angmom_xt_real(hfpsi(:,:,ii),hfpsi(:,:,jj),hfdpsi(:,:,:,jj))**2 
          ! |< k | j_y | -l >|^2 
          ME(2)= angmom_yt_imag(hfpsi(:,:,ii),hfpsi(:,:,jj),hfdpsi(:,:,:,jj))**2  
          ! |< k | j_z |  l >|^2 
          ME(3)= angmom_z_real( hfpsi(:,:,ii),hfpsi(:,:,jj),hfdpsi(:,:,:,jj))**2  

          ! The factor two is because of the presence of the four terms in 
          ! ME_{kl}, which are pair-wise equal when J_i = J_j.
          ME = 2 * ME 

          fi = rho_can(ii)/2. ; fj = rho_can(jj)/2.

          J2(:,it) = J2(:,it) + ME * fi*(1-fj)

          if(inversetemp.eq.-1) then
            dfdE = fj - fi
            if(abs(dfdE).gt.0) then
                dfdE = dfdE/(spenergies(ii) - spenergies(jj))
            endif
          else
            if(abs(spenergies(ii) - spenergies(jj)).gt.1d-8) then
              dfdE = (fj - fi)/(spenergies(ii) - spenergies(jj))
            else                
              dfdE = fi**2 * inversetemp                                     &
              &            * exp(inversetemp*(spenergies(ii)-FermiEnergy(it)))
            endif
          endif

          Belyaev(:,it) = Belyaev(:,it) + ME * dfde
        enddo
      enddo
      si = si +   HFBlocks(b)
    enddo
    !  Sum for the total
    J2(:,3) = sum(J2(:,1:2),2) ; Belyaev(:,3) = sum(Belyaev(:,1:2),2)
  end subroutine 

  subroutine calcJ2andBelyaev_BCS
    !---------------------------------------------------------------------------
    ! Calculate the expectation value of J^2_mu in the many-body state, as well
    ! as the Belyaev moment of inertia. 
    ! 
    ! Full expressions for both, in the BCS case at finite temperature 
    ! can be found in
    ! 
    !   Y. Alhassid et al., Phys. Rev. C 72, 064326 (2005).
    ! 
    !
    !  <J_i J_j> = sum_[ kl > 0 ]  ME_{kl} * weight_{kl}
    !    I_ij    = sum_[ kl > 0 ]  ME_{kl} * weight^B_{kl}
    !   
    ! where
    !      ME_{kl} = < k|j_i| l> < l|j_j| k>  + < l|j_i| k> < k|j_j| l>
    !              + < k|j_i|-l> <-l|j_j| k>  + <-k|j_i| l> < l|j_j|-k>
    !    
    ! and
    ! 
    !     weight_{kl}  = [u^2_k u_l^2 + u_k v_k u_l v_l] * f_k (1-f_l)     (a)
    !                  + [v^2_k v_l^2 + u_k v_k u_l v_l] * (1-f_k) f_l     (b)
    !                  + [u_k^2 v_l^2 - u_k v_k u_l v_l] * f_k f_l         (c)
    !                  + [v_k^2 u_l^2 - u_k v_k u_l v_l] * (1-f_k) (1-f_l) (d)
    !
    !
    !    weight^B_{kl} = (u_k u_l + v_k v_l)**2 * (f_l - f_k)/(E_k - E_l)  (Ba)
    !                  + (u_k v_l - v_k u_l)**2 * (1-f_k-f_l)/(E_k + E_l)  (Bb)
    !
    ! with the additional caveat that 
    !
    !   (f_l - f_k)/(E_k - E_l) => - df_k/dE_k as E_k => E_l
    !
    ! - df_k/dE_k = beta exp(beta * E_k)* f_k^2
    !
    ! Due to conserved signature and time-simplexes, the only matrix elements
    ! that enter the calculation of ME_{kl} are of the form
    !
    !  Re < i | J_x T  | j > 
    !  Im < i | J_y T  | j >
    !  Re < i | J_z    | j > 
    !
    ! In addition, we only calculate the diagonal moments of the Inertia tensor.
    !---------------------------------------------------------------------------
    ! Mental note to self.
    ! - - - - - - - - - - - -
    ! This routine has been checked
    !  *) Belyaev & J2 are close to values produced by EV8 for Zr84, 
    !     though not identical due to the difference in derivatives.
    !  *) BCS formulas and HF formulas agree
    !       a) when the pairing strength is set to zero, at zero temperature
    !       b) at high temperature after the pairing phase transition 
    !  *) Both HF and BCS formulas tend towards the zero temperature value as   
    !     beta gets increased.
    !  *) The formulas correctly produce 0 for spherical configurations at
    !     zero temperature.
    !---------------------------------------------------------------------------
    integer       :: i,j, b, it, ii, jj, si
    real(KIND=dp) :: ME(3), uvi, uvj, ui, uj, vi, vj, fi, fj, fac
    real(KIND=dp) :: wa, wb, wc, wd, Ba, Bb, dfde

    J2 = 0  ;  Belyaev = 0
 
    si = 0  
    do b = 1, Blocks
      do i=1, HFBlocks(b)
        ii = si + i
        it = 1
        if(ii.gt.nwn) it = 2
        do j=1,HFblocks(b)
          jj = si + j
        
          ! BCS --------------------------------------------------------------
          ! |< k | j_x | -l >|^2  
          ME(1)= angmom_xt_real(hfpsi(:,:,ii),hfpsi(:,:,jj),hfdpsi(:,:,:,jj))**2 
          ! |< k | j_y | -l >|^2 
          ME(2)= angmom_yt_imag(hfpsi(:,:,ii),hfpsi(:,:,jj),hfdpsi(:,:,:,jj))**2  
          ! |< k | j_z |  l >|^2 
          ME(3)= angmom_z_real( hfpsi(:,:,ii),hfpsi(:,:,jj),hfdpsi(:,:,:,jj))**2  

          ! The factor two is because of the presence of the four terms in 
          ! ME_{kl}, which are pair-wise equal when J_i = J_j.
          ME = 2 * ME 

          vi  = BCSOccupations(ii)/2. ; vj  = BCSOccupations(jj)/2.
          ui  = 1 - vi                ; uj  = 1 - vj
          fi  = BCSf(ii)              ; fj  = BCSf(jj)
          uvi = ui*vi                 ; uvj = uj*vj         

          ! Take the square root, but take care for numerical errors producing 
          ! small negative numbers.
          if(uvi .gt. 0) then
            uvi = sqrt(uvi) 
          else
            uvi = 0
          endif
          if(uvj .gt. 0) then
             uvj = sqrt(uvj)
          else
             uvj = 0
          endif

          !-------------------------------------------------------------------
          ! <J^2>
          ! [u^2_k u_l^2 + u_k v_k u_l v_l] * f_k (1-f_l)                  (a)
          wa = (ui*uj + uvi*uvj) * fi * (1-fj)
          ! [v^2_k v_l^2 + u_k v_k u_l v_l] * (1-f_k) f_l                  (b)
          wb = (vi*vj + uvi*uvj) * (1-fi) * fj
          ! [u_k^2 v_l^2 - u_k v_k u_l v_l] * f_k f_l                      (c)
          wc = (ui*vj - uvi*uvj) * fi * fj
          ! [v_k^2 u_l^2 - u_k v_k u_l v_l] * (1-f_k) (1-f_l)              (d)
          wd = (vi*uj - uvi*uvj) * (1-fi) *(1-fj)

          J2(:,it) = J2(:,it) + ME * (wa+wb+wc+wd)
          !-------------------------------------------------------------------
          ! I_xx, I_yy and I_zz
          if(abs(BCSqps(ii) - BCSqps(jj)).gt.1d-5) then            
            dfde = (fj - fi)/(BCSqps(ii) - BCSqps(jj))
          else
            if(inversetemp.gt.0) then
              ! beta exp(beta * E_k)* f_k^2
              dfde = inversetemp * exp(inversetemp*BCSqps(ii)) * fi**2
            else
              dfde = 0 ! All QPS are unoccupied for T=0 calculation
            endif
          endif
          ! (u_k u_l + v_k v_l)**2 * (f_l - f_k)/(E_k - E_l)  (Ba)
          Ba = (ui*uj + vi*vj + 2*uvi*uvj) * dfdE 
          
          dfdE = (1 - fi - fj)/(BCSqps(ii) + BCSqps(jj))

          ! (u_k v_l - v_k u_l)**2 * (1-f_k-f_l)/(E_k + E_l)  (Bb)
          Bb = (ui*vj + uj*vi - 2*uvi*uvj) * dfde
    
          Belyaev(:,it) =  Belyaev(:,it) + ME * (Ba + Bb)
        enddo
      enddo
      si = si +   HFBlocks(b)
    enddo
    !  Sum for the total
    J2(:,3) = sum(J2(:,1:2),2) ; Belyaev(:,3) = sum(Belyaev(:,1:2),2)
  end subroutine calcJ2andBelyaev_BCS

  subroutine calcJ2andBelyaev_HFB
    !---------------------------------------------------------------------------
    ! Calculate the expectation value of J^2_mu in the many-body state, as well
    ! as the Belyaev moment of inertia. 
    ! 
    ! No calculation for the Belyaev moment of inertia at finite temperature 
    ! right now.
    !---------------------------------------------------------------------------
   
    integer       :: i,j, b, it, ii, jj, si, N,k, kk, l, ll, sb
    real(KIND=dp) :: ME(3), uvi, uvj, ui, uj, vi, vj, fi, fj
    real(KIND=dp) :: wa, wb, wc, wd, Ba, Bb, dfde, fac, degen

    real(KIND=dp) :: jx(nwt,nwt), jy(nwt,nwt), jz(nwt,nwt)
    real(KIND=dp) :: jx_can(nwt,nwt), jy_can(nwt,nwt), jz_can(nwt,nwt)

    real(KIND=dp) :: J20(nwt,nwt, 3), J11(nwt,nwt,3)

    J2 = 0  ;  Belyaev = 0
    jx = 0  ; jy = 0 ; jz = 0

    !---------------------------------------------------------------------------
    ! First, we calculate the full matrix elements of jx, jy and jz
    si = 0  
    do b = 1, Blocks
      N = HFBlocks(b)
      do i=1, N
        ii = si + i
        it = 1
        if(ii.gt.nwn) it = 2
        do j=1,N
          jj = si + j
        
          ! |< k | j_x | -l >|^2            
          jx(ii,jj)= angmom_xt_real(hfpsi(:,:,ii),hfpsi(:,:,jj),hfdpsi(:,:,:,jj)) 
          ! |< k | j_y | -l >|^2 
          jy(ii,jj)= angmom_yt_imag(hfpsi(:,:,ii),hfpsi(:,:,jj),hfdpsi(:,:,:,jj))
          ! |< k | j_z |  l >|^2 
          jz(ii,jj)= angmom_z_real( hfpsi(:,:,ii),hfpsi(:,:,jj),hfdpsi(:,:,:,jj))  
        enddo
      enddo

      jx_can(si+1:si+N, si+1:si+N) = &
      &  matmul(transpose(cantransfo(si+1:si+N, si+1:si+N)), &
      &                              jx(si+1:si+N, si+1:si+N))
      jy_can(si+1:si+N, si+1:si+N) = &
      &  matmul(transpose(cantransfo(si+1:si+N, si+1:si+N)), &
      &                              jy(si+1:si+N, si+1:si+N))

      jz_can(si+1:si+N, si+1:si+N) = &
      & matmul(transpose(cantransfo(si+1:si+N, si+1:si+N)), &
      &                              jz(si+1:si+N, si+1:si+N))

      jx_can(si+1:si+N, si+1:si+N) = &
      & matmul(jx_can(si+1:si+N, si+1:si+N), &
      &        cantransfo(si+1:si+N, si+1:si+N))

      jy_can(si+1:si+N, si+1:si+N) = & 
      & matmul(jy_can(si+1:si+N, si+1:si+N), &
      &        cantransfo(si+1:si+N, si+1:si+N))

      jz_can(si+1:si+N, si+1:si+N) = &
      & matmul(jz_can(si+1:si+N, si+1:si+N), &
      &        cantransfo(si+1:si+N, si+1:si+N))

      si = si + N
    enddo

    ! First, construct J20
    call calcJ20(jx,bogoliubov, j20(:,:,1)) ! J20_x with an added T
    call calcJ20(jy,bogoliubov, j20(:,:,2)) ! J20_y with an added T
    call calcJ20(jz,bogoliubov, j20(:,:,3)) ! J20_z

    !Then , construct J11
    call calcJ11(jx,bogoliubov, j11(:,:,1), +1) ! J20_x with an added T
    call calcJ11(jy,bogoliubov, j11(:,:,2), +1) ! J20_y with an added T
    call calcJ11(jz,bogoliubov, j11(:,:,3), +1) ! J20_z
  
    !---------------------------------------------------------------------------
    ! Then we calculate the expectation value of J^2, in the canonical basis
    !---------------------------------------------------------------------------
    si = 0  
    do b = 1, Blocks
      N = HFBlocks(b)
      do i=1, N
        ii = si + i
        it = 1
        if(ii.gt.nwn) it = 2
        do j=1,N
          jj = si + j

          ! Factors 1./2 due to time-reversal
          fac=  rho_can(ii)/2.*(1-rho_can(jj)/2.)-kappa_can(ii)*kappa_can(jj)
          ME(1) = 2*jx_can(ii,jj)**2 
          ME(2) = 2*jy_can(ii,jj)**2 
          ME(3) = 2*jz_can(ii,jj)**2  

          J2(:,it) = J2(:,it) + ME * fac
        enddo
      enddo
      si = si + N
    enddo
    J2(:,3) = sum(J2(:,1:2),2) 

    !---------------------------------------------------------------------------
    ! I also calculate some approximation for the collective angular momentum, 
    ! which I define as <J^2> without the contribution from the blocked qps. 
    ! To safely remove this contribution, I calculate this in the qp basis.
    ! 
    !  <J^2> = sum_{ab} |J^20_ab|^2
    !
    ! where the sum simply does not include the blocked qps. 
    !---------------------------------------------------------------------------
    if(inversetemp .lt. 0) then
      si = 0 ; sb = 0
      J2_coll = 0
      do b=1,Blocks
        N = hfbsizes(b)
        do i=1, N
          ii = si + i
          it = 1
          if(ii.gt.nwn) it = 2

          do j=1,N
            jj = j + si
            
            ME = J20(ii,jj,:)**2             
            J2_coll(:,it) = J2_coll(:,it) + ME          
          enddo
        enddo
        si = si +   N
        sb = sb + 2*N
      enddo
      J2_coll(:,3) = sum(J2_coll(:,1:2), 2)  
    endif
    !---------------------------------------------------------------------------
    ! Then the Belyaev moment of inertia in the ordinary sp. basis.
    !
    ! From expanding the many-body wave-function around the HFB minimum for 
    ! small rotational frequency omega, we get the following expression
    !
    !  I_{mm} = 2 \sum_{ab} (E_a + E_b)^{-1} |J^{20}|^2_{m,ab}
    !
    ! based on pg 131 in Ring and Schuck, equation 3.92.
    !
    ! For a statistical mixture (such as an EFA configuration), this formula
    ! doesn't capture everything and we have to generalize:
    !
    !  I_mm = 2 \sum_{ab} (E_a + E_b)^{-1} |J^{20}|^2_{m,ab} (1 - f_a - f_b)
    !       + 2 \sum_{ab} (E_a - E_b)^{-1} |J^{11}|^2_{m,ab} (f_b - f_a)
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    si = 0 ; sb = 0
    do b = 1, Blocks
      N = HFBlocks(b)

      do i=1, N
        ii = si + i
        it = 1
        if(ii.gt.nwn) it = 2
          
        do j=1,N
          jj = si + j
          !---------------------------------------------------------------------
          !         1- f_i - f_j 
          fac = 1 - configmatrix(sb+i) - configmatrix(sb+j)

          Belyaev(:,it) = Belyaev(:,it) + &
          &           2*fac*J20(ii,jj,:)**2 /(Qpenergies(ii) + Qpenergies(jj))    

          !---------------------------------------------------------------------
          !      f_j - f_i
          fac =  configmatrix(sb+j) - configmatrix(sb+i)

          if(abs(Qpenergies(ii) - Qpenergies(jj)) .gt. 1d-8) then
            Belyaev(:,it) = Belyaev(:,it) + &
            &           2* fac*J11(ii,jj,:)**2 /(Qpenergies(ii)-Qpenergies(jj))   
          elseif(inversetemp .gt. 0) then
            degen = inversetemp * configmatrix(sb+i)**2 *                      &
            &                                  exp(inversetemp * Qpenergies(ii))
            Belyaev(:,it) = Belyaev(:,it) +  2*J11(ii,jj,:)**2 * degen   
          endif 
          !---------------------------------------------------------------------
        enddo
      enddo
      si = si +   N
      sb = sb + 2*N
    enddo
    Belyaev(:,3) = sum(Belyaev(:,1:2),2)

  end subroutine calcJ2andBelyaev_HFB  

  subroutine calcJ20(j,bogo, j20)
    !---------------------------------------------------------------------------
    ! Small subroutine to calculate the matrix J20 in the Ring and Schuck 
    ! notation.
    !
    !  J20 = U^\dagger j V^* - V^\dagger j^t U^*
    !
    !---------------------------------------------------------------------------
    ! Because of time-reversal symmetry, one should exercice caution when using 
    ! this routine. The full U and V matrices are 
    !
    !   ( U^+ 0  )    and  ( 0   -V ^+)
    !   ( 0   U^-)         ( V^+  0   )
    !
    ! and only the U^+ and V^+ are stored in memory. 
    !
    ! So, if the single-particle matrix elements passed in are just 
    !   < a | j_mu | b > this routine actually returns - J^{20}_{a\bar{b}} 
    ! where \bar{b} is the time-reversed partner of b.  
    ! 
    ! If instead, the single-particle matrix elements that are passed in are
    !   < a | j_mu T | b >, the result of this routine 
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in)  :: j(nwt,nwt)
    real(KIND=dp), intent(in)  :: bogo(2*nwt, 2*nwt)
    real(KIND=dp), intent(out) :: j20(nwt,nwt)
    real(KIND=dp)              :: tmp(nwt,nwt)
    integer :: N, si, sb, b

    si = 0 ; sb = 0
    do b = 1, Blocks
      N = HFBlocks(b)

      ! U^\dagger j 
      j20(si+1:si+N, si+1:si+N) = &
      & matmul(transpose(bogo(sb+1:sb+N, sb+N+1:sb+2*N)),j(si+1:si+N,si+1:si+N))

      ! U^\dagger j V^*
      j20(si+1:si+N, si+1:si+N) = &
      & matmul( j20(si+1:si+N, si+1:si+N), (bogo(sb+N+1:sb+2*N, sb+N+1:sb+2*N)))

      ! V^\dagger j^t 
      tmp(si+1:si+N, si+1:si+N) = transpose(j(si+1:si+N, si+1:si+N))
      tmp(si+1:si+N, si+1:si+N) = &
      & matmul(transpose(bogo(sb+N+1:sb+2*N, sb+N+1:sb+2*N)),                  &
      &                             tmp(si+1:si+N,si+1:si+N))
      
      ! V^\dagger j^T U
      j20(si+1:si+N, si+1:si+N) = j20(si+1:si+N, si+1:si+N) - &
      &  matmul( tmp(si+1:si+N, si+1:si+N), (bogo(sb+1:sb+N, sb+N+1:sb+2*N)))

      si = si +   N
      sb = sb + 2*N
    enddo

  end subroutine calcJ20

  subroutine calcJ11(j,bogo, j11, s)
    !---------------------------------------------------------------------------
    ! Small subroutine to calculate the matrix J11 in the Ring and Schuck 
    ! notation.
    !
    !  J11 = U^\dagger j U - V^\dagger j^t V^*
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in)  :: j(nwt,nwt)
    real(KIND=dp), intent(in)  :: bogo(2*nwt, 2*nwt)
    real(KIND=dp), intent(out) :: j11(nwt,nwt)
    real(KIND=dp)              :: tmp(nwt,nwt)
    integer :: N, si, sb, b, s

    si = 0 ; sb = 0
    do b = 1, Blocks
      N = HFBlocks(b)

      ! U^\dagger j 
      j11(si+1:si+N, si+1:si+N) = &
      & matmul(transpose(bogo(sb+1:sb+N, sb+N+1:sb+2*N)),j(si+1:si+N,si+1:si+N))

      ! U^\dagger j U
      j11(si+1:si+N, si+1:si+N) = &
      & matmul( j11(si+1:si+N, si+1:si+N), (bogo(sb+1:sb+N, sb+N+1:sb+2*N)))

      ! V^\dagger j^t 
      tmp(si+1:si+N, si+1:si+N) = transpose(j(si+1:si+N, si+1:si+N))
      tmp(si+1:si+N, si+1:si+N) = &
      & matmul(transpose(bogo(sb+N+1:sb+2*N, sb+N+1:sb+2*N)),                  &
      &                             tmp(si+1:si+N,si+1:si+N))
      
      ! V^\dagger j^T V
      j11(si+1:si+N, si+1:si+N) = j11(si+1:si+N, si+1:si+N) + s * &
      &  matmul( tmp(si+1:si+N, si+1:si+N),(bogo(sb+N+1:sb+2*N, sb+N+1:sb+2*N)))

      si = si +   N
      sb = sb + 2*N
    enddo
  end subroutine calcJ11

  subroutine PrintMomentsofInertia()
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    1 format (21('-'), 'Rotational properties', 20('-'))
    2 format ('                Rigid rotor (hbar^2/MeV)')
    3 format ('              n', 14x ,'p',14x,'t ')
    4 format (' I_R X ', 3f15.7)
    5 format (' I_R Y ', 3f15.7)
    6 format (' I_R Z ', 3f15.7)
    7 format ('                Belyaev     (hbar^2/MeV) ')
    8 format (' I_B X ', 3f15.7)
    9 format (' I_B Y ', 3f15.7)
   10 format (' I_B Z ', 3f15.7)
   11 format ('                 J^2        (hbar^2)')
   12 format (' J2_X  ', 3f15.7)
   13 format (' J2_Y  ', 3f15.7)
   14 format (' J2_Z  ', 3f15.7)
   15 format (' J2_t  ', 3f15.7)
   16 format ('                 J^2_coll   (hbar^2)')

  100 format (60('-'))

    print 1
    print 2
    print 3
    print 4, Rigid(1,:)
    print 5, Rigid(2,:)
    print 6, Rigid(3,:)
    print *
    print 7
    print 3
    print  8, Belyaev(1,:)
    print  9, Belyaev(2,:)
    print 10, Belyaev(3,:)
    print *
    print 11
    print 12, J2(1,:)
    print 13, J2(2,:)
    print 14, J2(3,:) 
    print 15, sum(J2,1)
    print *
    print 16
    print *
    print 12, J2_coll(1,:)
    print 13, J2_coll(2,:)
    print 14, J2_coll(3,:)
    

  end subroutine PrintMomentsofInertia
end module momentsofinertia
