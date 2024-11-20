module fission_MOI
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
 ! Module governing the calculation of the collective inertia tensor for 
 ! fission calculations. 
 !
 ! Currently implemented: 
 !  (i) perturbative cranking approximation
 !
 ! References used for the construction of this module:
 !  (a) S. A. Giuliani and L. M. Robledo,
 !      Non-perturbative collective inertias for fission: a comparative study
 !      Physics Letters B 787, 134-140 (2018).
 !
 !  (b) A. Baran, J. A. Sheikh, J. Dobaczewski, W. Nazarewicz and A. Staszczak
 !      Quadrupole collective inertia in nuclear fission: cranking approximation
 !      Phys. Rev. C 84, 054321 (2011)
 !
 ! Attention: to the best of my understanding, there are errors in Baran et al. 
 !            that lead to almost, but not quite, the same results! The final
 !            implementation here is based on S.A. Giuliani et al.
 ! 
 !------------------------------------------------------------------------------
 ! Hephaestos keywords
 ! 
 ! PBROKEN : $PBROKEN
 !
 !------------------------------------------------------------------------------

  use geninfo
  use densities
  use parameterization
  use moments
  use timing

  implicit none

  !-----------------------------------------------------------------------------
  ! Multipole moments for which to construct the inertia tensor. 
  integer :: N_inertia            = 0
  integer, allocatable :: inertia_l(:) 
  integer, allocatable :: inertia_m(:)  

  !-----------------------------------------------------------------------------
  ! The full collective inertia tensor, obtained by including information 
  ! on ALL the multipole moments that were asked for  
  real(KIND = dp), allocatable :: collective_inertia(:,:)
  ! Intermediate matrices M^1 and M^3 that are needed for the calculation
  ! of the collective_inertia. Stored separately so it can be output for 
  ! people wanting to recalculate the collective inertia.
  real(KIND = dp), allocatable :: M1(:,:,:), M3(:,:,:)
  !-----------------------------------------------------------------------------
  
contains 

  subroutine read_inertia(file_number)
    !---------------------------------------------------------------------------
    ! Subroutine to read the &inertia/ namelist from the specified file (via the
    ! specified channel) or from STDIN if the variables are not present.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   file_number : optional integer. If present, read from (open) channel
    !                 with this number. If absent, read from STDIN.
    !---------------------------------------------------------------------------
    integer(dp), intent(in), optional   :: file_number   
    integer :: k
#if(USE_MPI>0)
    integer :: mpi_err
#endif

    NameList /inertia/ inertia_l, inertia_m

    ! Sanity check    
    if(N_inertia .lt. 0) then
      call stp('N_inertia cannot be negative.')
    else if (N_inertia .eq. 0) then
      ! do nothing
      return
    endif
  
    allocate(inertia_l(N_inertia)) ; inertia_l = -1
    allocate(inertia_m(N_inertia)) ; inertia_m = -1

    ! only the very first MPI rank reads input   
    if(MPI_RANK.eq.0) then 
      if(present(file_number)) then
        read (unit=file_number, nml=inertia)
      else
        read (unit=*, nml=inertia)
      endif

      ! Some sanity checks
      do k=1, N_inertia
        if(inertia_l(k) .eq. -1) then 
          call stp('Number of elements in inertia_l does not match N_inertia.')
        endif
        
        if(inertia_l(k) .gt. maxmoment) then
          call stp('Cannot compute inertia for Qlm with l > Maxmoment.')
        endif
        
        if(inertia_m(k) .eq. -1) then 
          call stp('Number of elements in inertia_m does not match N_inertia.')
        endif

        if(inertia_m(k) .gt. inertia_l(k)) then
          call stp('Cannot compute inertia for Qlm with m > l.')
        endif
      enddo
    endif
    
    ! ... and then broadcast to all ranks
#if(USE_MPI > 0)
    call MPI_Bcast(inertia_l, N_inertia, MPI_INTEGER,0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(inertia_m, N_inertia, MPI_INTEGER,0, MPI_COMM_WORLD, mpi_err)
#endif      
  end subroutine read_inertia

  subroutine print_collective_inertia()
    !---------------------------------------------------------------------------
    ! Print all entries in the collective inertia tensor.
    !---------------------------------------------------------------------------

    1 format (26('-'),' Collective inertia tensor', 27('-'))
    2 format (80('-'))
    3 format ('          I_Q',2i1, 2x)
    4 format ('     I_Q',2i1, 1x,'|', 1x, 99es15.5)
    6 format (15('-'))
   31 format ('         M1_Q',2i1, 2x)
   41 format ('    M1_Q',2i1, 1x,'|', 1x, 99es15.5)
   32 format ('         M3_Q',2i1, 2x)
   42 format ('    M3_Q',2i1, 1x,'|', 1x, 99es15.5)

    
   99 format ('  Conventions:' /, & 
   &          '    Collective variables: multipole moments Qlm = r^l Y_lm ,',/,&
   &          '                          calculated wrt COM = (', 3f8.3, ')',/,&
   &          '    Collective modes normalized with hbar = 1.'              /,&
   &          '    Units of collective inertias in MeV^{-1} b^{-l} [hbar^2].')

    character(len=120) :: header, sep
    character(len=16) :: tmp
    real(KIND=dp)     :: shiftx, shifty, shiftz
    integer :: i

    
    print 1
    print *
    shiftx = -meshx_shifted(1) + meshx(1)      
    shifty = -meshy_shifted(1) + meshy(1)      
    shiftz = -meshz_shifted(1) + meshz(1)      
    print 99, shiftx, shifty, shiftz
    print *

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Constructing the right size separator line
    sep = ''
    tmp = ''
    write(sep,6)    
    do i=1,N_inertia
      write(tmp,6)    
      sep = adjustl(trim(sep)//tmp)  
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Collective inertia tensor
    header = ''
    do i=1,N_inertia
      write(tmp, 3) inertia_l(i),inertia_m(i) 
      header = adjustl(trim(header)//tmp)
    enddo

    print *
    print *, '                ', header
    print *,sep
    do i=1, N_inertia
      print 4, inertia_l(i),inertia_m(i), collective_inertia(i,1:N_inertia)
    enddo
    print *,sep
    print *

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Collective inertia tensor, just the M1 matrix
    header = ''
    do i=1,N_inertia
      write(tmp, 31) inertia_l(i),inertia_m(i) 
      header = adjustl(trim(header)//tmp)
    enddo

    print *, ' neutrons        ', header   
    print *, sep 
    do i=1, N_inertia
      print 41, inertia_l(i),inertia_m(i), M1(i,1:N_inertia,1)
    enddo
    print *,sep
    print *
    print *, ' protons         ', header   
    print *,sep
    do i=1, N_inertia
      print 41, inertia_l(i),inertia_m(i), M1(i,1:N_inertia,2)
    enddo
    print *,sep
    print *
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Collective inertia tensor, just the M3 matrix
    header = ''
    do i=1,N_inertia
      write(tmp, 32) inertia_l(i),inertia_m(i) 
      header = adjustl(trim(header)//tmp)
    enddo

    print *, ' neutrons        ', header   
    print *, sep 
    do i=1, N_inertia
      print 42, inertia_l(i),inertia_m(i), M3(i,1:N_inertia,1)
    enddo
    print *,sep
    print *
    print *, ' protons         ', header   
    print *,sep
    do i=1, N_inertia
      print 42, inertia_l(i),inertia_m(i), M3(i,1:N_inertia,2)
    enddo
    print *,sep
    print *
        
    print 2
  end subroutine print_collective_inertia
  
  subroutine verify_COM_motion()
    !---------------------------------------------------------------------------
    !  This routine performs calculations of the collective moments of inertia
    !  for the motion of the z-coordinate of the center of mass.
    !
    !  The collective coordinate for species q is thus
    !    Q_q = z_q / A
    !  to which corresponds a collective momentum  (in our convention)     
    !    P_q = - i \nabla_z (*)
    !  
    !  One can show that, analytically, the collective inertia associated with 
    !  movement of the centre-of-mass should be the TOTAL mass of the nucleus
    ! 
    !   M'_{0} = A m
    !
    !  This is what is always presented in the literature. Note however the 
    !  little accent M', indicating that this IS NOT the collective inertia
    !  associated with P_q. Rather, it is the collective mass associated with 
    !     P'_q = hbar P_q
    !  i.e. the 'physical' convention for the momentum. 
    !
    ! We calculate the collective inertia related to COM motion here in multiple
    ! ways:
    !
    !   (a) Analytically, printed as 'A m'
    !   (b) By using the Belyaev formula for the momentum in our convention (*)
    !       (and multiplying by hbar^2 afterward)
    !   (c) By using the perturbative cranking formula starting from Q_q
    !       (and converting convention again afterward)
    !
    ! The results of (b) and (c) are not necessarily close to (a) however: 
    ! for typical Skyrme interactions the effective mass m^*/m is not equal
    ! to one, spoiling the correspondence of the perturbative formulation. 
    ! The origin lies in the absence of Galileian invariance of the interaction 
    ! in the perturbative treatment, see 
    !
    !    K. Wen, and T. Nakatsukasa,  http://arxiv.org/abs/2112.13317
    ! 
    ! for a discussion. Ideally, we would calculate an "average" effective mass
    ! as these authors do and use it to correct our results.
    !
    !---------------------------------------------------------------------------
    
    use functional
    use densities
    
    1 format (' Pushing model                 M_0 (MeV/c^2)')
    2 format (25x, 'neutrons        protons          total', / &
    &         1x, 80('-'))
    3 format ('   Belyaev M_0      | ', 3f16.5)
    4 format ('   Pert. cranking   | ', 32x,f16.5)
    7 format (1x, 80('-'),/, &
    &         '   Am               | ', 3f16.5, /, 80('-'))
    
    real(KIND=dp), allocatable :: NablaMElements(:,:,:,:)
    real(KIND=dp) :: mat(2,2), Pmat(2,2), totalmass
      
    real(KIND=dp) :: Psp(nwt,nwt), Qsp(nwt,nwt), hbar
    real(KIND=dp) :: P20(nwt,nwt), Q20(nwt,nwt)
    
    ! Calculate hbar to make its use consistent
    hbar =  sqrt(hbm(1) * 2  * 0.5 * sum(nucleonmass))
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculate single-particle matrix elements of nabla_z 
    NablaMElements = compNablaMElements()
    Psp            =  NablaMElements(3,1,:,:) 
    if(pairingtype.eq.2) then
      ! Attention, in a HFB calculation these are the matrix elements in the 
      ! canonical basis, while Bogoliubov refers to the HF basis. So, we 
      ! transfer back to the HFbasis
      Psp = matmul(matmul(cantransfo, Psp), transpose(cantransfo))
    endif
    
    Psp(1:nwn,1:nwn)         =  Psp(1:nwn,1:nwn)         * hbar
    Psp(nwn+1:nwt,nwn+1:nwt) =  Psp(nwn+1:nwt,nwn+1:nwt) * hbar
        
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculate single-particle matrix elements of z-c.o.m. coordinate
    Qsp            = sqrt(4*pi/3) * Qlm_spme(1,0,.false.) * 10 
                    ! Q10 = sqrt(3/4pi) * z 
                    ! and the routine Qlm_spme uses units of b^1/2 => factor 10
    Qsp(1:nwn,1:nwn)         =      Qsp(1:nwn,1:nwn)        /(neutrons+protons)
    Qsp(nwn+1:nwt,nwn+1:nwt) =      Qsp(nwn+1:nwt,nwn+1:nwt)/(protons +neutrons)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TODO: Calculate the average effective mass for the pushing model.
    !       This used to work, but now fails thanks to the changes of the 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Calculate average effective mass 
    !
    !   m^*/m_q  = 2m/N_q * int d^3r   rho_q(r) hbar^2/2m^*_q(r)
    ! 
    ! with hbar^2/2m^*_q =  hbar^2/2m_q(r) + F_Nm_Nm(r)
    ! 
    ! which is the logical generalization from Eq. (30) in 
    !    K. Wen, and T. Nakatsukasa,  http://arxiv.org/abs/2112.13317
    ! and an explicit factor of hbar^2.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     avg_effmass = sum(F%F_Nm_Nm(:,:) * R%D_I_I(:,:) , 1) * dv
    !     avg_effmass(1) = avg_effmass(1) / (hbm(1) * neutrons)
    !     avg_effmass(2) = avg_effmass(2) / (hbm(2) * protons)
    !     avg_effmass    = avg_effmass + 1
    !     avg_effmass    = 1/avg_effmass
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    print 1
    print 2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! perform the summations
    select case(pairingtype)
    case(0) ! HF
      Pmat                     =  Ksum_Mij_HF(Psp, Psp, 1, 1,(/1,3/))
      mat                      =  Ksum_Mij_HF(Qsp, Qsp, 1, 1,(/1,3/))
    case(1) ! BCS
      Pmat                     =  Ksum_Mij_BCS(Psp, Psp, 1, 1,(/1,3/))
      mat                      =  Ksum_Mij_BCS(Qsp, Qsp, 1, 1,(/1,3/))
    case(2) !HFB
      P20 = calc_Q20(Psp, bogoliubov, 1)        
      Pmat                     =  Ksum_Mij(P20, P20, 1, 1,(/1,3/))
      Q20 = calc_Q20(Qsp, bogoliubov, 1)        
      mat                      =  Ksum_Mij(Q20, Q20, 1, 1,(/1,3/))
    end select
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculating the inertia parameters
    totalmass   = 1.0d0/(sum(mat(1,:))) * sum(mat(2,:)) * 1.0d0/(sum(mat(1,:)))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Printing 
    print 3,  Pmat(1,:), sum(Pmat(1,:))
    print 4,  totalmass * hbar**2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Analytical result
    print 7,  neutrons * nucleonmass(1), protons*nucleonmass(2), &
    &         neutrons * nucleonmass(1)+ protons*nucleonmass(2) 

    print *

  end subroutine verify_COM_motion

  subroutine calc_collective_inertia()
    !---------------------------------------------------------------------------
    ! Calculate the collective inertia tensor for the set of multipole moments
    ! Qlm determined in inertia_l and inertia_m. This tensor in the perturbative
    ! cranking approximation is given by
    !
    !     M_c =  M_1^{-1} M_3 M_1^{-1}    (from Giuliani and Robledo)
    !
    ! where the matrices M_n are determined by
    !
    !                          Q^{20}_{i,ab} Q^{20}_{j,ab} 
    ! M_{n,ij} = Re sum_{ab}  ---------------------------
    !                               (E_a + E_b)^n
    !
    ! where the sum is over all quasiparticle states and E_a and E_b are 
    ! quasiparticle energies. The collective coordinates Q_{i} are given by
    ! the multipole moments Q_lm as defined in the arrays inertia_l 
    ! and inertia_m.
    !
    ! Steps:
    !  (1) Calculate all the single-particle matrix elements of the Qlm
    !      in the HF-basis with routine Qlm_spme
    !  (2) Transform these matrix elements to the qp basis
    !      with routine calc_Q20
    !  (3) Sum the matrix elements, weighted with the appropriate power of 
    !      the quasiparticle energies, using Ksum_Mij
    !  (4) Invert M_1 with Lapack routines
    !  (5) Obtain M_c for each species. The total inertia is M_t = M_n + M_p
    !          
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Note: this routine is not yet ready to deal with blocked HFB vacua!
    !
    !---------------------------------------------------------------------------
    ! Explicit declaration of the external linear algebra routines
    external :: dsytrf, dsytri

    real(KIND=dp), allocatable :: Mat(:,:,:,:), Qsp(:,:,:), Q20(:,:,:)
    real(KIND=dp), allocatable :: work(:), M1_inv(:,:)
    integer :: i, j, la, lb, l, m, info, lwork
    integer, allocatable :: ipiv(:)
        
    call start_timer(T_collective_MOI)
        
    if(.not.allocated(collective_inertia)) then
      allocate(collective_inertia(N_inertia, N_inertia))
    endif
    collective_inertia = 0
    
    allocate(Mat(N_inertia, N_inertia, 2,2)) ;  Mat   = 0.0d0
    allocate(Qsp(nwt,nwt,N_inertia))         ;  Qsp = 0.0d0
    
    if(pairingtype.eq.2) then
      allocate(Q20(nwt,nwt,N_inertia))         ;  Q20 = 0.0d0
    endif  
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Step 1 & 2: construct relevant sp matrices
    do i=1, N_inertia
      l = inertia_l(i)
      m = inertia_m(i)

      ! Calculate all relevant single-particle matrix elements              
      Qsp(:,:,i) = Qlm_spme(l,m,.false.) ! Hardcoded to consider only real parts
                                         ! at the moment
      if(pairingtype.eq.2) then
        ! Transform to the quasiparticle basis if needed
        Q20(:,:,i) = calc_Q20(Qsp(:,:,i), bogoliubov, l)
      endif
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Step 3: we calculate the sums for every combination of collective DOF
    do i=1,N_inertia
      la = inertia_l(i)
      do j=1, N_inertia
        lb = inertia_l(j)
        ! Perform the sums to obtain M_k for k=1,3
        select case(pairingtype)
        case(0)
          ! HF summation
          Mat(i,j,:,:) = Ksum_Mij_HF(Qsp(:,:,i), Qsp(:,:,j), &
          &                                                      la, lb,(/1,3/))
        case(1)
          ! BCS summation
          Mat(i,j,:,:) = Ksum_Mij_BCS(Qsp(:,:,i), Qsp(:,:,j), &
          &                                                      la, lb,(/1,3/))
        case(2)
          ! HFB summation
          Mat(i,j,:,:) = Ksum_Mij(Q20(:,:,i), Q20(:,:,j), la, lb,  (/1,3/))
        end select
      enddo
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Constructing explicitly the matrices M_1 and M_3 for ease of reading
    if(.not. allocated(M1)) allocate(M1(N_inertia, N_inertia,3))
    if(.not. allocated(M3)) allocate(M3(N_inertia, N_inertia,3))
    M1(:,:,1:2) = Mat(:,:,1,1:2) ; M1(:,:,3) = sum(M1(:,:,1:2),3)
    M3(:,:,1:2) = Mat(:,:,2,1:2) ; M3(:,:,3) = sum(M3(:,:,1:2),3)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Step 4: use LAPACK routines to invert M1
    allocate(M1_inv(N_inertia, N_inertia))
    M1_inv = M1(:,:,3)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Ask for a workspace size
    allocate(work(1), ipiv(N_inertia))
    call dsytrf('U', N_inertia, M1_inv,N_inertia, ipiv, work,-1, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Factorize M1
    call dsytrf('U', N_inertia, M1_inv,N_inertia,ipiv,work,lwork,info)
    ! Invert M1
    if(MAXVAL(abs(M1_inv)).gt.1e-15) then
      call dsytri('U', N_inertia, M1_inv, N_inertia,ipiv,work, info)
      if(info.ne.0) then
        call stp('Problem for DSYTRI during the calculation of collective inertia.')
      endif
    endif
    deallocate(work, ipiv)

    ! Note that after DSYTRI, only the top half of M1 is guaranteed to be right
    ! Thus, we populate the other half here to avoid any surprises
    do i=1,N_inertia
      do j=i+1,N_inertia
        M1_inv(j,i) = M1_inv(i,j)
      enddo
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Step 5: calculate cranking tensor for every isospin
    !             M_c = 1/4 M1^{-1} M3 M1^{-1}
    !         and sum the results
    !             M_t = M_n + M_p
    collective_inertia = matmul(matmul(M1_inv, M3(:,:,3)), M1_inv)

    deallocate(Mat, M1_inv, Qsp)
    if(pairingtype.eq.2) deallocate(Q20)
    call stop_timer(T_collective_MOI)

  end subroutine calc_collective_inertia
  
  function Ksum_Mij_HF(Qa, Qb, la, lb, Ks) result (Ksum)
    !---------------------------------------------------------------------------
    ! Perform the sum over holes and particles for the matrices of the form
    !
    !                              Q_{i,ml} Q_{j,lm} 
    ! M_{n,ij} = 2 sum_m sum_l  ---------------------------
    !                               (e_m - e_l)^n
    ! 
    ! where Q_{i/j, ml} are single-particle matrix elements and the e_m/l
    ! are single-particle energies. The m are empty sp states, the l are 
    ! occupied states.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     Qa, Qb    : two-quasiparticle representations of both multipole 
    !                 operators
    !     Ks        : set of powers to use in the inverted calculation
    ! Output:
    !     Ksum      : result of the summations, array with the size of Ks
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Qa(:,:), Qb(:,:)
    real(KIND=dp), allocatable:: Ksum(:,:)
    integer, intent(in)       :: Ks(:), la, lb
 
    integer :: i, j, k, it, itb , Nk, trash
    real(KIND=dp) :: num, denom 
 
    Nk = size(Ks)
    allocate(Ksum(Nk,2)) ; Ksum = 0.0d0

    ! If parity is conserved, there is a parity selection rule    
$PCONSERVED if(mod(la,2) .ne. mod(lb,2)) return    
    if(la .eq. lb) trash = 0           
    !         ^---- trash statement to stop compilator complaints
    
    do i=1,nwt
      ! Loop over full HF states
      if(abs(rho_can(i)) .lt. 0.5) cycle
     
      it = 1 ;  if(i.gt. nwn) it = 2
      do j=1,nwt
        ! Loop over empty HF states
        if(abs(rho_can(j)) .gt. 0.1) cycle
        
        itb = 1 ;  if(j.gt. nwn) itb = 2
        if(it .ne. itb) cycle

        num   = Qa(i,j) * Qb(i,j) 
        do k=1,Nk
         denom      = (spenergies(j) - spenergies(i))**Ks(k)
         Ksum(k,it) = Ksum(k,it) + num/denom
        enddo
      enddo
    enddo    
    
    ! Factor two for the time-reversal partners
$TR    Ksum = 2*Ksum

    ! Explicit factor 2 Eq. 3.89
    Ksum = 2*Ksum
  end function Ksum_Mij_HF
  
  function Ksum_Mij(Qa, Qb, la, lb,  Ks) result(Ksum)
    !---------------------------------------------------------------------------
    ! Perform the relevant sums over the quasiparticle space for the calculation
    ! of the collective inertia, i.e.
    ! 
    !  M_k = sum_ij (<0| Qa | ij >< ij | Qb | 0 >)/[(Ei + Ej)^k]
    !
    ! This routine bunches the summations for the same multipole moments with
    ! all different powers k.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     Qa, Qb    : two-quasiparticle representations of both multipole 
    !                 operators
    !     Ks        : set of powers to use in the inverted calculation
    ! Output:
    !     Ksum      : result of the summations, array with the size of Ks
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Qa(:,:), Qb(:,:)
    real(KIND=dp), allocatable:: Ksum(:,:)
    integer, intent(in)       :: Ks(:), la, lb
    
    integer :: Na, N2a, Ta, Ba, Bb, Tb, Nb, N2b, ita, itb, k, i,j, Nk
    integer :: sai, sbi, sab, sbb, trash
    real(KIND=dp) :: num, denom
    
    logical :: blocked

    Nk = size(Ks)
    allocate(Ksum(Nk,2)) ; Ksum = 0.0d0

    ! If parity is conserved, there is a parity selection rule    
$PCONSERVED if(mod(la,2) .ne. mod(lb,2)) return     
    if(la .eq. lb) trash = 0           
    !         ^---- trash statement to stop compilator complaints

    sai = 0 ; sab = 0
    do Ba=1,8,2
      Na = HFBlocks(Ba) ; if(Na.eq.0) cycle
      N2a= HFBlocks(Ba+1)
      Ta = Na + N2a
      ita = 1 ; if(Ba.gt.4) ita=2
      
      sbi = 0 ; sbb = 0
      do Bb=1,8,2
        Nb = HFBlocks(Bb) ; if(Nb.eq.0) cycle
        N2b= HFBlocks(Bb+1)
        Tb = Nb + N2b
        itb= 1 ; if(Bb.gt.4) itb=2
        
        !Gain some CPU time
        if(ita.ne.itb) then
          sbi = sbi +   Tb
          sbb = sbb + 2*Tb
          cycle
        endif
        
        do i=1,Ta
          do j=1,Tb
            
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Don't include the contribution from the blocked qps and their
            ! partner qps.
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            blocked= .false.
            if(allocated(blocked_qps)) then
              do k=1,size(blocked_qps)
                if(sai+i.eq.blocked_qps(k) .or.  sbi+j.eq.blocked_qps(k)) then
                  blocked = .true.
                endif
                if(sai+i.eq.partner_qps(k) .or.  sbi+j.eq.partner_qps(k)) then
                  blocked = .true.
                endif
              enddo              
            endif
            if(blocked) cycle
          
            num = Qa(sai+i,sbi+j) * Qb(sai+i,sbi+j)
            do k=1,Nk
              denom = (qpenergies(sab+Ta+i) + qpenergies(sbb+Tb+j))**Ks(k)
              Ksum(k,ita) = Ksum(k,ita) + num/denom
            enddo
          enddo
        enddo
        sbi = sbi +   Tb
        sbb = sbb + 2*Tb
      enddo
      sai = sai +   Ta
      sab = sab + 2*Ta
    enddo
    
    ! Factor two for the time-reversal partners
$TR    Ksum = 2*Ksum
    
  end function Ksum_Mij
    
  function Ksum_Mij_BCS(Qa, Qb, la, lb,  Ks) result(Ksum)
    !---------------------------------------------------------------------------
    ! Perform the relevant sums over the quasiparticle space for the calculation
    ! of the collective inertia in the case of a BCS calculation, i.e.
    ! 
    !  M_k = sum_ij (<i|Qa|j><j|Q^\dagger_b|i >)/[(Ei + Ej)^k] {eta^+_ij}^2
    !
    ! where 
    !  (*) the <|Q|> are single-particle matrix elements
    !  (*) the Ei and Ej are a BCS quasiparticle energies
    !  (*) eta^+_ij = u_i v_j + u_j v_i 
    !  (*) the sum runs over all single-particle states
    !
    ! This expression is Eq. 58 in 
    ! 
    !  A. Baran et al,  Phys. Rev. C 84, 054321 (2011).
    !   
    ! This routine bunches the summations for the same multipole moments with
    ! all different powers k.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     Qa, Qb    : two-quasiparticle representations of both multipole 
    !                 operators
    !     Ks        : set of powers to use in the inverted calculation
    ! Output:
    !     Ksum      : result of the summations, array with the size of Ks
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Qa(:,:), Qb(:,:)
    real(KIND=dp), allocatable:: Ksum(:,:)
    integer, intent(in)       :: Ks(:), la, lb
    real(KIND=dp)             :: num, denom, eta, ui, vi, uj, vj
    integer                   :: i,j, it, itb, k, Nk, trash

    Nk = size(Ks)
    allocate(Ksum(Nk,2)) ; Ksum = 0.0d0

    ! If parity is conserved, there is a parity selection rule    
$PCONSERVED if(mod(la,2) .ne. mod(lb,2)) return     
    if(la .eq. lb) trash = 0           
    !         ^---- trash statement to stop compilator complaints

    do i=1,nwt
      call uv_from_occupation(BCSoccupations(i), ui, vi)
     
      it = 1 ;  if(i.gt. nwn) it = 2
      do j=1,nwt
        call uv_from_occupation(BCSoccupations(j), uj, vj)
        
        itb = 1 ;  if(j.gt. nwn) itb = 2
        if(it .ne. itb) cycle

        eta   = ui * vj + vi * uj  
        num   = Qa(i,j) * Qb(i,j) * eta**2
        do k=1,Nk
         denom      = (BCSqps(i) + BCSqps(j))**Ks(k)
         Ksum(k,it) = Ksum(k,it) + num/denom
        enddo
      enddo
    enddo    
    ! Factor two for the time-reversal partners
    Ksum = 2*Ksum
    
  end function Ksum_Mij_BCS
  
  function calc_Q20(Qsp, bogo, l) result(Q20)
    !---------------------------------------------------------------------------
    ! Function that calculates the two-quasiparticle matrix for a multipole 
    ! moment operator in the case of a HFB calculation.
    ! 
    !  Q20 = U^\dagger Q V^* - V^\dagger Q^t U^*
    !
    ! where Q is the matrix of single-particle matrix elements of the 
    ! multipole operator.
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     Qsp       : single-particle matrix elements 
    !     Bogo      : Bogoliubov transformation from the sp basis to the 
    !                 quasiparticle basis
    !     l         : ell of the multipole moment
    ! Output:
    !     Q20       : two-quasiparticle component of Qlm
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Note: 
    !
    !  * This routine assumes that the single-particle matrix elements 
    !    Qsp are (i) real and (ii) correspond to an operator that conserves
    !    signature symmetry. 
    ! 
    !  * If time-reversal is conserved, the routine assumes IN ADDITION that 
    !    the operator conserves time-reversal symmetry as well, i.e.
    !                 Q_ab = Q_{\bar{a} \bar{b}}
    !    where \bar{x} indicates a time-reversal partner. 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: Qsp(nwt,nwt), bogo(2*nwt,2*nwt)
    integer, intent(in)        :: l
    real(KIND=dp)              :: Q20(nwt,nwt), tmp(nwt,nwt)
    real(KIND=dp), allocatable :: U(:,:), V(:,:)    
    integer                    :: N, N2, T, si, sb, B, Tp, Tm

    Q20 = 0.0d0

    !---------------------------------------------------------------------------
    ! This situation is:
    !  (a) it concerns a multipole moment with even parity; or
    !  (b) parity is broken 
    ! In both cases we can simply multiply matrices straightforwardly.
$PCONSERVED    if(mod(l,2).eq.0) then
      si = 0 ; sb = 0
      do B=1,8,2
        N = HFBlocks(B) ; if(N.eq.0) cycle
        N2= HFBlocks(B+1)
        T = N + N2

        U = bogo(sb+  1:sb+  T,sb+ T+1:sb+2*T)
        V = bogo(sb+T+1:sb+2*T,sb+ T+1:sb+2*T)
        
        Q20(si+1:si+T,si+1:si+T) = matmul(transpose(U),Qsp(si+1:si+T,si+1:si+T))
        ! U^dagger Q V^*
        Q20(si+1:si+T,si+1:si+T) = matmul(Q20(si+1:si+T,si+1:si+T), V)
        
        ! Q V^* 
        tmp(si+1:si+T,si+1:si+T) = matmul(Qsp(si+1:si+T,si+1:si+T), V)
        ! V^\dagger Q^t
        tmp(si+1:si+T,si+1:si+T) = transpose(tmp(si+1:si+T,si+1:si+T))
        ! V^dagger Q^t U^*
        tmp(si+1:si+T,si+1:si+T) = matmul(tmp(si+1:si+T,si+1:si+T), U)

        ! Put both parts together      
        ! Default formula
$NTR        Q20(si+1:si+T,si+1:si+T) = + Q20(si+1:si+T,si+1:si+T) &
$NTR        &                          - tmp(si+1:si+T,si+1:si+T)

        ! Attention to the extra sign incurred when time-reversal is 
        ! conserved
$TR        Q20(si+1:si+T,si+1:si+T) = - Q20(si+1:si+T,si+1:si+T) &
$TR        &                          - tmp(si+1:si+T,si+1:si+T)

        si = si +   N +   N2
        sb = sb + 2*N + 2*N2
      enddo
$PCONSERVED    endif
    
$PBROKEN return
    
    !---------------------------------------------------------------------------
    ! The following situation concerns only the case for a multipole moment 
    ! with odd l and parity is conserved.
    !
    ! Using straightforward notation, this case corresponds to
    !
    !  U = ( U+ 0 )  V = ( V+ 0 )  Q = (0   Q+-)
    !      ( 0  U-)      ( 0  V-)      (Q-+ 0  )
    !
    ! And so
    !
    ! U^dagger Q V^*    = (  0              U+^dagger Q+- V-  )
    !                     ( U-^dagger Q-+ V+          0       )
    !
    !
    ! V^\dagger Q^t U^* = (  0              V+^dagger Q-+ U-  )
    !                     ( V-^dagger Q+- U+          0       )
    !
    if(mod(l,2).eq.1) then
      si = 0 ; sb = 0
      do B=1,8,4 ! Just an isospin loop ...
        ! Size of the positive parity block
        Tp = HFBlocks(B) + HFBLocks(B+1)
        ! Size of the negative parity block
        Tm= HFBlocks(B+2) + HFBLocks(B+3)
  
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Uplus : U in the positive parity block
        U = bogo(sb+        1:sb+  Tp     , sb+  Tp   +1:sb+2*Tp     )
        ! Vmin  : V in the negative parity block
        V = bogo(sb+2*Tp+Tm+1:sb+2*Tp+2*Tm, sb+2*Tp+Tm+1:sb+2*Tp+2*Tm)    
        
        ! U^+\dagger Q V^-
        Q20(si+1:si+Tp,si+Tp+1:si+Tp+Tm) = &
          &  matmul(matmul(transpose(U),Qsp(si+1:si+Tp,si+Tp+1:si+Tp+Tm)), V)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Umin : U in the negative parity block
        U = bogo(sb+2*Tp+1:sb+2*Tp+Tm      ,sb+2*Tp+Tm+1:sb+2*Tp+2*Tm)
        ! Vplus  : V in the positive parity block
        V = bogo(sb+ Tp+1:sb+2*Tp          ,sb+  Tp   +1:sb+2*Tp     )    

        tmp(si+1:si+Tp,si+Tp+1:si+Tp+Tm) = &
          &  matmul(matmul(transpose(V),Qsp(si+1:si+Tp,si+Tp+1:si+Tp+Tm)), U)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ! Putting both together
        ! Default formula
$NTR     Q20(si+1:si+Tp,si+Tp+1:si+Tp+Tm) = + Q20(si+1:si+Tp,si+Tp+1:si+Tp+Tm) &
$NTR     &                                  - tmp(si+1:si+Tp,si+Tp+1:si+Tp+Tm) 
!       ! Extra time-reversal sign when T is conserved
$TR     Q20(si+1:si+Tp,si+Tp+1:si+Tp+Tm) = - Q20(si+1:si+Tp,si+Tp+1:si+Tp+Tm) &
$TR     &                                  - tmp(si+1:si+Tp,si+Tp+1:si+Tp+Tm)
        
        ! And then we abuse hermeticity of these operators
        Q20(si+Tp+1:si+Tp+Tm, si+1:si+Tp) = &
        &                            transpose(Q20(si+1:si+Tp,si+Tp+1:si+Tp+Tm))
        
        si = si +   Tm  +   Tp
        sb = sb + 2*Tm  + 2*Tp
      enddo
    endif

  end function calc_Q20

  function Qlm_spme(l,m,Imaginary) result(me)
    !---------------------------------------------------------------------------
    ! Routine to calculate single-particle matrix elements of a given multipole
    ! moment in the (a) HARTREE-FOCK basis and (b) in units of b^(ell/2) with
    ! 1 b = 100 fm^2.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     l,m       : integers, 
    !     Imaginary : to calculate < i | Re Q_lm | j > or < i |Im Q_lm | j >
    !                 logical
    ! Output:
    !     me        : single-particle matrix elements of the multipole moment
    !                 in the canonical basis, real(nwt,nwt), in units of 
    !                 b^(ell/2).
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Notes:
    ! *) does not allow for the calculation of matrix elements with odd values 
    !    of m yet; those connect single-particle states with different 
    !    signature, for the traditional definition of spherical harmonics at 
    !    least.
    ! *) This routine shifts the multipole moments with respect to the center
    !    of mass of the nucleus, since it uses the meshx/y/z_shifted arrays to
    !    calculate the numerical values of the spherical harmonics on the mesh.
    !---------------------------------------------------------------------------

    integer, intent(in)        :: l, m
    logical, intent(in)        :: Imaginary
    real(KIND=dp)              :: me(nwt,nwt)
    
    real(KIND=dp), allocatable         :: SpherHarmMesh(:,:,:,:,:,:)
    real(KIND=dp), allocatable, target :: Qlm(:)
    real(KIND=dp), pointer             :: harm_3D(:,:,:)
    
    integer :: i, j, Bi, Bj, Ni, Nj, si, sj, im, k

    if(mod(m,2) .eq. 1) then
      print *, 'Qlm_spme does not allow for the calculation of multipole matrix elements with odd m YET.'
      stop
    endif
    if(Imaginary) then 
      print *, 'Qlm_spme does not allow for the calculation of imaginary multipole matrix elements YET.'
      stop
    endif
    if((m .lt. 0) .or. (m .gt. l)) then
      print *, 'Invalid value of m in Qlm_spe.'      
      stop
    endif
    
    ! Initialize to zero
    me = 0.0d0

    ! Generate the necessary spherical harmonic
    allocate(SpherHarmMesh(nx,ny,nz,0:l,0:l,2), Qlm(nx*ny*nz))
    call GenSphericalHarmonics(l,nx,ny,nz, &
    &                          meshx_shifted,meshy_shifted,meshz_shifted,    & 
    &                          SpherHarmMesh,quantisationaxis,secondaryaxis) 

    if(Imaginary) then
      im  = 2
    else
      im  = 1
    endif
    
    harm_3D(1:nx,1:ny,1:nz) => Qlm
    harm_3D                 = SpherHarmMesh(:,:,:,l,m,im)
    deallocate(SpherHarmMesh)

    !--------------------------------------------------------------------------- 
    ! Loop over the neutron single-particle states
    !---------------------------------------------------------------------------
    si = 0
    do Bi = 1, 4
      Ni =  HFBlocks(Bi) ; if(Ni.eq.0) cycle
      do Bj = 1, 4
        Nj = HFBlocks(Bj) ; if(Nj.eq.0) cycle

        if(bj.eq.1) then
          sj = 0
        else
          sj = sum(HFblocks(1:Bj-1))
        endif

        ! Parity selection rule: if l = even, then only single-particle states
        !                        of identical parity (and signature) contribute    
$PCONSERVED        if(mod(l,2) .eq. 0 .and.  Bi .ne. Bj)        cycle
        ! Parity selection rule: if l = odd, then only single-particle states
        !                        of different parity contribute
$PCONSERVED        if(mod(l,2) .eq. 1 .and. abs(Bi - Bj).ne.2 ) cycle
        ! (THIS is only needed if parity conserved of course)
$PBROKEN if( Bi .ne. Bj ) cycle

        do i=1,NI    
         do j=1,NJ
          ! me = Int d^3r Sum_sigma psi^*_i(r,sigma) psi_j(r,sigma) Qlm(r)
          me(si+i,sj+j) = 0
          do k=1,4
            me(si+i,sj+j) = me(si+i,sj+j) & 
            &                    + sum(HFpsi(:,k,si+i)*HFpsi(:,k,sj+j)*Qlm(:))
          enddo
          ! All these matrix elements are real if 
          ! (i)  we consider only real multipole moments
          ! (ii) time simplex is conserved
          ! 
          ! which means the matrix we store them in is symmetric
          me(si+i,sj+j) = me(si+i,sj+j) * dv 
          me(sj+j,si+i) = me(si+i,sj+j)
         enddo
        enddo
      enddo
      si = si + NI
    enddo

    !--------------------------------------------------------------------------- 
    ! Loop over the proton single-particle states
    !---------------------------------------------------------------------------
    si = nwn
    do Bi = 5, 8
      Ni =  HFBlocks(Bi) ; if(Ni.eq.0) cycle
      do Bj = 5, 8
        Nj = HFBlocks(Bj) ; if(Nj.eq.0) cycle
        sj = sum(HFblocks(1:Bj-1))

        ! Parity selection rule: if l = even, then only single-particle states
        !                        of identical parity (and signature) contribute    
$PCONSERVED        if(mod(l,2) .eq. 0 .and.  Bi .ne. Bj)        cycle
        ! Parity selection rule: if l = odd, then only single-particle states
        !                        of different parity contribute
$PCONSERVED        if(mod(l,2) .eq. 1 .and. abs(Bi - Bj).ne.2 ) cycle
        ! (THIS is only needed if parity conserved of course)
$PBROKEN if( Bi .ne. Bj ) cycle

        do i=1,NI    
         do j=1,NJ
          ! me = Int d^3r Sum_sigma psi^*_i(r,sigma) psi_j(r,sigma) Qlm(r)
          me(si+i,sj+j) = 0
          do k=1,4
            me(si+i,sj+j) = me(si+i,sj+j) & 
            &                    + sum(HFpsi(:,k,si+i)*HFpsi(:,k,sj+j)*Qlm(:))
          enddo
          ! All these matrix elements are real if 
          ! (i)  we consider only real multipole moments
          ! (ii) time simplex is conserved
          ! 
          ! which means the matrix we store them in is symmetric
          me(si+i,sj+j) = me(si+i,sj+j) * dv 
          me(sj+j,si+i) = me(si+i,sj+j)
         enddo
        enddo
      enddo
      si = si + NI
    enddo
    
    ! Rescale with the units of b^(ell/2) with 1 b = 100 fm^2.
    me = me/(100**(l/2.0))

  end function Qlm_spme

end module fission_MOI

! Code zoo
