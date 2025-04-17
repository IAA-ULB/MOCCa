module HFB_direct
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

  use wavefunctions
  use parameterization

  implicit none

  real(KIND=dp) :: blockJ = 1.0

contains 

  function ConstructConfiguration(Bogo, Eqp, blocks, blocktype, blockconf,     &
  &                               tag_overlaps, blocked_qp, partner_qp,        &
  &                               qp_overlap)  result(R)
    !---------------------------------------------------------------------------
    ! The generalized density matrix in a HFB calculation is given by
    !
    ! R = ( rho      kappa  )
    !     (-kappa^*  1-rho^*)
    !
    !   = W^T  ( f_1  0    ....  0    0     ....  0      )  W
    !          ( 0    f_2  ....  0    0     ....  0      )
    !          ( 0    0    ....  f_n  0     ....  0      )
    !          ( 0    0    ....  0    1-f_1 ....  0      )
    !          ( 0    0    ....  0    0     ....  0      )
    !          ( 0    0    ....  0    0     ....  1-f_n  )
    !
    ! where B is the Bogoliubov matrix and the matrix in the middle I call the 
    ! configuration matrix. 
    !
    ! For a simple HFB ground state, we have
    !
    !               f_i = 0    for all i=1,N
    !
    ! If we block any given quasiparticle k, we have
    ! 
    !               f_k = 1 
    !               f_i = 0    for all i not = k
    !
    ! If we are looking at calculations at finite temperature, we have 
    !
    !               f_i = (1+ exp(beta E^qp_i))^{-1}
    ! 
    ! where beta is the inverse temperature.         
    !
    ! This routine constructs the configuration matrix for a HFB calculation on 
    ! the basis of a great many relevant options, and (if needed) returns 
    ! indices of blocked qps, partner qps and the relevant overlaps. When this
    ! is done, we can construct the density and anomalous density matrices and 
    ! go on with life. 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input :
    !   Bogo      : Bogoliubov transformation obtained by diagonalisation of the 
    !               HFB Hamiltonian
    !   Eqp       : quasiparticle energies obtained by diagonalisation of the 
    !               HFB Hamiltonian
    !   blocks    : sizes of the symmetry blocks of the HFB hamiltonian
    !   blocktype : type of blocking to perform (1-6)
    !   blockconf : in which symmetry blocks to excite qps, for 
    !               blocktype = 3,4
    !   tag_overlaps: overlaps of the HFBasis wavefunctions with a tagging spwf
    !
    ! Output:
    !   blocked_qp: indices of the blocked quasiparticles. These are indexed
    !               from 1 to N, where N is the size of the single-particle 
    !               space, NOT from 1 to 2*N.
    !
    !   partner_qp: indices of the partner quasiparticles, i.e. the
    !               quasiparticles that are closely related to the blocked_qps.
    !
    !   qp_overlap: overlap between the time-reverse of blocked_qp(i) and 
    !               partner_qp(i)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Note: this functions assumes the ordering of the Bogoliubov transformation
    !       as used in the diagonalisation procedure of HFB_direct module. This
    !       is NOT the same as the ordering assumed in the rest of the program, 
    !       see the routine reorganise_matrices.
    !---------------------------------------------------------------------------

    integer, intent(in)          :: BlockType
    integer, intent(in)          :: blocks(4) 
    integer, intent(in)          :: blockconf(:)
    integer, allocatable         :: blocked_qp(:), partner_qp(:)
    real(KIND=dp), allocatable, intent(out) :: qp_overlap(:)
    real(KIND=dp), intent(in)    :: Eqp(:), Bogo(:,:), tag_overlaps(:)
    real(KIND=dp), allocatable   :: R(:)
       
    integer                      :: N, N2,B, sb, i, NB, j, k,qblock, ind, si, bi
    real(KIND=dp)                :: compare, occ, qpmin, blockoverlap, maxover,O
    integer                      :: toblock(4), qpb, indover
    logical                      :: found

    N = size(Eqp) 
    allocate(R(N)) ;  R = 0
    qpb = -1

    !---------------------------------------------------------------------------
    ! Construct the DEFAULT configuration, corresponding to all positive energy
    ! quasiparticles.
    sb = 0 ; si = 0

    do B=1,4
        N = blocks(B) ; if (N.eq. 0) cycle
        do i=1,N
            if(inversetemp .gt. 0.0_dp) then
                !---------------------------------------------------------------
                ! At finite temperature, things can get partially occupied and
                ! we are dealing with a statistical mixture.
                occ = exp(inversetemp * Eqp(si+i))
                occ = 1.0/(1 + occ)
                R(sb+N+i) = 1.0_dp - occ
                R(sb  +i) =          occ
            else
                !---------------------------------------------------------------
                !Completely empty or full, we want pure HFB states.
                R(sb+N+i) = 1.0_dp
                R(sb  +i) = 0.0_dp
            endif
        enddo
        si = si +   N
        sb = sb + 2*N
    enddo
    !---------------------------------------------------------------------------
    occ = 0
    select case(Blocktype)
    case(1,2,7)
        ! Full blocking
        occ = 1.0_dp
    case(3,4)
        ! EFA blocking
        occ = 0.5_dp
    case(5)
        ! Spherical-like blocking
        occ = 1/(2*blockJ +1)
    end select
    !---------------------------------------------------------------------------
    ! Modify this default configuration when needed.
    select case(Blocktype)
    case(0)
        !-----------------------------------------------------------------------
        ! No blocking asked for. 
    case(1,3,5)
$NTR    if(blocktype.eq.3) then
$NTR     call stp('Can not do EFA blocking when T is broken.')
$NTR    endif
        !-----------------------------------------------------------------------
        ! The user asked for a specific configuration that needs to be 
        ! identified. The array blockconf now contains the indices in the 
        ! HF-basis.
        NB = size(blockconf)
        if(.not.allocated(blocked_qp)) then
          allocate(blocked_qp(NB)) ; blocked_qp = 1
        endif
        if(.not.allocated(partner_qp)) then
          allocate(partner_qp(NB)) ; partner_qp = 1
        endif
        if(.not.allocated(qp_overlap)) then
          allocate(qp_overlap(NB)) ; qp_overlap = 1.0d0
        endif        

        blocked_qp = 0
        do j=1,NB
          compare = 0.0
          ind     = 0

          ! Check which block the requested index is in.
          call Identify(blockconf(j),blocks, bi, qblock)
          !-------------------------------------------------------------------
          ! Look for the column in the second half of the eigenvectors with
          ! the largest overlap with asked for state, but taking care not to 
          ! select the same qp twice
          sb = 0; si = 0
          do B=1,4
            N = blocks(B) ; if (N.eq.0) cycle
            if(B.eq.qblock) then
              do i=N+1,2*N
                ! removing previously selected qps from the comparison
                found = .false.
                do k=1,j
                  if(si+i-N .eq.blocked_qp(k)) found=.true.
                enddo              
                if(found) cycle
              
                if(Bogo(sb+bi,sb+i)**2 .gt. compare) then
                  compare = Bogo(sb+bi+N,sb+i)**2 
                  ind     = i
                endif
              enddo
            endif
            si = si +  N
            sb = sb +2*N
          enddo 
          !-------------------------------------------------------------------
          ! Change the occupation of this particular qp
          sb = 0 ; si =0 
          do B=1,4
              N = blocks(B)
              if(qblock.eq.B) then
                  R(sb+ind-N)       = occ
                  R(sb+ind)         = 1 - occ
                  ! Save which one we blocked
                  blocked_qp(j)     = si+ind-N
              endif
              sb = sb + 2*N
              si = si + N
          enddo    
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          ! Then we look for the closest thing to a time-reversal partner. 
          sb = 0 ; si =0 
          do B=1,4,2
              N = blocks(B)
              N2= blocks(B+1)
              if(qblock.eq.B .or. qblock.eq.B+1) then
                 call find_partner(N,N2,si,blocked_qp(j), &
                 &                 Bogo(sb+1:sb+2*N+2*N2,sb+1:sb+2*N+2*N2),  &
                 &                 partner_qp(j),qp_overlap(j))
              endif
              sb = sb + 2*N + 2*N2
              si = si +   N +   N2
          enddo    
        enddo
    case(2,4)
        !-----------------------------------------------------------------------
        ! The user asked for a the lowest configuration of a specific type.
        ! In this case, blockconf contains the number of qp excitations to  
        ! construct in every block.
        toblock = blockconf(1:4)

        if(blockconf(5).ne.0) then
          do i = 1, blockconf(5)
            qpmin = 10000000
            si = 0 ;  sb = 0
            do B=1,4,2
              N = blocks(B) ; if (N.eq.0) cycle
              N2= blocks(B+1) 
              if(Eqp(sb+N+toblock(B)+1) .lt. qpmin) then
                qpmin = Eqp(sb+N+toblock(B)+1)
                qpb   = B
              endif
              si = si +   N +   N2
              sb = sb + 2*N + 2*N2
            enddo
            toblock(qpb) = toblock(qpb) + 1

            !--------------------------------------------------------
            ! I'm not entirely sure why I coded this before....
            ! if(blocktype.eq.4) toblock(qpb+1) = toblock(qpb+1) +1 

          enddo
        endif

        NB = sum(toblock)
        if(allocated(blocked_qp)) then
          deallocate(blocked_qp)
          deallocate(partner_qp)
        endif

        if(.not.allocated(blocked_qp)) then
          allocate(blocked_qp(NB)) ; blocked_qp = 0
          allocate(partner_qp(NB)) ; partner_qp = 0
          allocate(qp_overlap(NB)) ; qp_overlap = 0.0d0
        endif        

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  For every block, we flip the required number of qps.  
        sb = 0 ; si = 0 ; ind = 0
        do B=1,4
          N = blocks(B) ; if (N.eq.0) cycle
          do j=1,toblock(B)
            ! The qps are ordered in energy from the diagonalization
            !  So we simply flip the first ones
            R(sb + N + j ) = 1 - occ
            R(sb     + j ) =     occ

            ! Saving the one we flipped
            ind = ind + 1
            blocked_qp(ind) = si + j
          enddo
          si = si +   N
          sb = sb + 2*N
        enddo
              
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Then we look for the closest thing to a time-reversal partner. 
        sb = 0 ; si = 0 ; ind = 0
        do B=1,4,2
           N = blocks(B)   ; if(N.eq.0) cycle    
           N2= blocks(B+1) 

           do j=1,toblock(B)
              call find_partner(N,N2,si,blocked_qp(j), &
                   &                 Bogo(sb+1:sb+2*N+2*N2,sb+1:sb+2*N+2*N2),  &
                   &                 partner_qp(j),qp_overlap(j))
           enddo
           si = si +   N +  N2 
           sb = sb + 2*N +2*N2         
        enddo
    case (7)
        !-----------------------------------------------------------------------
        ! We have a tagging state in memory, so we will look which quasiparticle
        ! has the largest overlap with said state.
        NB = sum(blockconf)
        if(allocated(blocked_qp)) then
          deallocate(blocked_qp)
          deallocate(partner_qp)
        endif

        if(.not.allocated(blocked_qp)) then
          allocate(blocked_qp(NB)) ; blocked_qp = 0
          allocate(partner_qp(NB)) ; partner_qp = 0
          allocate(qp_overlap(NB)) ; qp_overlap = 0.0d0
        endif

        B = 0
        do i=1,4
          if(blockconf(i) .eq. 1) B = i 
        enddo
        if(B.eq.0) return
        if(B.gt.1) then
          si =   sum(blocks(1:B-1))
          sb = 2*sum(blocks(1:B-1))
        else
          si = 0
          sb = 0
        endif
        N = blocks(B)
        N2 = N ! TO BE CHANGED

        ! ... and now calculate the overlap of the U and V parts of the
        ! quasiparticles with the tagging spwfs
        maxover = -10
        indover =   0

        ! First scan the U-matrix
        do i=1,N
          O = 0
          do j=1, N
            O = O + Bogo(sb    +j,sb+N+i) * tag_overlaps(si+j)
            ! Please read the note on organisation of the Bogo matrix in
            !  this routine before you start tinkering with this statement!
          enddo
          O = abs(O)
          if(O .gt. maxover) then
              maxover = O
              indover = i
          endif
        enddo

        ! ... and then the V-matrix
        do i=1, N
          O = 0
          do j=1, N
            O =  O + Bogo(sb+N+j,sb+N+i) * tag_overlaps(si+j)
            ! Please read the note on organisation of the Bogo matrix in
            !  this routine before you start tinkering with this statement!
          enddo
          O  = abs(O)
          if(O .gt. maxover) then
              maxover = O
              indover = i
          endif
        enddo

        ! Make the selection!
        R(sb +  N + indover ) = 1 - occ
        R(sb      + indover ) =     occ

        blocked_qp(1) = si + indover
        blockoverlap  = maxover

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Then we look for the closest thing to a time-reversal partner.
        sb = 0 ; si = 0 ; ind = 0
        do B=1,4,2
           N = blocks(B)   ; if(N.eq.0) cycle
           N2= blocks(B+1)

           if(blockconf(B) .eq. 1 .or. blockconf(B+1) .eq. 1) then
              call find_partner(N,N2,si,blocked_qp(1), &
                   &                 Bogo(sb+1:sb+2*N+2*N2,sb+1:sb+2*N+2*N2),  &
                   &                 partner_qp(1),qp_overlap(1))
           endif
           si = si +   N +  N2
           sb = sb + 2*N +2*N2
        enddo

    end select

  end function ConstructConfiguration
  
  subroutine find_partner(N,N2,si,qp,bogo,partner,overlap)
    !---------------------------------------------------------------------------
    ! Find the partner-qp (= almost time-reversal partner) within the right 
    ! symmetry block. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input: 
    !    N, N2 : size of the two symmetry blocks related by a linear,
    !            antihermitian symmetry. 
    !    si    : matrix indices of the HFB transformation
    !    qp    : index of the original blocked qp
    !    bogo  : part of the Bogoliubov transformation in the right subblocks.
    !            Size     : 2*N+2*N2
    !            Subblock : sb+1:sb+2*N+2*N2, sb+1:sb+2*N+2*N2  
    !
    ! Output:
    !    partner: index of the partner-qp
    !    overlap: overlap of the partner-qp with the time-reversed original qp
    !
    !---------------------------------------------------------------------------
    integer, intent(in)       :: N, N2, si, qp
    real(KIND=dp), intent(in) :: bogo(:,:)
    integer, intent(out)      :: partner
    real(KIND=dp), intent(out):: overlap

    real(KIND=dp), allocatable :: tr_qp(:)
    real(KIND=dp)              :: overl
    integer                    :: i, column, offset
    
    allocate(tr_qp(2*N+2*N2))

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Inside the direct-HFB routines, the Bogoliubov matrix is 
    ! ordered somewhat differently from the rest of the program. 
    ! Remember
    !
    !       (  V^*+  U+    0     0   )
    !  W =  (  U^*+  V+    0     0   )
    !       (  0     0     V^*-  U-  )
    !       (  0     0     U^*-  V-  )
    !   
    !                ^           ^
    !                |           | 
    !               (1)         (2)
    !
    ! because of the diagonalisation in subblocks. 
    !
    ! Hence, a "blocked" qp is located in (1), while its possible
    ! time-reversal partners are located in (2).
    !
    ! Note that this particular routine is not yet ready for 
    ! blocked quasiparticles with signature -i.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if(qp - si .lt. N) then
      ! our blocked qp is in the first of the two related blocks
      column = qp-si+N2
      offset = N+N2
  
      tr_qp(       1:  N     ) =   0
      tr_qp(  N   +1:  N+  N2) =   0
      tr_qp(  N+N2+1:2*N+  N2) = + Bogo(  1:N   ,column)
      tr_qp(2*N+N2+1:2*N+2*N2) = - Bogo(N+1:N+N2,column)  
                               ! Sign due to time-reversal
    else
      column = qp-si+N2+N
      offset = 0

      tr_qp(       1:  N     ) = + Bogo(  N+N2+1:2*N+  N2,column)
      tr_qp(  N   +1:  N+  N2) = - Bogo(2*N+N2+1:2*N+2*N2,column) 
      tr_qp(  N+N2+1:2*N+  N2) =   0
      tr_qp(2*N+N2+1:2*N+2*N2) =   0
                               ! Sign due to time-reversal
    endif
    
    overlap = -100000d0
    do i=1,N+N2
       overl = sum(tr_qp(:) * Bogo(1:2*N+2*N2,offset+i))

       if(abs(overl) .gt. overlap) then
         overlap = abs(overl)
         partner = si + i 
       endif
    enddo 
    return
  end subroutine find_partner

  subroutine FindFermi_secant(H, blocks, targetparticles, config, Bogo, Eqp,   & 
            &                 lambda, maxhfbiter, blocktype, blockconf,        &
            &                 tag_overlaps, blocked_qp, partner_qp, qp_overlap,&
            &                 ifail)
      !-------------------------------------------------------------------------
      ! Subroutine that diagonalizes the HFB hamiltonian (repeatedly) to find  
      ! the correct Fermi energy that fixes the average number of particles.
      ! The routine only solves this for one particular isospin.
      !
      ! Input
      !   H               : HFB hamiltonian, without Fermi energy
      !   blocks          : Sizes of the symmetry blocks that can be used to 
      !                     simplify the problem.
      !   targetparticles : Average number of particles to target. 
      !   lambda          : Initial guess for the Fermi energy
      !   maxhfbiter      : Maximum number of iterations to perform
      !
      ! Output
      !   config   : Configuration matrix of the final solution
      !   Eqp      : Quasiparticle energies of the final solution
      !   Bogo     : Bogoliubov transformation that diagonalizes H
      !   Lambda   : Final fermi energy
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in)    :: H(:,:), targetparticles, tag_overlaps(:)
      real(KIND=dp), intent(out)   :: config(:), Bogo(:,:), Eqp(:)
      real(KIND=dp), intent(inout) :: lambda
      integer, intent(in)          :: blocks(4), maxhfbiter, blocktype
      integer, intent(in)          :: blockconf(:)
      integer, intent(out)         :: ifail
      integer, allocatable         :: blocked_qp(:), partner_qp(:)
      real(KIND=dp), intent(out), allocatable :: qp_overlap(:)

      real(KIND=dp)                :: df, dn(2), particles
      integer                      :: iter

      ! Initialization
      df = 0 ; dn = 0.0
      do iter=1, maxHFBiter

        particles = & 
        &   diagbyblock(H,blocks, config, Bogo,Eqp,lambda, blocktype,blockconf,&
        &               tag_overlaps,blocked_qp, partner_qp, qp_overlap, ifail)
        ! Return if we do not want to readjust the Fermi energy
        if(MaxHFBiter.eq.1) return
        !-----------------------------------------------------------------------
        ! Readjust the Fermi energy based on the number of particles.
        ! We use the secant method.
        dn(2) = dn(1)
        dn(1) = particles - targetparticles

        if(abs(dn(1)).lt.pairing_prec) return

        if(iter.eq.1) then
          ! We try lambda + 0.1 for the first iteration
          lambda = lambda + 0.1
          df     =          0.1
        else
          df     = - dn(1) * df/(dn(1) - dn(2))
  
          if(abs(df).gt.1.0) df = 0.1 * df/abs(df)
          lambda = lambda + df
        endif
      enddo
  end subroutine FindFermi_secant

  function diagbyblock(H, blocks, config, Bogo,Eqp,lambda, blocktype,blockconf,&
  &                 tag_overlaps, blocked_qp, partner_qp, qp_overlap, ifail)   &
  &        result(particles)
      !-------------------------------------------------------------------------
      ! Routine that diagonalizes, block by block, a HFB Hamiltonian that is 
      ! passed in. It does the low-level work for all the high-level routines
      ! in this moodule. 
      !
      ! Input
      !   H        : HFB hamiltonian, without Fermi energy
      !   blocks   : Sizes of the symmetry blocks that can be used to simplify 
      !              the problem.
      !   lambda   : Fermi energy.
      !   blocktype | Options for the construction of the configuration matrix
      !   blockconf | see, contructConfiguration.
      !   tag_overlaps: overlaps between HFBasis and the tagging spwf
      ! 
      ! Output
      !   config    : Configuration matrix of the final solution
      !               ( see  constructconfiguration )
      !   Eqp       : Quasiparticle energies of the final solution
      !   Bogo      : Bogoliubov transformation that diagonalizes H
      !   Lambda    : Final fermi energy
      !   blocked_qp| outputs of the constructconfiguration routine regarding 
      !   partner_qp| blocking possibilities.
      !   qp_overlap|
      !   ifail     : 0 if succesful diagonalisation + fermi energy found
      !              1 otherwise
      !   particles : total number of particles for this input
      !
      !
      ! Note that this routine is designed to be called for FOUR symmetry 
      ! blocks, i.e. for ONE nucleon species. 
      !
      !-------------------------------------------------------------------------
      external :: DSYEV

      real(KIND=dp), intent(in)    :: H(:,:), tag_overlaps(:)
      real(KIND=dp), intent(out)   :: config(:), Bogo(:,:), Eqp(:)
      real(KIND=dp), intent(inout) :: lambda
      integer, intent(in)          :: blocks(4), blocktype
      integer, intent(in)          :: blockconf(:)
      integer, allocatable         :: blocked_qp(:), partner_qp(:)
      real(KIND=dp), allocatable, intent(out) :: qp_overlap(:)

      real(KIND=dp), allocatable   :: eigen(:), work(:), A(:,:)
      real(KIND=dp)                :: particles
      integer                      :: sb, si, N, B, i, ifail, lwork

      allocate(eigen(2*sum(blocks)))
      !-----------------------------------------------------------------------
      ! a) Diagonalization of the HFB Hamiltonian by block. 
      si = 0 ; sb = 0
      do B=1,4
        N = blocks(B) ; if(N .eq. 0) cycle

        ! Construct the blocks of H including the Fermi energy
        allocate(A(2*N,2*N)) ; A = 0
        
        A = H(sb+1:sb+2*N, sb+1:sb+2*N)
        do i=1,N
          A(i  ,i  ) = A(i  , i  ) - lambda
          A(i+N,i+N) = A(i+N, i+N) + lambda
        enddo
                        
        ! Diagonalize every block
        lwork = -1; allocate(work(1))
        call DSYEV( 'V', 'U', 2*N, A, 2*N, eigen(sb+1:sb+2*N),work,lwork,ifail)
        lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
        call DSYEV( 'V', 'U', 2*N, A, 2*N, eigen(sb+1:sb+2*N),work,lwork,ifail)
        deallocate(work)

        Bogo(sb+1:sb+2*N, sb+1:sb+2*N) = A

        if(ifail.ne.0) then
          print *, 'WARNING: diagon failed in subroutine DiagByBlock.'
          print *, '         Problematic block B = ', B
          deallocate(A, eigen)
          particles = 0.0
          return
        endif

        Eqp(sb+1:sb+2*N) = eigen(sb+1:sb+2*N)
        deallocate(A)
        ! Indices for the next block
        si = si +   N
        sb = sb + 2*N
      enddo
      !-----------------------------------------------------------------------
      ! b) We construct the configuration matrix that was asked for
      config = ConstructConfiguration(Bogo, Eqp, blocks, blocktype, blockconf, &
      &                       tag_overlaps, blocked_qp, partner_qp, qp_overlap)
      !-----------------------------------------------------------------------
      ! c) Count the total number of particles that we have.
      si = 0 ; sb = 0
      particles   = 0
      do B=1,4                
          N = Blocks(B) ;  if(N .eq. 0) cycle 
          !-------------------------------------------------------------------
          ! Calculate the number of particles in here  
          ! Sum_i rho_ii =  Sum_ii   U   f U^{\dagger} + V^{*}(1 - f)V^{T}
          !          sum_(ij>N) f_(j) V^*_ij V^T_ji = sum_ij f_(j) V^*_ij V_ij
          !        + sum_(ij<N) f_(j) U^*_ij U^T_ji = sum_ij f_(j) U^*_ij U_ij         
          do i=1,N
              particles = particles                                          &
              &         +   config(sb+N+i)*sum(bogo(sb+N+1:sb+2*N, sb+N+i)**2)            
          enddo
          do i=1,N
              particles = particles                                          &
              &         +   config(sb  +i)*sum(bogo(sb  +1:sb+  N, sb+N+i)**2)            
          enddo

          ! indices
          si = si +  N
          sb = sb +2*N
      enddo
      ! When Time-reversal is conserved, we need an extra factor of two
$TR   particles = 2 * particles                 
      deallocate(eigen)
  end function diagbyblock

  subroutine FindFermi_Brent(H, blocks, targetparticles, config, Bogo, Eqp,    & 
   &                         lambda, maxhfbiter, blocktype, blockconf,         &
   &                         tag_overlaps, blocked_qp, partner_qp, qp_overlap, &
   &                         ifail)
      !-------------------------------------------------------------------------
      ! Subroutine that diagonalizes the HFB hamiltonian (repeatedly) to find  
      ! the correct Fermi energy that fixes the average number of particles.
      ! The routine only solves this for one particular isospin.
      !
      ! This particular subroutine employs Brents method to fix the Fermi 
      ! energy. See
      ! 
      ! https://en.wikipedia.org/wiki/Brent%27s_method
      ! 
      ! which combines bisection, secant method and inverse quadratic 
      ! interpolation.The original source is probably
      ! R. P. Brent (1973), "Chapter 4: An Algorithm with Guaranteed Convergence
      ! for Finding a Zero of a Function", Algorithms for Minimization without
      ! Derivatives, Englewood Cliffs, NJ: Prentice-Hall,  
      !
      ! Input
      !   H        : HFB hamiltonian, without Fermi energy
      !   blocks   : Sizes of the symmetry blocks that can be used to simplify 
      !              the problem.
      !   lambda   : Initial guess for the Fermi energy
      !   maxhfbiter: Maximum number of iterations to perform
      !   targetparticles: Average number of particles to target. 
      !   tag_overlaps: overlaps between HFBasis and the tagging spwf
      !
      ! Output
      !   config   : Configuration matrix of the final solution
      !   Eqp      : Quasiparticle energies of the final solution
      !   Bogo     : Bogoliubov transformation that diagonalizes H
      !   Lambda   : Final fermi energy
      !-------------------------------------------------------------------------
      ! This routine is very heavily inspired/copy-pasted by the routines 
      ! implemented in MOCCa by M. Bender. 
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in)    :: H(:,:), targetparticles, tag_overlaps(:)
      real(KIND=dp), intent(out)   :: config(:), Bogo(:,:), Eqp(:)
      real(KIND=dp), intent(inout) :: lambda
      integer, intent(in)          :: blocks(4), maxhfbiter, blocktype
      integer, intent(in)          :: blockconf(:)
      integer, intent(out)         :: ifail  
      real(KIND=dp), allocatable, intent(out) :: qp_overlap(:)

      real(KIND=dp)                :: InitialBracket(2), FA, FB, N
      integer                      :: idir = 0 , idirsig = 1, FailCount
      logical                      :: Success
      integer, allocatable         :: blocked_qp(:), partner_qp(:)
 
      if ( targetparticles .lt. 0.1_dp ) then
        lambda = -10000.0_dp
      endif
      !-------------------------------------------------------------------------
      ! STEP 1: set up an initial bracket
      !-------------------------------------------------------------------------
      N = diagbyblock(H, blocks, config, Bogo,Eqp,lambda, blocktype,blockconf, &
      &               tag_overlaps, blocked_qp, partner_qp, qp_overlap, ifail)
      N = N - targetparticles
      ! Check if this guess for lambda is good enough
      if(abs(N).lt.pairing_prec) return

      ! Use present Fermi energy as starting point and check the direction
      ! where the zero of <N>-N0 can be expected.
      ! If <N>-N0 <  0, search at higher values.
      ! If <N>-N0 >= 0, search at lower  values.
      ! Initialize InitialBracket(it,1) = A, InitialBracket(it,2) = B with 
      ! present  Fermi energy.
      ! "dir" is the label of the InitialBracket(it,idir) that has to be moved,
      ! "idirsig" is the sign of steps needed to go into that direction.

      InitialBracket(:) = lambda
      if ( N .lt. 0.0_dp ) then 
        idir   =  2 ;  idirsig =  1
      else 
        idir   =  1 ;  idirsig = -1
      endif

      ! Try to find a boundary that brackets the Fermi energy in the direction 
      ! into which the Fermi energy has to be changed.
      FailCount = -1 ;  Success = .false.

      do while(.not. Success)
        FailCount = FailCount + 1

        ! update moving boundary and recalculate particle numbers at both.
        InitialBracket(idir) = &
        &                 InitialBracket(idir) + idirsig * 0.01_dp*(FailCount+1)

        FA = diagbyblock(H,blocks,config,Bogo,Eqp,InitialBracket(1),blocktype, &
        &     blockconf, tag_overlaps, blocked_qp, partner_qp, qp_overlap,ifail)
        FB = diagbyblock(H,blocks,config,Bogo,Eqp,InitialBracket(2),blocktype, &
        &     blockconf, tag_overlaps, blocked_qp, partner_qp, qp_overlap,ifail)
        FA = FA - targetparticles ; FB = FB - targetparticles

        ! check if N(epsilon_F) is a monotonically growing function.
        ! It should be, but who knows, pigs may fly ...
        if ( FB .lt. FA ) then 
          print '(" : Warning N(eps_F) decreases ")'
          print '(" A = ",f13.8," FA = ",f14.8," B = ",f13.8," FB = ",f14.8)', &
               & InitialBracket(1),FA+N, InitialBracket(2),FB+N
        endif

        ! diagnostic printing for convergence analysis (usually commented out)
!        print '(" Bracketing ",i4,1l2,(2(f13.8,es16.7)))',        &
!              & FailCount,Success,InitialBracket(1),FA, InitialBracket(2),FB          
        ! code failure (Fermi energy has changed by 30 MeV)
        if (Failcount .gt. 300) then
          print '(/," A = ", f13.8, " FA = ",1es12.4,              &
               &    " B = ", f13.8, " FB = ",1es12.4)',            &
               &     InitialBracket(1),FA,InitialBracket(2),FB 
          ifail = 1
          return
          !stop 'FindFermiBrent: Search for InitialBracket failed.'
        endif
        ! check if root is bracketed for isospin it after the update
        if( FA*FB .lt. 0.0_dp ) then 
            ! Correct Bracket found!
            Success = .true.
        endif
      enddo
      !-------------------------------------------------------------------------
      ! STEP 2: call the routine for the actual bisection
      !-------------------------------------------------------------------------
      call BrentBisection(lambda,N,InitialBracket(1), InitialBracket(2),FA,FB, &
      &                   maxHFBIter,H,blocks, targetparticles, config, Bogo,  & 
      &                   Eqp, blocktype, blockconf, tag_overlaps, blocked_qp, &
      &                   partner_qp, qp_overlap)
  
  end subroutine FindFermi_brent

  subroutine BrentBisection(lambda,particles, X1,X2,FX1,FX2,Depth, H, blocks,  & 
    &                       targetparticles, config, Bogo, Eqp, blocktype,     &
    &                       blockconf, tag_overlaps, blocked_qp, partner_qp,   &
    &                       qp_overlap)
    !---------------------------------------------------------------------------
    ! This routine searches for the Fermi energy
    ! by Brent's methods https://en.wikipedia.org/wiki/Brent%27s_method
    ! which combines bisection, secant method and inverse quadratic 
    ! interpolation. The original source is probably
    ! R. P. Brent (1973), "Chapter 4: An Algorithm with Guaranteed Convergence
    ! for Finding a Zero of a Function", Algorithms for Minimization without
    ! Derivatives, Englewood Cliffs, NJ: Prentice-Hall, 
    !---------------------------------------------------------------------------
    ! see pages 1188 - 1189 of http://apps.nrbook.com/fortran/index.html
    ! W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery,
    ! Numerical Recipes in Fortran in Fortran 90, Second Edition (1996).
    !---------------------------------------------------------------------------
    ! Input
    !   H        : HFB hamiltonian, without Fermi energy
    !   blocks   : Sizes of the symmetry blocks that can be used to simplify 
    !              the problem.
    !   maxhfbiter: Maximum number of iterations to perform
    !   targetparticles: Average number of particles to target. 
    !
    ! Output
    !   config   : Configuration matrix of the final solution
    !   Eqp      : Quasiparticle energies of the final solution
    !   Bogo     : Bogoliubov transformation that diagonalizes H
    !   Lambda   : Final fermi energy
    !   Particles: Final number of particles
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)    :: H(:,:), targetparticles, tag_overlaps(:)
    real(KIND=dp), intent(out)   :: config(:), Bogo(:,:), Eqp(:), lambda
    real(KIND=dp), intent(out)   :: particles
    integer, intent(in)          :: blocks(4), blocktype
    integer, intent(in)          :: blockconf(:)
    integer, intent(in)          :: Depth
    real(KIND=dp), intent(in)    :: X1 , X2, FX1 , FX2 
    integer, allocatable         :: blocked_qp(:), partner_qp(:)
    real(KIND=dp), allocatable, intent(out) :: qp_overlap(:)

    real(KIND=dp)                :: A , B, C , FA, FB , FC
    real(KIND=dp)                :: D , E, S , P  , Q , R 
    real(KIND=dp)                :: Num , Tol , XM 
    real(KIND=dp)                :: eps = 1.d-15
    integer                      :: FailCount, ifail
    logical                      :: Found

    A  = X1 ; B  = X2 
    FA = FX1; FB = FX2
    Found = .false.    
    if (A .eq. B) then 
      !-------------------------------------------------------------------------
      ! This signals that FA = FB is zero within the tolerance.
      ! Either near-converged HFB or HF case of completely broken-down pairing
      ! which also satisfies FA = FB = 0 within an interval. The 
      ! latter case cannot be handled by the algorithm below.
      !-------------------------------------------------------------------------
      Found = .true.
    endif

    C = B ; FC = FB 
    E = -1000000 ; D = -1000000
  
    FailCount = -1

    do while(.not.Found) 
      FailCount = FailCount + 1
      if ( ( FB .gt. 0.0_dp .and. FC .gt. 0.0_dp ) .or. & 
         & ( FB .lt. 0.0_dp .and. FC .lt. 0.0_dp ) )  then
        C  = A     ;  FC = FA
        D  = B - A ;  E  = D
      endif
      if ( abs(FC) .lt. abs(FB) ) then
        A  = B ;  FA = FB
        B  = C ;  FB = FC
        C  = A ;  FC = FA
      endif
      !-------------------------------------------------------------------------
      ! Convergence check
      ! Note (W.R.): I have tightened convergence a bit compared to the values
      !              in MOCCa by M.B. 
      !-------------------------------------------------------------------------
      Tol  = 2.0_dp * eps * abs(B) + 0.05_dp * Pairing_prec
      XM   = 0.5_dp * (C-B)
      !----------------------------------------------------------------
      ! Note: the tolerance is on the precision of the Fermi energy,
      ! NOT the nearness of the particle number to the targeted value.
      !----------------------------------------------------------------
      if ( abs(XM) .le. Tol .or. FB .eq. 0.0_dp ) then
        Lambda =  B
        Found = .true. 
        cycle
      endif
      if ( abs(E) .ge. Tol .and. abs(FA) .gt. abs(FB) ) then
        S = FB/FA
        if ( A .eq. C ) then
          P = 2.0_dp * XM * S
          Q = 1.0_dp - S
        else
          Q = FA/FC
          R = FB/FC
          P = S * (2.0_dp * XM * Q * (Q-R) & 
                  &    - (B-A)*(R-1.0_dp))
          Q = (Q-1.0_dp)*(R-1.0_dp)*(S-1.0_dp)
        endif
        if ( P .gt. 0.0_dp ) Q = -Q
        P = abs(P)
        if (2.0_dp * P .lt. min(3.0_dp*XM*Q - abs(Tol*Q),abs(E*Q))) then
          E = D
          D = P / Q
        else
          D = XM
          E = D 
        endif
      else
        D = XM
        E = D 
      endif
      A  = B 
      FA = FB
      B  = B + merge(D,sign(Tol,XM),abs(D) .gt. Tol)    
  
      !-------------------------------------------------------------------------
      ! B is present best guess for the fermi energy, FB the corresponding 
      ! particle number.
      !-------------------------------------------------------------------------
      Num = diagbyblock(H, blocks, config, Bogo,Eqp,B, blocktype,blockconf,    &
      &                 tag_overlaps, blocked_qp, partner_qp, qp_overlap, ifail)
      FB  = Num - targetparticles

      !-------------------------------------------------------------------------
      ! diagnostic printing for convergence analysis (usually commented out)
      !-------------------------------------------------------------------------
      ! NOTE: B is the the best guess for the zero of F. A has been the previous
      ! "closest" interval boundary that is not updated after Found
      ! is set to .true. The actual zero might therefore be outside the 
      ! interval [A,B]. If so, the true zero is typically closer to B than 
      ! A is to B.
      !-------------------------------------------------------------------------
      ! Note further: as A and B are swapped from time to time, B might be 
      ! smaller than A when printed here
      !-------------------------------------------------------------------------
!      print '(" BrentBisection ",i4,(1l2,2(f13.8,es16.7),f14.8))',    &
!           & FailCount, Found,A,FA,B,FB,Num
!     
      if ( FailCount .gt. Depth ) then
        print '(/," Warning: BrentBisection did not converge after ",i4," iterations")', & 
        &      FailCount
      endif
    enddo
    ! Output
    Lambda    = B ; particles = FB
    
   !---------------------------------------------------------------------------
   ! output
   !---------------------------------------------------------------------------
   if ( abs(FA) .lt. abs(FB) ) then
     FA = diagbyblock(H, blocks, config, Bogo,Eqp,A, blocktype,blockconf,      &
        &             tag_overlaps, blocked_qp, partner_qp, qp_overlap, ifail)
     FA = FA - targetparticles
     lambda = A
   else
     FB = diagbyblock(H, blocks, config, Bogo,Eqp,B, blocktype,blockconf,      &
        &             tag_overlaps, blocked_qp, partner_qp, qp_overlap, ifail)
     FB = FB - targetparticles
     lambda = B
   endif
  end subroutine BrentBisection

  subroutine Identify(i, blocks, bi, qblock)
    !---------------------------------------------------------------------------
    ! Given a state in the HF-basis (labelled by i), we find the corresponding 
    ! symmetry block (bi) and the corresponding index (bi) in said block.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !       i     :  label of the sp state in the HF-basis
    !       blocks:  size of the symmetry blocks
    ! Output:
    !       bi    : index (=row) in the matrix for this particular HF-state
    !       qblock: symmetry block to which this HF-state belongs
    !---------------------------------------------------------------------------
    integer, intent(in)  :: i, blocks(4)
    integer, intent(out) :: bi, qblock    
        
    integer :: sb, N, B

    sb = 0; bi = 0; qblock = 0
    do B=1,4
        N = Blocks(B) ; if (N.eq.0) cycle

        if( i .gt. sb .and. i.le.sb+N) then
            qblock = B
            bi     = i - sb
            return
        endif
        sb = sb + N
    enddo
  end subroutine Identify

  subroutine reorganise_matrices(Bogo,qpe,c)
    !---------------------------------------------------------------------------
    ! The diagonalisation of the HFB Hamiltonian is not performed with the 
    ! same block structure as is adopted in the rest of the program. This 
    ! routine takes the quantities obtained by diagonalisation and reorders
    ! them so that they can be used in the rest of the code.
    !
    !            Input                        Output
    !  H : 
    !  B :   Bogo transform 
    !---------------------------------------------------------------------------
    ! Naively writing down the HFB hamiltonian, we have 
    ! (if there is a linear, antihermitian conserved symmetry)  
    ! 
    !
    !       (  h+  0     0     d+- )
    !  H =  (  0   h-    d-+   0   )
    !       (  0  -d-+  -h+    0   )
    !       ( -d+- 0     0    -h-  ) 
    !
    !
    !  but we diagonalize in practice the two submatrices
    !
    !  H+ =  (  h+     d+-)     H- = ( h-     d-+ )
    !        (  -d+-  -h- )          ( d-+   -h+  )
    !
    !  as
    !         H = ( H+ 0 )
    !             ( 0  H-)
    !     
    ! Hence, the Bogoliubov transformation in memory is structured as
    !        
    !       (  V^*+  U+    0     0   )
    !  W =  (  U^*+  V+    0     0   )
    !       (  0     0     V^*-  U-  )
    !       (  0     0     U^*-  V-  )
    !   
    !  which we need to correct by moving things around to
    !
    !       (  V^*+   0     U+  0   )
    !  W =  (  0      V^*-  0   U-  )      (*)
    !       (  0      U^*-  0   V-  )
    !       (  U^*+   0     V+  0   )
    !
    ! Note that this reordering is necessary too for the
    !  (i)  QPenergies
    !  (ii) Configmatrix
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Almost(!) all of these matrices are ordered by increasing value of 
    ! quasi-particle energy. So, in the Bogoliubov transformation at (*),
    ! the columns correspond to the following ordering of qp energies   
    !
    !  -E_-N ... -E_-1 , -E_+N ... -E_+1, E_+1, ... E_+N, E_-1 .... E_-N
    !
    !  where E_+/-1 .... E_+/-N are increasing sequences, and I draw attention
    !  to the fact that the first half of quasi-particles is organised in 
    !  reverse order from the second half. 
    !
    ! The EXCEPTION to this rule is the configuration matrix configmatrix. That
    ! is constructed in blocks of increasing ABSOLUTE quasi-particle energy
    ! meaning that in that matrix the ordering is  
    !   
    ! -E_+1 ... -E_+N , -E_-1 ... -E_-N, E_+1, ... E_+N, E_-1 .... E_-N
    !
    ! meaning that the order has been reversed for the first half. This is 
    ! convenient for the construction of the pairing matrices, as then the  
    ! full generalized density matrix in quasi-particle representation 
    !
    ! looks like
    !          ( f_+1  0                     ......                    0  )
    !    C =   ( 0    ....                                                ) 
    !          (            f_-1                                          ) 
    !          (                  .....                                   )
    !          (                         1-f_+1                           ) 
    !          (                                  ......                  ) 
    !          (                                          1-f_-1          )
    !          ( 0                                                 .......)
    !
    !
    ! or more specifically:   C(i) = 1 - C(i+N+N2) as is often used here 
    ! and elsewhere in the code.
    !
    ! These warnings about matrix ordering of the first half of the Bogoliubov
    ! transformation is a little bit academical: the code has been constructed
    ! such that only the right half of the Bogoliubov transformation in (*), 
    ! their quasi-particle energies and the FULL configuration matrix enter
    ! the generalized density matrix, which is the only quantity affecting the
    ! rest of the program. 
    !---------------------------------------------------------------------------
    integer                      :: si, sb, B, N, N2
    real(KIND=dp), intent(inout) :: Bogo(:,:) , qpe(:), c(:)
    real(KIND=dp), allocatable   :: temp(:,:), tempqe(:), tempc(:) 
 
    !---------------------------------------------------------------------------
    ! Explicit allocation statements to satisfy high-level optimisations by 
    ! recent CRAY compilers.    
    allocate(temp(size(Bogo,1), size(Bogo,2)))
    allocate(tempqe(size(qpe)))
    allocate(tempc(size(c)))

    temp = Bogo ;  tempqe = qpe   ; tempc = c
    Bogo = 0    ; qpe     = 0.0d0 ; c     = 0.0

    si = 0 ; sb = 0
    do B=1,8,2
      N = HFBlocks_global(B) ; N2 = HFBlocks_global(B+1)

      !-------------------------------------------------------------------------
      ! Moving the first block

      ! First N eigenvalues of  H+ 
      ! (negative qp energies generally, but not always)
      Bogo(sb+       1:sb  +N   ,sb+1:sb+N) = &
      &                                        temp(sb+  1:sb+  N,sb+  1:sb+N)
      Bogo(sb+N+2*N2+1:sb+2*N+2*N2,sb+1:sb+N) = &
      &                                        temp(sb+N+1:sb+2*N,sb+  1:sb+N)

      ! Second set of N eigenvalues of H+
      Bogo(sb+       1:sb+  N     ,sb+N+  N2+1:sb+2*N+  N2) = &
      &                                        temp(sb+  1:sb+  N,sb+N+1:sb+2*N)
      Bogo(sb+N+2*N2+1:sb+2*N+2*N2,sb+N+  N2+1:sb+2*N+  N2) = &
      &                                        temp(sb+N+1:sb+2*N,sb+N+1:sb+2*N)
      !-------------------------------------------------------------------------
      ! Moving the second block (which doesn't exist if T is conserved)

      ! First N2 eigenvalues of H-
      ! (negative qp energies generally, but not always)
      Bogo(sb+N+1:sb+N+2*N2,sb+N+1:sb+N+N2) = &
      &                          temp(sb+2*N+1:sb+2*N+2*N2,sb+2*N+1:sb+2*N+N2)
      ! Second set of N2 eigenvalues of H-
      Bogo(sb+N+1:sb+N+2*N2,sb+2*N+N2+1:sb+2*N+2*N2) = &
      &                       temp(sb+2*N+1:sb+2*N+2*N2,sb+2*N+N2+1:sb+2*N+2*N2)

      !-------------------------------------------------------------------------
      ! First N eigenvalues of H+
      c(sb+     1:sb+N    )  = tempc (sb  +1:sb+  N)
      ! Second set of N eigenvalues of H+
      c(sb+N+N2+1:sb+2*N+N2) = tempc (sb+N+1:sb+2*N) 
      !-------------------------------------------------------------------------
      ! First set of N2 eigenvalues of H-
      c(sb+N+     1:sb+N+N2  )   = tempc (sb+2*N+1   :sb+2*N+  N2)
      ! Second set of N2 eigenvalues of H-
      c(sb+2*N+N2+1:sb+2*N+2*N2) = tempc (sb+2*N+N2+1:sb+2*N+2*N2)

      !-------------------------------------------------------------------------
      qpe  (sb+     1:sb+N)      = tempqe(sb    +1:sb+N   ) 
      qpe  (sb+N+   1:sb+N+N2)   = tempqe(sb+2*N+1:sb+2*N+N2) 

      qpe  (sb+  N+N2+1:sb+2*N+  N2) = tempqe(sb+  N   +1:sb+2*N   )
      qpe  (sb+2*N+N2+1:sb+2*N+2*N2) = tempqe(sb+2*N+N2+1:sb+2*N+2*N2)

      !-------------------------------------------------------------------------
      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
    enddo
  end subroutine reorganise_matrices

end module HFB_direct
