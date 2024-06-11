module wavefunctions
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
 !
 ! Module containing the single-particle wave-functions (spwfs for short) and 
 ! many routines to calculate their properties.
 !
 !==============================================================================
 ! Some notes on the current MPI implementation
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! All large arrays, meaning those with a dimension of at least
 !
 !    mesh-points * number of wavefunctions = mv * nwt
 !
 ! are spread between MPI ranks in a naive block-distribution. These are the 
 ! spwfs, their derivatives and momentum_updates.
 !
 ! All other quantities with smaller dimensions (such as the spenergies etc.)
 ! are broadcast to all MPI ranks. Their memory requirements are comparatively
 ! small and simply having access to them everywhere will hopefully avoid 
 ! future bugs.
 !==============================================================================
 ! TODO: 
 ! 1. The routines transfer_psi and its colleagues offer an interface that is
 !    too complicated. The argument send_rank can be figured out INSIDE the 
 !    routine; it is superfluous. The argument wave would be easier to use if
 !    it referred to the GLOBAL index of the spwf.
 !==============================================================================
 ! Hephaestos keywords:
 ! 
 !    N2  : $N2
 !    N3  : $N3
 !
 !
 !    Symmetry properties of the spwfs
 !
 !            |  c  |    X   Y   Z
 !    -----------------------------------------
 !    Block 1 |  1  |   $SX11  $SY11  $SZ11 
 !            |  2  |   $SX12  $SY12  $SZ12
 !            |  3  |   $SX13  $SY13  $SZ13
 !            |  4  |   $SX14  $SY14  $SZ14
 !    -----------------------------------------
 !    Block 2 |  1  |   $SX21  $SY21  $SZ21 
 !            |  2  |   $SX22  $SY22  $SZ22
 !            |  3  |   $SX23  $SY23  $SZ23
 !            |  4  |   $SX24  $SY24  $SZ24
 !    -----------------------------------------
 !    Block 3 |  1  |   $SX31  $SY41  $SZ21 
 !            |  2  |   $SX32  $SY42  $SZ22
 !            |  3  |   $SX33  $SY43  $SZ23
 !            |  4  |   $SX34  $SY44  $SZ24
 !    -----------------------------------------
 !    Block 4 |  1  |   $SX41  $SY41  $SZ21 
 !            |  2  |   $SX42  $SY42  $SZ22
 !            |  3  |   $SX43  $SY43  $SZ23
 !            |  4  |   $SX44  $SY44  $SZ24
 !    -----------------------------------------
 ! 
 !==============================================================================
 use derivatives
 use nil8
 use timing

 implicit none
 
 !------------------------------------------------------------------------------
 ! Array containing the spwfs and their derivatives: for ease of use in density
 ! and derivative calculations, these are stored with spatial (nx*ny*nz points)
 ! and spin indices (real/imaginary parts of spin up/down) separately.
 !
 ! These spwfs are called the Hartree-Fock basis throughout the code (hence
 ! the name HFBasis), but they are not guaranteed to be the basis that
 ! diagonalises the sphamiltonian. An extra unitary transformation might be
 ! required among them to obtain the physical HF basis.
 !
 ! Her-order derivative tensors are stored in lexicographical order to cut down
 !  on the number of indices and wasted computation.
 !            1    2    3    4    5    6    7    8    9    10
 ! 1st order: Dx   Dy   Dz
 ! 2nd order: Dxx  Dxy  Dxz  Dyy  Dyz  Dzz
 ! 3rd order: Dxxx Dxxy Dxxz Dxyy Dxyz Dxzz Dyyy Dyyz Dyzz Dzzz
 real(KIND=dp), allocatable, target ::      HFPsi(:,:,:)
 real(KIND=dp), allocatable, target ::   HFdPsi(:,:,:,:)!First order derivatives
 real(KIND=dp), allocatable, target ::  HFddPsi(:,:,:,:)!Second order derivatives
 real(KIND=dp), allocatable, target :: HFdddPsi(:,:,:,:)!Third order derivatives
 !------------------------------------------------------------------------------
 ! A copy of the HFPsi array, to be redistributed across MPI ranks in a
 ! 2D block-cyclic distribution with row and column blocking factors
 ! BLOCK_FACTOR_ROW and BLOCK_FACTOR_COLUN, respectively.
 real(KIND=dp), allocatable         :: HFPSI_2D(:,:)
 !------------------------------------------------------------------------------
 ! Array containing the values of the spwfs in the Canonical basis
 ! and their derivatives.
 real(KIND=dp), allocatable, target::     CANPsi(:,:,:)
 real(KIND=dp), allocatable, target::  CANdPsi(:,:,:,:)!First order derivatives
 real(KIND=dp), allocatable, target:: CANddPsi(:,:,:,:)!Second order derivatives
 real(KIND=dp), allocatable, target::CANdddPsi(:,:,:,:)!Third order derivatives
 !---------------------------------------------------------------------------
 ! Pointer to which basis is supposed to be used to calculate the densities
 ! Based on pairingtype
 !  (0) HF  => use the HF basis
 !  (1) BCS => use the HF basis
 !  (2) HFB => Use the canonical basis
 real(KIND=dp), pointer ::      DenPsi(:,:,:)
 real(KIND=dp), pointer ::   DendPsi(:,:,:,:)
 real(KIND=dp), pointer ::  DenddPsi(:,:,:,:)
 real(KIND=dp), pointer :: DendddPsi(:,:,:,:)
 !------------------------------------------------------------------------------
 ! Store the change in the spwfs from last iteration for momentum
 ! (Stored here such that they can be basis-transformed by other modules)
 real(KIND = dp), allocatable, target :: Momentum_Updates(:,:,:)  
 !------------------------------------------------------------------------------
 ! Number of the blocks with the same quantum numbers that divide up the 
 ! the single-particle wavefunctions.
 ! 
 ! Any of the possible symmetry combinations give rise to at most eight 
 ! different symmetry blocks.
 ! 
 !   2 for protons <-> neutrons
 !   2 for a hermitian, linear operator      ( parity      in EV8) 
 !   2 for an antihermitian, linear operator ( z-signature in EV8)
 ! 
 ! If symmetries are not conserved, we can simply eliminate blocks by setting 
 ! their size to zero.
 !
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! Explanation of the bookkeeping
 !   --------------------------
 ! HFBlocks       : number of spwfs in a given symmetry-block stored LOCALLY, 
 !                  i.e. on the current MPI rank
 ! HFBlocks_global: number of spwfs in a given symmetry-block stored GLOBALLY,
 !                  i.e. across all MPI ranks
 ! nwn, nwp       : TOTAL number of neutron/proton spwfs across all MPI ranks
 ! nwt            : TOTAL number of wavefunctions across all MPI ranks
 ! nwt_local      : LOCAL number of wavefunctions on the current MPI rank
 ! rank_map       : identifies the MPI rank that holds a given spwf within 
 !                  the global communicator
 ! spwf_map       : identifies the index of a (locally stored) spwf in the 
 !                  TOTAL calculation. I.E. this maps
 !                       spwf  1,  2, 3, ....,  nwt_local
 !                             |   |  |          |
 !                       spwf  X,  Y, Z, ....., AAA , AAA+1, ..... ,nwt
 ! spwf_inverse   : the LOCAL index of a given spwf on its MPI rank.
 !                     spwf   1, 2, 3,...., X,  X+1, ......, nwt 
 !                            |  |  |       |   |
 !                            1  2  3, ..., X,  1 , 2, .....
 !                  this will have many repeating entries; starting from a
 !                  global index X, one needs both (!)
 !                       rank_map(X) and spwf_inverse(X)
 !                  to uniquely find its location in memory
 !------------------------------------------------------------------------------
 integer, parameter   :: Blocks                  = 8  ! This can always be fixed
 integer              :: HFBlocks(Blocks)        = 0
 integer              :: HFBlocks_global(Blocks) = 0
 integer              :: nwn = 6, nwp = 6, nwt = 12, nwt_local =12
 integer, allocatable :: rank_map(:), spwf_map(:), spwf_inverse(:)
 ! The number of MPI ranks assigned to each symmetry block
 integer              :: ranks_per_block(Blocks)
 !------------------------------------------------------------------------------
 ! Properties of the single-particle wave-functions with regard to reflections
 ! of the axes. Note that these are properties of the LOCALLY stored spwfs, 
 ! i.e. they are allocated as sx(4,nwt_loc)
 integer, allocatable :: sx(:,:), sy(:,:), sz(:,:)
 ! Except for this one particular spwf: the one we use to estimate the maximal
 ! s.p. energy available on the grid in evolution.f90. When doing an MPI 
 ! calculation, this needs to be implemented separately such that all ranks
 ! use the same set symmetries to propagate it.
 integer              :: sx_max(4), sy_max(4), sz_max(4)
 !------------------------------------------------------------------------------
 ! Single-particle energies, 
 ! Either:
 ! (i)  diagonal elements of the single-particle hamiltonian: 
 !      \langle psi_i | h | psi_i \rangle
 ! (ii) eigenvalues of the single-particle hamiltonian when restricted to the
 !      subspace being iterated
 real(KIND=dp), allocatable :: spenergies(:)
 ! Complete single-particle hamiltonian 
 real(KIND=dp), allocatable :: sphamil(:,:)
 ! Dispersions of the spwfs with respect to h
 real(KIND=dp), allocatable :: dispersions(:)
 ! expectation values of the single-particle hamiltonian in the canonical basis
 real(KIND=dp), allocatable :: canenergies(:)
 !------------------------------------------------------------------------------
 ! Single-particle expectation values of Parity in the HF basis and in the 
 ! canonical basis
 real(KIND=dp), allocatable :: P_hf(:), P_can(:)
 !------------------------------------------------------------------------------
 ! Flag indicating whether or not to print advanced properties of the spwfs
 ! DURING the iterations. Their properties are calculated and printed for the 
 ! first and final iteration, but it can save quite some CPU time if they are 
 ! not calculated during the iterations in some conditions. Since this 
 ! information is not useful in the vast majority of cases, so we set this to 
 ! false by default. 
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! Advanced properties:
 !   - expectation values of different angular momentum operators
 !   - r^2 
 logical :: print_adv_spwf_properties = .false.
 !------------------------------------------------------------------------------
 ! Angular momentum properties of the spwfs in
 !  (i)   the ordinary basis, i.e. the spwfs in storage: spwf_[...]
 !  (ii)  the Hartree-Fock basis                       :   HF_[...]
 !  (iii) the canonical basis                          :  can_[...]
 ! "Ordinary" Jx, Jy, Jz 
 real(KIND=dp), allocatable :: spwf_J(:,:,:), hf_J(:,:), can_J(:,:)
 ! Squared   Jx^2, Jy^2, Jz^2
 real(KIND=dp), allocatable :: spwf_J2(:,:,:), HF_J2(:,:), can_J2(:,:)
 ! With an extra time-reversal operator JxT, JyT, JzT
 ! Both real and imaginary parts
 real(KIND=dp), allocatable :: spwf_JTR(:,:,:), HF_JTR(:,:), can_JTR(:,:)
 real(KIND=dp), allocatable :: spwf_JTI(:,:,:), HF_JTI(:,:), can_JTI(:,:)
 ! Total angular momentum "quantum number", i.e. the number J such that 
 !  J (J+1) = J^2_x +  J^2_y + J^2_z
 real(KIND=dp), allocatable :: spwf_JJ(:), HF_JJ(:), can_JJ(:)
 ! 
 ! Spin properties, i.e. < S_mu >
 real(KIND=dp), allocatable :: spwf_spin(:,:,:), hf_spin(:,:), can_spin(:,:)
 ! And values with a T-operator mixed in
 real(KIND=dp), allocatable :: spwf_STR(:,:,:) , hf_STR(:,:) , can_STR(:,:)
 real(KIND=dp), allocatable :: spwf_STI(:,:,:) , hf_STI(:,:) , can_STI(:,:)
 !
 ! Remark: when we diagonalise the sp hamiltonian explicitly, we are happy with 
 ! calculating only the diagonal matrix elements. When not diagonalising the 
 ! sphamiltonian explicitly, we are in need of the full matrices if we want to
 ! print information in the actual HF basis.  
 !------------------------------------------------------------------------------
 !------------------------------------------------------------------------------
 ! Oscillator frequencies to use for the initialization with a Nilsson  
 ! hamiltonian.
 real(KIND=dp) :: osc_freq(3) = (/ 0.2125, 0.2125, 0.175 /)
 !------------------------------------------------------------------------------
 ! Tell Tantalus to either 
 !  (i)  diagonalise the sp hamiltonian the ordinary way, i.e. using an
 !       iterative scheme
 !  (ii) to stop caring about the diagonalisation of the sphamiltonian
 !       and simply care about the space spanned by the spwfs.
 logical                    :: diagsphamil = .false.
 real(KIND=dp), allocatable :: HFtransfo(:,:)
 !------------------------------------------------------------------------------
 ! Use (or not) the more efficient implementation of the two-basis method
 logical :: efficientHFB = .false.
 !------------------------------------------------------------------------------
 ! The full matrix elements of <r^2> in the HF and canonical basis
 real(KIND=dp), allocatable, target :: spwf_r2_hf(:,:), spwf_r2_can(:,:)
 ! Important note: in many types of calculations, only the diagonal elements of
 ! this matrices will be calculated.
 !------------------------------------------------------------------------------
 ! Procedure to call to orthonormalize the s.p. wavefunctions in HFPSI.
 ! The code offers several strategies, hence the need for a procedure pointer.
 ! This pointer is assigned based on the input parameter 'ortho_strategy' of the
 ! evolution module.
 !
 ! 1. GramSchmidt :  (modified) Gram-Schmidt process
 ! 2. Loewdin     :  Loewdin (symmetric) orthonormalisation [NOT IMPLEMENTED]
 ! TODO: CORRECT THIS DOCUMENTATION
 !
 ! Note: this input parameter is not case-sensitive.
 procedure(GramSchmidt), pointer :: Orthonormalize
 !------------------------------------------------------------------------------
 ! Procedure used to initialise a bunch of spwfs. 
 ! Currently available:
 !     nilsson     : lowest-energy states of a simple Nilsson hamiltonian
 !                   used for finite nuclei.
 !     randomspwfs : random values, used for pasta calculations.
 ! The functioning of these routines is pretty particular, so please have 
 ! a look at their documentation. 
 ! The code currently offers no option to change the assignment of this pointer
 ! at runtime. 
 !------------------------------------------------------------------------------
#if(PASTA == 0)
 procedure(nilsson), pointer :: initialise_wavefunctions => nilsson
#else
 procedure(nilsson), pointer :: initialise_wavefunctions => randomspwfs
#endif
contains 

  subroutine ReadWFdata(file_number)
    !---------------------------------------------------------------------------
    ! Read the number of single-particle neutron and proton wave-functions.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   file_number : optional integer. If present, read from (open) channel
    !                 with this number. If absent, read from STDIN.
    !---------------------------------------------------------------------------

    integer(dp), intent(in), optional   :: file_number   
#if(USE_MPI>0)
    integer                             :: mpi_err
#endif

    namelist /wfs/ nwn, nwp, osc_freq, print_adv_spwf_properties, &
    &              max_spwf_per_rank, max_drop_ranks

    ! Only the first MPI rank reads input
    if(MPI_rank .eq. 0) then
      if(present(file_number)) then
        read(unit=file_number, nml = wfs)
      else
        read(unit=*, nml = wfs)
      endif
    endif
    
#if(USE_MPI > 0)
    ! Broadcasting all information
    call MPI_Bcast(nwn           , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(nwp           , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(max_drop_ranks, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(osc_freq      , 3, MPI_REAL8  , 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(print_adv_spwf_properties, 1, MPI_LOGICAL  , 0,             &
    &                                                   MPI_COMM_WORLD, mpi_err)

    call MPI_BCAST(max_spwf_per_rank ,1,MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
#endif
    ! Bookkeeping for all MPI ranks
    nwt = nwn + nwp
  end subroutine ReadWFdata

  subroutine loadbalance(blocks_global,balancing,blocks_local,spwf_map,        &
  &                                                       rank_map,spwf_inverse)
    !---------------------------------------------------------------------------
    ! Balance the loading of large arrays across MPI ranks in a 1D fashion.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   blocks_global : integer (8)
    !                   TOTAL size of all symmetry blocks, i.e. the total number
    !                   of spwfs in each block.
    !
    !   balancing     : integer.
    !                   0 => load balance purely on a spwf-by-spwf basis
    !                   1 => balance the load HFBLock-wise, i.e. assign 
    !                        spwfs per symmetry block. The number of MPI ranks
    !                        should be a multiple of the actually non-zero
    !                        symmetry blocks.
    ! Output:
    !   blocks_local  : integer(8)
    !                   LOCAL size of the symmetry blocks, i.e. the total number
    !                   of spwfs in each block FOR THIS MPI RANK.
    !   spwf_map      : integer(:)
    !                   mapping of the spwfs on this MPI rank to the whole
    !                   calculation
    !   rank_map      : integer(nwt)
    !                   mapping of the ranks for each spwf, i.e. spwf X is 
    !                   stored on MPI_RANK rank_map(X). 
    !   spwf_inverse  : integer(nwt)
    !                   index of the wavefunction defined by a global index on
    !                   the LOCAL MPI rank
    !---------------------------------------------------------------------------
    integer, intent(in)  :: balancing
    integer, intent(in)  :: blocks_global(blocks)
    integer, intent(out) :: blocks_local(blocks)
    integer, intent(out), allocatable :: spwf_map(:),rank_map(:),spwf_inverse(:)

    integer              :: B, activeblocks, blocks_per_rank, drop
    integer              :: block_count, i, offset, Nspwf, already_assigned
#if(USE_MPI>0)
    integer, external    :: numroc
    integer              :: mpi_err, dims(2), k, j, info, xsize
    integer              :: loc_psi, Bmax(1), local_ind,N
    integer, allocatable :: map_1D(:,:), map_2d(:,:),team(:), local_count(:)
    real(KIND=dp)        :: remaining, frac
#endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! The total number of spwfs to load-balance across different MPI ranks
    ! Note: not simply set to nwt to keep some flexibility...
    Nspwf  = sum(blocks_global)
    allocate(rank_map(Nspwf), spwf_inverse(Nspwf))
    rank_map     = 0 ; spwf_inverse = 0 ; blocks_local = 0; ranks_per_block = 0
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Required for printing assignments
    if(allocated(MPI_2D_COORDINATES))    deallocate(MPI_2D_COORDINATES)
    allocate(MPI_BLOCK_ASSIGNMENTS(NPROCS)); mpi_BLOCK_ASSIGNMENTS = 0
    allocate(MPI_2D_COORDINATES(NPROCS,2)) ; mpi_2D_coordinates    = 0
    ! Count the number of active symmetry blocks (blocks with non-zero spwfs)
    activeblocks = 0
    do B=1,8
      if(blocks_global(B) .ne. 0) then
        activeblocks = activeblocks + 1
        ! ensure that every active block gets
        !  (1) sufficient processes such that no process should go above max_spwf_per_rank
        !  (2) at least one attributed process
        !
        ! Important note: the code does not strictly enforce max_spwf_per_rank because
        ! of the rounding to nearest integer below; max_spwf_per_rank should be understood
        ! more as a rough guideline.
        ranks_per_block(B) = max(1,ceiling(blocks_global(B)/(1.0d0*max_spwf_per_rank)))
      endif
    enddo
    already_assigned = sum(ranks_per_block)

    select case(balancing)
    case (0)
      !-------------------------------------------------------------------------
      ! Full on load-balancing: each symmetry block gets assigned a bunch of
      ! MPI ranks, which divide equally the number of spwfs among them.
      !
      ! Todo: figure out if we can do more clever things by estimating the work
      !       to be done as a function of the number of spwfs in a given block
      !-------------------------------------------------------------------------

#if(USE_MPI==0)
      call stp('Balancing_strategy = 0 is not compatible with sequential calculations.')
#else
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Assign workload quadratically
      do B=1,8
          frac      = BLOCKS_GLOBAL(B)**2 / (1.0d0*sum(BLOCKS_GLOBAL**2))
          remaining = NPROCS - already_assigned
          ranks_per_block(B) = ranks_per_block(B) + NINT(frac * remaining)
      enddo
      ! This weighting might end up with a total number of ranks that is somewhat
      ! less than NPROCS, since we are dealing with the integer division. I solve
      ! this in the ad-hoc way of simply 
      !    1) adding the missing number of procs to the block that is largest.
      ! or 2) subtracting the superflous nprocs from the block with most nprocs
      if(Nprocs .gt. sum(ranks_per_block)) then
        Bmax = maxloc(Blocks_global)
        ranks_per_block(Bmax(1)) = ranks_per_block(Bmax(1)) + NPROCS - sum(ranks_per_block)
      elseif(Nprocs .lt. sum(ranks_per_block)) then
        Bmax = maxloc(ranks_per_block)
        ranks_per_block(Bmax(1)) = ranks_per_block(Bmax(1)) + NPROCS - sum(ranks_per_block)
        if(ranks_per_block(Bmax(1)) .le. 0) then 
          call stp('Invalid load balancing detected!')
        endif
      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Create a new communicator dealing with each symmetry block
      MPI_SYM_BLOCK = 1
      do B=1,8
        N = blocks_global(B); if(N.eq.0) cycle
        if(MPI_RANK +1.le. sum(ranks_per_block(1:B))) then
          if(MPI_RANK+1 .gt. sum(ranks_per_block(1:B-1))) then
            MPI_SYM_BLOCK = B
            ! MPI_SYM_BLOCK is the symmetry block this particular process
            ! is assigned to.
          endif
        endif
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Split the MPI_COMM_WORLD communicator
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, MPI_SYM_BLOCK, MPI_RANK, MPI_COMM_BLOCK, MPI_ERR)
      ! MPI_BLOCK_RANK is the MPI RANK of this particular process within its
      !                assigned symmetry block.
      call MPI_COMM_RANK(MPI_COMM_BLOCK, MPI_BLOCK_RANK, MPI_ERR)
      ! MPI_BLOCK_NPROCS is the total number of MPI ranks assigned to this
      !                symmetry block.
      call MPI_COMM_SIZE(MPI_COMM_BLOCK, MPI_BLOCK_NPROCS, MPI_ERR)
      ! MPI_BLOCK_SIZE is the number of spfs in this symmetry block;
      !   redundant information of course, but nice to have
      MPI_BLOCK_SIZE = blocks_global(MPI_SYM_BLOCK)

      call MPI_ALLGATHER(MPI_SYM_BLOCK,1,MPI_INT,MPI_BLOCK_ASSIGNMENTS,1,MPI_INT,MPI_COMM_WORLD,mpi_err)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Constructing the BLACS 1D and 2D layouts
      ! NOTE: this cannot be accomplished with the blacs_gridinit subroutine 
      !       because it does not support multigridding, i.e. grids that are 
      !       split by symmetry blocks such as we attempt here.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Default context containing ALL processes
      call blacs_get( 0, 0, blacs_cntxt)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! 1D context
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      blacs_cntxt_1D = blacs_cntxt
      allocate(map_1D(1,MPI_BLOCK_NPROCS))
      allocate(team(MPI_BLOCK_NPROCS))
      ! every member of the MPI_COMM_BLOCK sends its reank into the array team
      call MPI_ALLGATHER(MPI_RANK,1,MPI_INT,team,1,MPI_INT,MPI_COMM_BLOCK,mpi_err)
      map_1D(1,:) = team
      CALL BLACS_GRIDMAP (blacs_cntxt_1D, team , 1,  1, MPI_BLOCK_NPROCS)
      CALL BLACS_GRIDINFO(blacs_cntxt_1D,NROW_1D,NCOL_1D,MYROW_1D,MYCOL_1D)
      ! Blocking factor for the 1D distribution is engineered such that all 
      ! processes get one single contiguous block of memory. This simplifies
      ! the bookkeeping down below....
      BLOCK_FACTOR_1D = ceiling(MPI_BLOCK_SIZE/(1.0d0*MPI_BLOCK_NPROCS))
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! 2D context
      call blacs_get(0, 0, blacs_cntxt_2D)
      ! We try to determine an optimal 2D layout by repeatedly querying MPI_DIMS_CREATE
      ! until it does not return an error. For this to work, we first need to
      ! change the MPI error handling...
      CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN,mpi_err)

      mpi_err = 1
      if(MPI_BLOCK_NPROCS .gt. 5) then
        ! First try to get an (almost) square grid by setting dims(1) ~ sqrt(MPI_BLOCK_NPROCS)
        dims(1) = NINT(sqrt(1.0d0*MPI_BLOCK_NPROCS))
        dims(2) = 0
        call MPI_DIMS_CREATE(MPI_BLOCK_NPROCS,2, dims, mpi_err)
        drop = 0
        do while((MPI_ERR .ne. 0) .and. (drop+1 .lt. max_drop_ranks))
          drop = drop + 1
          call MPI_DIMS_CREATE(MPI_BLOCK_NPROCS-drop,2, dims, mpi_err)
          ! note: this will always succeed if MPI_BLOCKS_NPROCS-drop = 1,
          ! hence no need to test if drop >= MPI_BLOCK_NPROCS
        enddo
        if(mpi_err .ne. 0) then
          ! Explicitly set dims to 0 for just the default answer of MPI_DIMS_CREATE
          dims = 0
          call MPI_DIMS_CREATE(MPI_BLOCK_NPROCS,2, dims, mpi_err)
        endif

        ! Reset error handling to be fatal
        CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL,mpi_err)
      else
        dims = 0
        call MPI_DIMS_CREATE(MPI_BLOCK_NPROCS,2, dims, mpi_err)
      endif
      ! Size of the 2D team has been determined
      MPI_BLOCK_NPROCS_2D = dims(1) * dims(2)

      allocate(map_2D(dims(1),dims(2))); map_2D = 0
      K=  sum(ranks_per_block(1:MPI_SYM_BLOCK-1))
      do i=1,dims(1)
       do j=1,dims(2)
        map_2D(i,j) = K
        K = K+1
       enddo
      enddo
      CALL BLACS_GRIDMAP (blacs_cntxt_2D, map_2D, dims(1),dims(1), dims(2))
      CALL BLACS_GRIDINFO(blacs_cntxt_2D,NROW_2D,NCOL_2D,MYROW_2D,MYCOL_2D)
      call MPI_ALLGATHER(MYROW_2D,1,MPI_INT,MPI_2D_COORDINATES(:,1),1,MPI_INT,MPI_COMM_WORLD,mpi_err)
      call MPI_ALLGATHER(MYCOL_2D,1,MPI_INT,MPI_2D_COORDINATES(:,2),1,MPI_INT,MPI_COMM_WORLD,mpi_err)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! From the BLACS context, we now construct SCALAPACK descriptors
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! The 1D distribution is a column-cyclic one: each MPI rank gets complete
      ! columns of the spwf-matrix.
      CALL DESCINIT(desc_psi_1D,4*mv,MPI_BLOCK_SIZE,4*mv,                      &
      &             BLOCK_FACTOR_1D,0,0, blacs_cntxt_1d,4*mv,info)
      if(info.ne.0) then 
        call stp('Problem with DESCINIT call for psi_1D.')
      endif
      ! Scalapack knows exactly what amount of spwfs that are stored!
      loc_psi = NUMROC(MPI_BLOCK_SIZE,BLOCK_FACTOR_1D, MYCOL_1D,0, NCOL_1D)
      !  ------> this is a shorthand for use in the rest of this routine below
      blocks_local(MPI_SYM_BLOCK) = loc_psi
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! The 2D distribution of spwfs is a block-cyclic one with blocking factors
      ! decided by the user
      xsize = NUMROC(4*mv, block_factor_row, MYROW_2D,0, NROW_2D)
      if(xsize .le. 0) xsize = 1 ! things will get allocated, but this process 
                                 ! should not be participating in calculations
      if(MYROW_2D .ne. -1) then
        CALL DESCINIT(desc_psi_2D,4*mv,MPI_BLOCK_SIZE,    &
        &             block_factor_row, block_factor_col, &
        &             0,0,blacs_cntxt_2d, xsize ,info)
        if(info.ne.0) then
          call stp('Problem with DESCINIT call for psi_2D.')
        endif
      else
        ! Indicate that this particular process is not part of this context
        desc_psi_2d = -1
      endif
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! ... and similar for the descriptor of matrices in spwf x spwf space
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      xsize = NUMROC(MPI_BLOCK_SIZE, block_factor_row, MYROW_2D,0, NROW_2D)
      if(xsize .le. 0) xsize = 1 ! things will get allocated, but this process 
                                 ! should not be participating in calculations

      if(MYROW_2D .ne. -1) then
        CALL DESCINIT(desc_mat_2D,MPI_BLOCK_SIZE,MPI_BLOCK_SIZE,    &
        &             block_factor_row, block_factor_col, &
        &             0,0,blacs_cntxt_2d, xsize ,info)
        if(info.ne.0) then
          call stp('Problem with DESCINIT call for mat_2D.')
        endif
      else
        ! Indicate that this particular process is not part of this context
        desc_mat_2d = -1
      endif
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Constructing the bookkeeping on each process
      ! For this, each process needs to know the amount of spwfs each of his
      ! partners is storing, hence the MPI_ALLREDUCE
      allocate(spwf_map(loc_psi), local_count(0:MPI_BLOCK_NPROCS-1))
      local_count = 0
      local_count(MPI_BLOCK_RANK) = loc_psi
      call MPI_ALLREDUCE(MPI_IN_PLACE, local_count, MPI_BLOCK_NPROCS,     & 
      &                  MPI_INTEGER, MPI_SUM, MPI_COMM_BLOCK, mpi_err)

      ! The bookkeeping below works because we give all ranks one contiguous
      ! sets of spwfs because of our specific choice for the 1D blocking factor
      spwf_map  = 0
      local_ind = 0
      N      = MPI_BLOCK_SIZE
      offset = sum(blocks_global(1:MPI_SYM_BLOCK-1))
      do i=1,N
        if( i .gt. sum(local_count(0:MPI_BLOCK_RANK-1)) ) then
          if(MPI_BLOCK_RANK .ne. MPI_BLOCK_NPROCS) then
             if (i .le. sum(local_count(0:MPI_BLOCK_RANK))) then
                local_ind                = local_ind + 1
                spwf_map(local_ind)      = i + offset
                spwf_inverse(i + offset) = local_ind
                rank_map(i + offset)     = MPI_RANK
             endif
          else
            local_ind                    = local_ind + 1
            spwf_map(local_ind)          = i + offset
            spwf_inverse(i + offset)     = local_ind
            rank_map(i + offset)         = MPI_RANK
          endif
        endif
      enddo
      !-------------------------------------------------------------------------
#endif
    case (1)
      !-------------------------------------------------------------------------
      ! Balancing per symmetry block: each MPI rank gets one or more symmetry
      ! blocks to account for.
      !-------------------------------------------------------------------------
      if(activeblocks .ge. NPROCS) then
        ! More symmetry blocks than MPI ranks, i.e. we assign each rank
        ! one or more entire symmetry blocks
        if(mod(activeblocks, NPROCS) .ne. 0) then
          call stp('Incompatible number of MPI ranks for balancing_strategy = 1.')
        endif
        blocks_per_rank = activeblocks/NPROCS

        block_count = -1 ! unintuitive starting point: first block will be '0'
        do B=1,8
          if(blocks_global(B) .eq. 0) cycle
          block_count = block_count + 1
          if( block_count / blocks_per_rank .eq. MPI_rank) then
            ! attention, INTEGER division in the line above
            blocks_local(B) = blocks_global(B)
          endif
        enddo

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Find the first non-zero size in blocks_local
        do B=1,8
          if(blocks_local(B) .ne. 0) exit
        enddo
        offset = sum(blocks_global(1:B-1))

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Calculate the spwf <-> rank mapping and its inverse
        allocate(spwf_map(sum(blocks_local)))

        do i=1,sum(blocks_local)
          spwf_map(i)            = offset + i
          spwf_inverse(offset+i) = i
          rank_map(offset+i)     = MPI_RANK
        enddo
      else
        ! More ranks than blocks
        call stp('There are MPI ranks than symmetry blocks. Pick a different load-balancing strategy')
      endif
      !-------------------------------------------------------------------------
    case DEFAULT
      !-------------------------------------------------------------------------
      call stp('Unknown type of load balancing.')
      !-------------------------------------------------------------------------
    end select

#if(USE_MPI>0)
  ! these calls to allreduce is valid since we took care to zero things above
  call MPI_ALLREDUCE(MPI_IN_PLACE, rank_map, sum(blocks_global),     & 
  &                  MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpi_err)
  call MPI_ALLREDUCE(MPI_IN_PLACE, spwf_inverse, sum(blocks_global), & 
  &                  MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpi_err)

  call print_loadbalancing_information()
#endif

end subroutine loadbalance

#if(USE_MPI > 0)
  subroutine print_loadbalancing_information()
    !---------------------------------------------------------------------------------------
    ! Print detailed information on the balancing of spwfs between the different MPI ranks,
    ! BEFORE allocating any kind of memory.
    !---------------------------------------------------------------------------------------
    use vectors, only : memory_for_densities, memory

    integer              :: rank, B, row, col
    integer              :: mpi_err, tcount, neutron_ranks, proton_ranks, si, N
    integer, allocatable :: spwf_count(:)
    real(KIND=dp)        :: spwf_mem_local, den_mem, pot_mem, spwf_mem_2D

    1 format  (30('-'), ' MPI load balancing ', 30('-'))
    2 format  (' number of processes = ', i7)
    3 format ( '   Matrix blocking factors : ', i4, ' x ' i4)
    4 format  ('1D Layout')
    5 format  ('     B = ', i1, ' has ',  i4, ' MPI ranks for ', i7, ' spwfs in total.')

    6 format  ('2D Layout')
    7 format  ('     B = ', i1, ' has ',  i4, ' x ', i4, ' = ', i4' ranks (dropped = ', i4,')')

    9 format  (' Memory requirements')
   10 format  ('    Densities     = ', f10.3 , ' GB')
   11 format  ('    Potentials    = ', f10.3 , ' GB')
   12 format  ('    Spwfs (total) = ', f10.3 , ' GB')
   13 format  ('    Spwfs (max)   = ', f10.3 , ' GB')

   14 format ( ' Detailed information')
   15 format ( '   RANK  |  SYM_BLOCK    P     Q   #SPWFS | spwf_mem (GB) spwf_mem 2D (GB)')
   16 format ( 3x, i4, 2x, '|', 2x, i4, 6x, i4, 2x, i4, 3x, i4, 3x, '|', 3x, f10.3, 3x, f10.3 )

   99 format ( '--------------------------------------------------------------')

    tcount = sum(HFBlocks)
    if(MPI_rank .eq. 0) allocate(spwf_count(NPROCS))
    call MPI_gather(tcount,1,MPI_INTEGER,spwf_count,1,MPI_Integer, &
    &                      0,MPI_COMM_WORLD, mpi_err)

    if(MPI_RANK .eq. 0) then
      print 1
      print 2, NPROCS
      print 3, block_factor_row, block_factor_col

      print 99
      print 4
      print 99
      do B=1,8
        if(HFBLOCKS_GLOBAL(B) .eq. 0) cycle
        print 5, B, ranks_per_block(B), HFBlocks_global(B)
      enddo

      neutron_ranks = sum(ranks_per_block(1:4))
      proton_ranks  = sum(ranks_per_block(5:8))
      print 99
      print 6
      si = 0
      do B=1,8
        N = ranks_per_block(B)
        if(N.eq.0) cycle
        row = MAXVAL(MPI_2D_COORDINATES(si+1:si+N,1))+1
        col = MAXVAL(MPI_2D_COORDINATES(si+1:si+N,2))+1
        print 7, B, row, col, row*col, ranks_per_block(B) - row*col
        si = si + N
      enddo

      print 99
      print 9

      den_mem =transform_memory(memory_for_densities()) ! one set of densities
      ! several sets of potentials: F_in, F_out + memory * (Potential_iterates, Potential_updates)
      pot_mem = den_mem * 2 * (1 + memory)
      print 10, den_mem
      print 11, pot_mem
      print 12, transform_memory(memory_wavefunctions(sum(HFBlocks_global)))
      print 13, transform_memory(memory_wavefunctions(maxval(spwf_count)))

      print 99
      print 14
      print 15
      print 99
      do rank=1, NPROCS
        B = MPI_BLOCK_ASSIGNMENTS(rank)
        si = sum(ranks_per_block(1:B-1))
        N = ranks_per_block(B)
        row = MAXVAL(MPI_2D_COORDINATES(si+1:si+N,1))+1
        col = MAXVAL(MPI_2D_COORDINATES(si+1:si+N,2))+1
        spwf_mem_local = transform_memory(memory_wavefunctions(spwf_count(rank)))
        spwf_mem_2D    = transform_memory(memory_wavefunctions_2D(HFBLOCKS_GLOBAL(B), MPI_2D_COORDINATES(rank,1), row, MPI_2D_COORDINATES(rank,2), col))
        print 16, rank,B , &
        &          MPI_2D_COORDINATES(rank,1), MPI_2D_COORDINATES(rank,2), &
        &          spwf_count(rank), spwf_mem_local, spwf_mem_2D
      enddo
      print 99
      deallocate(spwf_count)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    end subroutine print_loadbalancing_information
#endif

  subroutine iniwavefunctions(ininx,ininy, ininz, ininwn, ininwp)   
    !---------------------------------------------------------------------------
    ! Build a set of initial wavefunctions in an EV8-like box, achieved through
    ! repeated calls to (pointer) subroutine initialise_wavefunctions.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Also initialized:
    !  *) Diagonal matrix elements of <h> = spenergies
    !  *) A default value for the sphamil as a diagonal matrix
    !  *) A default value for the hf_transfo as a trivial identity matrix
    ! 
    ! Not initialized here:
    !  *) Delta for the gaps. Since this module can not know what kind of 
    !     pairing is needed, it cannot correctly guess a structure. 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !    ininx/y/z : number of mesh points in 1/8th of the box
    !    ininwn    : number of neutron wavefunctions
    !    ininwp    : number of proton wavefunctions
    !---------------------------------------------------------------------------

    integer                   :: i,j
    integer, intent(in)       :: ininx, ininy, ininz, ininwn, ininwp
    integer                   :: ininwt
    integer, allocatable      :: kparz(:)
    integer                   :: nshells_ev
#if(USE_MPI > 0)
    integer                   :: xs, ys, B, wave, p,q, wave_global
    integer, external         :: NUMROC, INDXG2L, INDXG2P
#endif
    ininwt = ininwn + ininwp

    ! a) Generating the nilsson wave-functions in an EV8-box   
    if(allocated(spwf_map)) deallocate(spwf_map)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    ! First call of initialise_wavefunctions: this does not yet construct the 
    ! actual spwfs, but does get the parities in the array kparz.
    nshells_ev=max(11,int(1.5d0*max(ININWN,ININWP)**(1.d0/3.d0)))
    call initialise_wavefunctions(HFPsi,                                &
    &              kparz,spenergies,nshells_ev,nshells_ev-1,            &
    &              ININWT,ININWP,ININWN,floor(neutrons),floor(protons), &
    &              ININX,ININY,ININZ,dx,osc_freq,spwf_map)

    ! Based on the information in kparz, we construct the correct symmetry 
    ! properties and initialize the GLOBAL sizes of the symmetry blocks.
    hfblocks_global = 0
    do i=1,ININWN
        if(kparz(i) .gt. 0) HFBlocks_global(1) = HFBlocks_global(1) +1
        if(kparz(i) .lt. 0) HFBlocks_global(3) = HFBlocks_global(3) +1
    enddo
    do i=ININWN+1,ININWT
        if(kparz(i) .gt. 0) HFBlocks_global(5) = HFBlocks_global(5) +1
        if(kparz(i) .lt. 0) HFBlocks_global(7) = HFBlocks_global(7) +1
    enddo
    ! then, we are capable of figuring out the way to balance the spwfs among
    ! the different MPI ranks. 
    call loadbalance(HFblocks_global,balancing_strategy,& 
    &                        HFblocks,spwf_map,rank_map, spwf_inverse)
    ! now each MPI rank knows which spwfs it should grab and can allocate 
    ! the required space. 
    allocate(HFPSI(ININX*ININY*ININZ,4,sum(HFblocks))); hfpsi = 0.0d0

    ! and finally, we can call initialise_wavefunctions a second time in order 
    ! actually the requested spwfs in this partical process.
    call initialise_wavefunctions(HFPsi,                                &
    &              kparz,spenergies,nshells_ev,nshells_ev-1,            &
    &              ININWT,ININWP,ININWN,floor(neutrons),floor(protons), &
    &              ININX,ININY,ININZ,dx,osc_freq,spwf_map)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! b) and now we go on to populate more symmetry information
    allocate(sx(4,sum(hfblocks)), sy(4,sum(hfblocks)), sz(4,sum(hfblocks)))
    do i=1, HFBlocks(1)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
    enddo

    do i=HFBlocks(1) + 1,HFBlocks(1) + HFBlocks(3)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
    enddo

    do i=HFBlocks(1) + HFBlocks(3)+1,HFBlocks(1) + HFBlocks(3) +HFBlocks(5)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
    enddo
    do i=HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + 1,                      &
    &       HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + HFBLocks(7)
        sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1 
        sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! c) and perform some other initializations
    allocate(dispersions(ININWT)) ; dispersions  = 0

#if(USE_MPI == 0) 
    allocate(hftransfo(nwt,nwt), sphamil(nwt,nwt))
    do i=1, nwt
      hftransfo(i,i) = 1.0d0
      do j=i+1,nwt
        hftransfo(i,j) = 0.0d0 
        hftransfo(j,i) = 0.0d0 
      enddo
    enddo

    do i=1, nwt
      sphamil(i,i) = spenergies(i)
      do j=i+1,nwt
        sphamil(i,j) = 0.0d0 
        sphamil(j,i) = 0.0d0 
      enddo
    enddo
#else 
    B = MPI_BLOCK_SIZE
    ! Both hftransfo and sphamil get distributed on the 2D layout
    xs = NUMROC(B,BLOCK_FACTOR_ROW,MYROW_2D,0,NROW_2D)
    ys = NUMROC(B,BLOCK_FACTOR_COL,MYCOL_2D,0,NCOL_2D)

    allocate(hftransfo(xs,ys), sphamil(xs,ys))
    hftransfo = 0.0d0 ; sphamil   = 0.0d0
    
    do wave=1,MPI_BLOCK_SIZE
      ! Identify which of the ranks has information on this particular spwf
      p = INDXG2P(wave, BLOCK_FACTOR_ROW, MYROW_2D, 0, NROW_2D)
      q = INDXG2P(wave, BLOCK_FACTOR_COL, MYCOL_2D, 0, NCOL_2D)
      ! ... and to which local matrix element the (global) diagonal matrix
      !     element corresponds.
      i = INDXG2L(wave, BLOCK_FACTOR_ROW, MYROW_2D, 0, NROW_2D)
      j = INDXG2L(wave, BLOCK_FACTOR_COL, MYCOL_2D, 0, NCOL_2D)
      ! ... but we also need to know what is the GLOBAL index of the spwf for 
      !     the spenergies array
      wave_global = INDXG2L(wave,BLOCK_FACTOR_1D, 0, MYCOL_1D, NCOL_1D) &
      &           + sum(HFBLOCKS_global(1:MPI_SYM_BLOCK-1))
      if((MYROW_2D .eq. p) .and. (MYCOL_2D .eq. q)) then
        if(i.ne.0 .and. j.ne.0) then
          sphamil(i,j)   = spenergies( + wave_global)
          hftransfo(i,j) = 1.0d0
        endif
      endif
    enddo
#endif
  end subroutine iniwavefunctions

  subroutine randomspwfs(psi,par,spe,nshells_even, nshells_odd, &
  &                      nw,nwp,nwn,neut,prot,mx,my,mz,dx,osc_freq,map)
     !--------------------------------------------------------------------------
     ! Generate a set of single-particle wavefunctions randomly. 
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! Input:
     !    nw, nwp, nwn : total/proton/neutron number of spwfs to construct
     !    map          : If allocated, this routine will construct the spwfs.
     !                   Otherwise unused.
     !    mx,my,mz,dx,osc_freq, nshells_even, nshells_odd, neut, prot :
     !            -> dummy arguments to make this routines calling signature
     !               identical to that of nilsson.
     ! Output:
     !    psi: a set of wavefunctions with random values. Only initialised if 
     !         the map input is allocated.
     !         ATTENTION: this set is not orthonormal at all. 
     !    spe: a guess (trivial in this routine) of the single-particle energies
     !         of these random states. 
     !    par: the parity quantum numbers of the spwfs
     !--------------------------------------------------------------------------
     real(KIND=dp), allocatable, intent(inout) :: psi(:,:,:), spe(:)
     real(KIND=dp), intent(in)                 :: dx, osc_freq(3)
     integer, intent(inout), allocatable       :: par(:)
     integer, intent(in)        :: nw,nwn,nwp, mx, my, mz, neut, prot
     integer, intent(in)        :: nshells_even, nshells_odd
     integer, intent(in), allocatable :: map(:)
     
     integer :: s
     integer, allocatable :: seed(:)

     if(allocated(par)) deallocate(par)
     if(allocated(spe)) deallocate(spe)
     allocate(par(nw),spe(nw))
     
     ! We just take spwfs of positive and negative in 50/50 proportion
     par(            1:nwn/2    ) = +1
     par(nwn/2      +1:nwn      ) = -1
     par(nwn        +1:nwn+nwp/2) = +1
     par(nwn  +nwp/2+1:nw       ) = -1
     
     spenergies = 100
     
     if(allocated(map)) then
      call random_seed(size=s)
      allocate(seed(s))
      seed = 987654321 + MPI_RANK * 123456789 ! Seed value needs to depend on
                                              ! MPI_RANK; if not, we all ranks 
                                              ! will generate the same sequence
                                              ! and we will run in trouble with
                                              ! orthonormalisation
      call random_seed(put=seed)
      call random_number(psi) ! randomize
     endif
  end subroutine randomspwfs

  !subroutine planewaves(psi,par,spe,nwaves,nwavesp,nwavesn, &
  !  &                                      mx,my,mz,dx,map)
  !  
  !  real(KIND=dp), allocatable :: psi(:,:,:)
  !  integer, intent(inout)     :: par(:), spe(:)
  !  integer, intent(in)        :: nwaves, nwavesp, nwavesn, mx, my, mz, dx
  !  integer, intent(in), allocatable :: map(:)
  !  
  !  real(KIND=dp)        :: kx, ky, kz
  !  integer              :: nmaxx, nmaxy, nmaxz
  !  integer, allocatable :: waves(:,:)
  !  
  !  !---------------------------------------------------------------------------
  !  ! wavenumber multiplicators
  !  kx = pi/(2*mx*dx)
  !  ky = pi/(2*my*dx)
  !  kz = pi/(2*mz*dx)
  !
  !  nmaxx = max(11,int(1.5d0*max(ININWN,ININWP)**(1.d0/3.d0)))
  !  nmaxy = max(11,int(1.5d0*max(ININWN,ININWP)**(1.d0/3.d0)))
  !  nmaxz = max(11,int(1.5d0*max(ININWN,ININWP)**(1.d0/3.d0)))!

  !  ! First even parity
  !  allocate(waves(nmaxx,nmaxy,nmaxz,4))

    ! then odd parity
  
  !end subroutine planewaves()

  subroutine deriveHF()
    !---------------------------------------------------------------------------
    ! Derives all of the single-particle wave-functions in the basis in memory
    ! (which is not always the actual HF basis)
    !---------------------------------------------------------------------------
    integer :: wave,k
    
    call start_timer(T_derivatives)

    if(.not.allocated(HFdPsi)) then
        allocate(HFdPsi(nx*ny*nz,3,4,nwt_local))
        allocate(HFddPsi(nx*ny*nz,6,4,nwt_local))
    endif

$N3    if(.not.allocated(HFdddpsi)) then
$N3        allocate(HFdddPsi(nx*ny*nz,10,4,nwt_local))
$N3    endif

#if(USE_Periodic==0)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Original Lagrange mesh derivatives
    do wave=1,nwt_local
        do k=1,4
$N2        call Derive_tot(HFPsi(:,k,wave), sx(k,wave), sy(k,wave), sz(k,wave),&
$N2        &                                           HFdPsi(:,:,k,wave),     &
$N2        &                                           HFddPsi(:,:,k,wave))

$N3        call Derive_tot(HFPsi(:,k,wave), sx(k,wave), sy(k,wave), sz(k,wave),&
$N3        &                                           HFdPsi(:,:,k,wave),     &
$N3        &                                           HFddPsi(:,:,k,wave),    &
$N3        &                                           HFdddPsi(:,:,k,wave))
        enddo
    enddo
    
#else    
    !---------------------------------------------------------------------------
    ! NS: for periodic boundary conditions the spwfs can no longer be derived
    !     component-by-component. Real and imaginary parts of spin-up and 
    !     spin-down are connected!
    ! Not ready for N3LO calculations
    !---------------------------------------------------------------------------
    do wave=1,nwt_local
        do k=1,2
$N2        call derive_tot_periodic(HFPsi(:,   (2*k-1):2*k,wave),              &
$N2        &                                sx((2*k-1):2*k,wave),              &
$N2        &                                sy((2*k-1):2*k,wave),              &
$N2        &                                sz((2*k-1):2*k,wave),              &
$N2        &                          HFdPsi(:,:,(2*k-1):2*k,wave),            &
$N2        &                            HFddPsi(:,:,(2*k-1):2*k,wave))

$N3        call stp('N3LO functionals not yet supported for periodic BCs.')
        enddo
    enddo
#endif

    call stop_timer(T_derivatives)
  end subroutine DeriveHF

  
  subroutine derive_extra_spwfs(extraspwfs)
      !-------------------------------------------------------------------------
      ! Derives all of the single-particle wave-functions that were added as
      ! "bonus". Useful if these spwfs are evolved separately from the rest.
      !-------------------------------------------------------------------------
  
      integer, intent(in) :: extraspwfs(8)
      integer :: wave,k, B, si, N
      
      si = 0
      do B=1,8
        N = HFBlocks(B) ; if(N.eq.0) cycle
        
        do wave=si+N-extraspwfs(B)+1, si+N
          do k=1,4
$N2        call Derive_tot(HFPsi(:,k,wave), sx(k,wave), sy(k,wave), sz(k,wave),&
$N2        &                                           HFdPsi(:,:,k,wave),     &
$N2        &                                           HFddPsi(:,:,k,wave))

$N3        call Derive_tot(HFPsi(:,k,wave), sx(k,wave), sy(k,wave), sz(k,wave),&
$N3        &                                           HFdPsi(:,:,k,wave),     &
$N3        &                                           HFddPsi(:,:,k,wave),    &
$N3        &                                           HFdddPsi(:,:,k,wave))
          enddo
        enddo
        si = si + N
      enddo
  
  end subroutine derive_extra_spwfs
  
  subroutine deriveCan()
    !---------------------------------------------------------------------------
    ! Derives all of the single-particle wave-functions in the canonical basis.
    !---------------------------------------------------------------------------
    integer :: wave,k

    call start_timer(T_derivatives_can)

    if(allocated(CanPsi)) then
      if(.not.allocated(CANdPsi)) then
          allocate( CANdPsi(nx*ny*nz,3,4,nwt_local))
          allocate(CANddPsi(nx*ny*nz,6,4,nwt_local))
      endif
    endif

$N3    if(allocated(CanPsi)) then
$N3       if(.not.allocated(CANdddpsi)) then
$N3         allocate(CandddPsi(nx*ny*nz,10,4,nwt_local))
$N3       endif
$N3    endif

    if(allocated(CanPsi)) then
      do wave=1,nwt_local
        do k=1,4
$N2        call Derive_tot(CANPsi(:,k,wave), sx(k,wave),sy(k,wave), sz(k,wave),&
$N2        &                                           CANdPsi(:,:,k,wave),    &
$N2        &                                           CANddPsi(:,:,k,wave))

$N3        call Derive_tot(CANPsi(:,k,wave), sx(k,wave),sy(k,wave), sz(k,wave),&
$N3        &                                           CANdPsi(:,:,k,wave),    &
$N3        &                                           CANddPsi(:,:,k,wave),   &
$N3        &                                           CANdddPsi(:,:,k,wave))

        enddo
      enddo
    endif
    call stop_timer(T_derivatives_can)
    
  end subroutine DeriveCan
  
  function OrderSpwfsISO(Isospin, canonical) result(Indices)
    !---------------------------------------------------------------------------
    ! Orders the wavefunctions within an isospin block. 
    !---------------------------------------------------------------------------
    integer, intent(in)        :: Isospin
    
    integer, allocatable       :: Indices(:)
    real(Kind=dp), allocatable :: Energies(:)
    integer                    :: i, nwf,  HolePos, ToInsertIndex
    real(Kind=dp)              :: ToInsert
    logical, intent(in), optional :: canonical
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Count the number of proton and neutron wavefunctions
    if (Isospin .eq. -1) then
        nwf = nwn
    else
        nwf = nwp
    endif
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Filling Energies & Indices
    allocate(Indices(nwf), Energies(nwf))
    do i=1,nwf
       Indices(i) = i 
    enddo

    if(.not. present(canonical)) then
      if(Isospin.eq.-1) then
          Energies = spenergies(1:nwn)
      else
          Indices  = Indices + nwn
          Energies = spenergies(nwn+1:nwn+nwp)
      endif
    elseif(canonical) then
      if(Isospin.eq.-1) then
          Energies = canenergies(1:nwn)
      else
          Indices  = Indices + nwn
          Energies = canenergies(nwn+1:nwn+nwp)
      endif
    else
      if(Isospin.eq.-1) then
          Energies = spenergies(1:nwn)
      else
          Indices  = Indices + nwn
          Energies = spenergies(nwn+1:nwn+nwp)
      endif
    endif
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Sort the energies
    do i=2,nwf
      !Make a hole at index i
      ToInsert = Energies(i)
      HolePos  = i
      ToInsertIndex = Indices(i)
      do while(ToInsert.lt.Energies(HolePos-1))
        !Move the hole one place down
        Energies(HolePos) = Energies(HolePos-1)
        Indices(HolePos) = Indices(HolePos-1)
        HolePos = HolePos - 1
        if(HolePos.eq.1.0_dp) exit
      enddo
      !Insert the energy at the correct place
      Energies(HolePos) = ToInsert
      Indices(HolePos)  = ToInsertIndex
    enddo

    ! Making sure these are not saved for a next call
    deallocate(Energies)
  end function OrderSpwfsISO
  
  function OrderSpwfsSym(block) result(indices)
    !---------------------------------------------------------------------------
    ! Sort the single-particle wave-functions in the given symmetry-block 
    ! by single-particle energy. Note, this routine works with the indices 
    ! LOCAL to any given MPI rank.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !      block   : integer
    !                index of the symmetry block to consider.
    ! Output:
    !      indices : allocatable, integer.
    !                indices of the spwfs, in increasing order
    !---------------------------------------------------------------------------
    integer, intent(in)        :: block
    integer, allocatable       :: Indices(:)
    real(KIND=dp), allocatable :: Energies(:)
    integer                    :: nwf, startind, HolePos, ToInsertIndex, i
    real(Kind=dp)              :: ToInsert
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Count the number of relevant wavefunctions
    nwf      = HFBlocks(block) 
    startind = sum(HFBlocks(1:block-1))
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Filling Energies & Indices
    allocate(Indices(nwf), Energies(nwf))
    do i=1,nwf
       Indices(i)  = startind + i 
       Energies(i) = spwf_map(i) 
    enddo
    
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Sort the energies
    do i=2,nwf
      !Make a hole at index i
      ToInsert = Energies(i)
      HolePos  = i
      ToInsertIndex = Indices(i)
      do while(ToInsert.lt.Energies(HolePos-1))
        !Move the hole one place down
        Energies(HolePos) = Energies(HolePos-1)
        Indices(HolePos) = Indices(HolePos-1)
        HolePos = HolePos - 1
        if(HolePos.eq.1.0_dp) exit
      enddo
      !Insert the energy at the correct place
      Energies(HolePos) = ToInsert
      Indices(HolePos)  = ToInsertIndex
    enddo
    deallocate(energies,indices)
  end function OrderSpwfsSym

  subroutine set_spwf_symmetries(sx, sy, sz, blocks)
    !---------------------------------------------------------------------------
    ! Routine that assigns the correct "signs" under symmetry transformation 
    ! to all single-particle wavefunctions. 
    !
    ! Hephaestos fills in all practical signs here, using the following key
    !
    !  [dollar sign] SXab
    !
    !  where 
    !      X  = X;Y;Z depending on the direction
    !      a  = block index, i.e. 1-4 depending on the symmetries of the spwf
    !      b  = component index, i.e. 1-4 depending on the spinor component
    !           we are dealing with
    !---------------------------------------------------------------------------
    ! Input: blocks 
    !        The number of spwfs in every symmetry block. 
    !        This is an input of this routine, because this subroutine could be 
    !        called for a subset of spwfs when reading from file.
    !
    ! Output: sx, sy, sz
    !        The signs with respect to x/y/z reflection of the components of 
    !        the spwf spinors. 
    !---------------------------------------------------------------------------
    integer, intent(out), allocatable :: sx(:,:), sy(:,:), sz(:,:)
    integer, intent(in)               :: blocks(8)
    integer                           :: N, i, B, si

    if(allocated(sx)) deallocate(sx)
    if(allocated(sy)) deallocate(sy)
    if(allocated(sz)) deallocate(sz)

    N = sum(blocks)
    allocate(sx(4,N), sy(4,N), sz(4,N))

    si = 0
    do B=1,8,4 ! This is simply an isospin loop
  
      ! Positive parity neutrons
      do i=si+1, si+blocks(B)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! Block 1: positive parity, signature = +i spwfs
        !
        ! EV8-values for these quantities, as a point of comparison.
        !
        ! sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        ! sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1 
        ! sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        ! sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          sx(1,i) = $SX11 ; sy(1,i) = $SY11 ; sz(1,i) = $SZ11
          sx(2,i) = $SX12 ; sy(2,i) = $SY12 ; sz(2,i) = $SZ12
          sx(3,i) = $SX13 ; sy(3,i) = $SY13 ; sz(3,i) = $SZ13
          sx(4,i) = $SX14 ; sy(4,i) = $SY14 ; sz(4,i) = $SZ14
      enddo
      do i=si+blocks(B)+1, si+blocks(B)+blocks(B+1)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  Block 2: positive parity, signature=-i spwfs
        !
        ! EV8-values for these quantities, as a point of comparison.
        !
        ! sx(1,i) = -1 ; sy(1,i) = +1 ; sz(1,i) = -1
        ! sx(2,i) = +1 ; sy(2,i) = -1 ; sz(2,i) = -1 
        ! sx(3,i) = +1 ; sy(3,i) = +1 ; sz(3,i) = +1
        ! sx(4,i) = -1 ; sy(4,i) = -1 ; sz(4,i) = +1
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          sx(1,i) = $SX21 ; sy(1,i) = $SY21 ; sz(1,i) = $SZ21
          sx(2,i) = $SX22 ; sy(2,i) = $SY22 ; sz(2,i) = $SZ22
          sx(3,i) = $SX23 ; sy(3,i) = $SY23 ; sz(3,i) = $SZ23
          sx(4,i) = $SX24 ; sy(4,i) = $SY24 ; sz(4,i) = $SZ24
      enddo
      ! Negative parity neutrons
      do i=si+sum(blocks(B:B+1))+1,si+sum(blocks(B:B+2))
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  Block 3: negative parity, signature = +i spwfs
        !
        ! EV8-values for these quantities, as a point of comparison.
        !
        !  sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        !  sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1 
        !  sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        !  sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          sx(1,i) = $SX31 ; sy(1,i) = $SY31 ; sz(1,i) = $SZ31
          sx(2,i) = $SX32 ; sy(2,i) = $SY32 ; sz(2,i) = $SZ32
          sx(3,i) = $SX33 ; sy(3,i) = $SY33 ; sz(3,i) = $SZ33
          sx(4,i) = $SX34 ; sy(4,i) = $SY34 ; sz(4,i) = $SZ34
      enddo
      do i=si+sum(blocks(B:B+2))+1, si+sum(blocks(B:B+3))
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  Block 4: negative parity, signature = -i spwfs
        !
        ! EV8-values for these quantities, as a point of comparison.
        !  sx(1,i) = -1 ; sy(1,i) = +1 ; sz(1,i) = +1
        !  sx(2,i) = +1 ; sy(2,i) = -1 ; sz(2,i) = +1 
        !  sx(3,i) = +1 ; sy(3,i) = +1 ; sz(3,i) = -1
        !  sx(4,i) = -1 ; sy(4,i) = -1 ; sz(4,i) = -1
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          sx(1,i) = $SX41 ; sy(1,i) = $SY41 ; sz(1,i) = $SZ41
          sx(2,i) = $SX42 ; sy(2,i) = $SY42 ; sz(2,i) = $SZ42
          sx(3,i) = $SX43 ; sy(3,i) = $SY43 ; sz(3,i) = $SZ43
          sx(4,i) = $SX44 ; sy(4,i) = $SY44 ; sz(4,i) = $SZ44
      enddo
      si = si + sum(blocks(B:B+3))
    enddo 

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Also set the symmetry properties of the max_spwf employed in evolution.f90
    ! These are simply a copy of those in the very first symmetry block
    sx_max(1) =  $SX11 ; sy_max(1) = $SY11 ; sz_max(1) = $SZ11
    sx_max(2) =  $SX12 ; sy_max(2) = $SY12 ; sz_max(2) = $SZ12
    sx_max(3) =  $SX13 ; sy_max(3) = $SY13 ; sz_max(3) = $SZ13
    sx_max(4) =  $SX14 ; sy_max(4) = $SY14 ; sz_max(4) = $SZ14

  end subroutine set_spwf_symmetries

!===============================================================================
! Orthonormalisation routines
! 
  subroutine GramSchmidt
    !---------------------------------------------------------------------------
    ! This subroutine uses a (modified) Gram-Schmidt scheme to orthonormalise 
    ! the spwfs in the array HFPsi. The orthonormalisation proceeds per symmetry
    ! block, as this saves precious CPU cycles.
    !
    ! In the interest of convergence speed, the orthogonalisation is done in 
    ! order of ascending single-particle energy if this is possible, i.e. if
    ! diagsphamil == .true..
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Note: this routine is written in a way such that it does not care about
    !       the spatial dimensions of the HFPSI array. This way, it can be 
    !       called both during a calculation and/or during its set-up when
    !       HFPsi can perhaps be defined on a smaller mesh. Since this routine
    !       cannot infer the value of dv in the latter case, it might be that
    !       the resulting spwfs are not completely normalized.
    !---------------------------------------------------------------------------
    integer  :: b, i,j,nw, mw, si, N
    integer  :: indices(maxval(HFBlocks))

    real(KIND=dp) ::  norm

#if(USE_MPI>0)
    if(balancing_strategy.ne.1) then
      call stp('Balancing_strategy should be 1 for GramSchmidt to work.')
    endif
#endif

    call start_timer(T_ortho)

    si = 0
    do b = 1, Blocks 
        N = HFBlocks(B) ; if(N.eq.0) cycle

        indices = 0
        if(diagsphamil) then
          indices(1:HFblocks(b)) = OrderSpwfsSym(b)
          ! Note, in the case of MPI calculations OrderSpwfsSym deals with
          !       indices that are local to an MPI-rank already.
        else
          do i=1, N
            indices(i) = si + i
          enddo
        endif

        do i = 1,N
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Normalize wave-function nw
            nw = indices(i)
            norm = sum(HFpsi(:,:,nw)**2) * dv
            HFPsi(:,:,nw) = (sqrt(1.0/norm)) * HFPsi(:,:,nw) 
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Then subtract the projection on \Psi_{nw} from all the following
            ! Spwf.
            ! Re(\Psi(\sigma)_{mw}) = Re(\Psi(\sigma)_{mw})
            !                - Re(< \Psi_{nw}|\Psi_{mw} >) Re(\Psi(\sigma)_{nw})
            !                + Im(< \Psi_{nw}|\Psi_{mw} >) Im(\Psi(\sigma)_{nw})
            ! Im(\Psi(\sigma)_{mw}) = Im(\Psi(\sigma)_{mw})
            !                - Re(< \Psi_{nw}|\Psi_{mw} >) Im(\Psi(\sigma)_{nw})
            !                - Im(< \Psi_{nw}|\Psi_{mw} >) Re(\Psi(\sigma)_{nw})
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Note that the imaginary part of the inproduct only needs to be 
            ! taken into account when there is no antilinear, hermitian 
            ! symmetry that is conserved.
            !
            ! THIS IS NOT IMPLEMENTED YET HOWEVER!
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! What is also missing is an orthonormalisation versus the spwf
            ! that are assumed to be present but not represented numerically.
            ! The MOCCa example is conserved time-reversal but broken signature.
            !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            do j= i+1, HFBlocks(b)
              mw = indices(j)    
              ! Real part of the inproduct
              norm = sum(HFpsi(:,:,nw)*HFpsi(:,:,mw)) * dv
              HFPsi(:,:,mw) = HFPsi(:,:,mw) - norm * HFPsi(:,:,nw)
            enddo
            norm = sum(HFpsi(:,:,nw)**2) * dv
        enddo
        si = si + N
    enddo
    call stop_timer(T_ortho)

  end subroutine GramSchmidt

  subroutine Cholesky_orthonormalisation()
    !---------------------------------------------------------------------------
    ! Orthonormalize the single-particle wavefunctions in the hfpsi array 
    ! through a Cholesky decomposition.
    ! 
    ! Let N be the real positive-definite matrix of the overlaps of the 
    ! single-particle wavefunctions. Then we can write
    !        N^T N = L L^T
    ! where L is a real lower triangular matrix with positive diagonal entries.
    !
    ! Then X = N (L^-1)^T is an orthonormal set of vectors, since 
    !
    !  X^T X = L^-1 N^T N [L^-1]^T 
    !        = L^-1 L L^T [L^-1]^T 
    !        = (L^-1 L) (L^-1 L)^T
    !        = identity matrix. 
    ! 
    ! Once L is constructed, one can obtain X easily through solving a set of 
    ! linear equations:
    !      X L^T = N 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Technical note: I am for now happily ignorant about any and all kind of
    !                 numerical instabilities that Cholesky factorisation implies.
    !                 I know these exist in theory, but have so far to the best
    !                 of my knowledge not encountered them while running this code.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !    None
    ! Output:
    !    None
    ! Side-effects:
    !    the contents of the array hfpsi get changed and now are orthonormal.
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable         :: overlaps(:,:)
#if(USE_MPI == 0)
    real(KIND=dp), pointer, contiguous :: wfs_reshape(:,:)
#endif
    integer                    :: N,si, B, info
#if(USE_MPI > 0)
    integer                    :: xs, ys
    integer, external          :: NUMROC
#endif

#if(USE_MPI > 0)
    if(MYROW_2D .eq.-1) then  ! this particular rank is not part of the 2D layout
        allocate(overlaps(1,1)); deallocate(overlaps) 
        ! allocate/deallocate to ensure the cray compiler does not complain 
        return
    endif
#endif

    call start_timer(T_ortho)
 
    si = 0
    do B=1,8
      N = HFBlocks(B) 

#if(USE_MPI == 0)
      if (N.eq.0) cycle   ! This can only be explicitly put here.
      allocate(overlaps(N,N))
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Reshaping the wavefunctions by pointer in order to make the BLAS calls
      ! as efficient and easy as possible. Note: the "contiguous" keyword for
      ! this pointer array is crucial to make this trick work without tripping
      ! boundary-checking by compilers. 
      wfs_reshape(1:4*mv,1:N) => hfpsi(1:mv,1:4,si+1:si+N)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Build overlaps within this symmetry block 
      call start_timer(T_norm_ortho)
      call dgemm('t','n',N,N,4*mv, dv,wfs_reshape,4*mv, wfs_reshape, 4*mv, &
      &                            0.0d0, overlaps,N)
      call stop_timer(T_norm_ortho)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Calculate the cholesky decomposition
      call start_timer(T_diag_ortho)
      call dpotrf('l',N,overlaps,N,info)
      ! overlaps now contains the factorisation L
      if(info.ne.0) then
        print *, 'Issue with the Cholesky decomposition.'
        print *, 'INFO = ', info
      endif
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Solve the linear equations
      !    X L^T = N
      call dtrsm('r','l','t','n',4*mv,N,1.0d0,overlaps,N,wfs_reshape,4*mv)
      call stop_timer(T_diag_ortho)
#else
      ! The MPI version cannot cycle over blocks of size 0, since the MPI
      ! ranks might not be relevant to the calculation but at the same time
      ! might be part of a given BLACS context and hence are required to 
      ! call all these SCALAPACK routines; if not, we get stuck eternally...
      ! if (N.eq.0) cycle   ! This can only be explicitly put here.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      call start_timer(T_norm_ortho)
      ! Build overlaps within this symmetry block
      xs = NUMROC(MPI_BLOCK_SIZE,BLOCK_FACTOR_ROW,MYROW_2D,0,NROW_2D)
      ys = NUMROC(MPI_BLOCK_SIZE,BLOCK_FACTOR_COL,MYCOL_2D,0,NCOL_2D)
      allocate(overlaps(xs,ys))
      call PDGEMM ('T', 'N', MPI_BLOCK_SIZE, MPI_BLOCK_SIZE, 4*mv, dv, &
      &             HFpsi_2d, 1, 1, desc_psi_2d, &
      &             HFpsi_2d, 1, 1, desc_psi_2d, &
      &             0.0d0,                       &
      &             overlaps, 1, 1, desc_mat_2d)
      call stop_timer(T_norm_ortho)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Calculate the cholesky decomposition
      call start_timer(T_diag_ortho)
      CALL PDPOTRF('l',MPI_BLOCK_SIZE,overlaps,1,1,desc_mat_2d, info)
      ! overlaps now contains the factorisation L
      if(info.ne.0) then
        print *, 'Issue with the Cholesky decomposition.'
        print *, 'INFO = ', info
      endif
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Solve the linear equations
      !    X L^T = N
      CALL PDTRSM('r','l','t','n',4*mv,MPI_BLOCK_SIZE,1.0d0,           &
      &                                overlaps, 1, 1, desc_mat_2d, &
      &                                hfpsi_2D, 1, 1, desc_psi_2d)
      call stop_timer(T_diag_ortho)
#endif
      deallocate(overlaps)
      si = si +N
    enddo
    call stop_timer(T_ortho)
  end subroutine Cholesky_orthonormalisation
!===============================================================================
! MPI-only routines
!===============================================================================
#if(USE_MPI > 0)
  subroutine transfer_1D_to_2D(A_1D, A_2D)
    !---------------------------------------------------------------------------
    ! Use the scalapack routine PDGEMR2D to copy the 1D-distributed array into
    ! the 2D-distributed version A_2D. If necessary, A_2D will get allocated.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input :
    !    A_1D : array (dimension (mv,4,X)), 1D-layout distributed
    ! Output:
    !    A_2D : array (dimension (mv*4,X)), copy of A_1D remapped with a pointer
    !           but distributed among processes in a 2D layout.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) , contiguous, target :: A_1D(:,:,:)
    real(KIND=dp), intent(out), allocatable        :: A_2D(:,:)
    ! Pointer to remap the (mv,4,X) array into a (4*mv,X) array
    real(KIND=dp), pointer, contiguous             :: A_1Dc(:,:)

    integer           :: xs, ys
    integer, external :: NUMROC

    call start_timer(T_transfer_psi)

    if(.not.allocated(A_2D) ) then
      ! Asking for the appropriate size of the A_2D matrix on this process
      ! (if MY_ROW2D = -1,then the node is not part of the 2D grid)
      if(MYROW_2D .ne. -1) then
        xs = NUMROC(          4*mv,block_factor_row,MYROW_2D,0,NROW_2D)
        ys = NUMROC(MPI_BLOCK_SIZE,block_factor_col,MYCOL_2D,0,NCOL_2D)
      else
        xs = 1
        ys = 1
      endif
      allocate(A_2D(xs,ys))
    endif

    A_1Dc(1:4*mv, 1:nwt_local) => A_1D
    call pdgemr2d(4*mv,MPI_BLOCK_SIZE,A_1Dc, 1,1, desc_psi_1D,                 &
     &                                A_2D  ,1,1, desc_psi_2D, blacs_cntxt_1D)

    call stop_timer(T_transfer_psi)

  end subroutine transfer_1D_to_2D
  
  subroutine transfer_2D_to_1D(A_2D, A_1D)
    !---------------------------------------------------------------------------
    ! Use the scalapack routine PDGEMR2D to copy the 2D-distributed array into
    ! the 1D-distributed version A_1D. This routine is the inverse of the one
    ! above, except it assumes all relevant arrays are allocated.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input :
    !    A_2D : array (dimension (4*mv,X)), 2D-layout distributed
    ! Output:
    !    A_1D : array (dimension (mv,4,X)), copy of A_2D remapped with a pointer
    !           but distributed among processes in a 1D layout.
    !--------------------------------------------------------------------------- 
    real(KIND=dp), intent(in)                        :: A_2D(:,:)
    real(KIND=dp), intent(inout), contiguous, target :: A_1D(:,:,:) 
    ! ^- INOUT attribute, since otherwise the compiler deallocates stuff
    ! Pointer for remapping 
    real(KIND=dp), pointer, contiguous :: A_1Dc(:,:)

    call start_timer(T_transfer_psi)

    ! Pointer remapping 
    A_1Dc(1:4*mv, 1:nwt_local) => A_1D
    call pdgemr2d(4*mv,MPI_BLOCK_SIZE,A_2D ,1,1, desc_psi_2D,                  &
    &                                 A_1Dc,1,1, desc_psi_1D, blacs_cntxt_1D)

    call stop_timer(T_transfer_psi)
  end subroutine transfer_2D_to_1D

#endif
!===============================================================================

  function TimeReverse(psi) result(Tpsi)
    !---------------------------------------------------------------------------
    ! Perform a time-reversal on the input spinor.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Psi(mv,4)
    real(KIND=dp)             :: Tpsi(mv,4)
    
    TPsi(:,1) =   Psi(:,3)
    TPsi(:,2) = - Psi(:,4)
    TPsi(:,3) = - Psi(:,1)
    TPsi(:,4) =   Psi(:,2)
    
  end function TimeReverse
  
  subroutine update_spwf_symmetries(fullmatrices)
      !-------------------------------------------------------------------------      
      ! Update/calculate all relevant expectation values of single-particle 
      ! wavefunctions, for both HF-basis and canonical basis.
      !
      ! Currently implemented:
      !      *   Parity: P_HF and P_can
      !
      ! Input: 
      !    NONE
      ! Output:
      !    NONE
      !-------------------------------------------------------------------------
      real(KIND=dp), allocatable :: full_P(:,:)
      logical, intent(in)        :: fullmatrices
      integer                    :: i

      if(.not. allocated(P_HF))  allocate(P_HF(nwt))
#if(PASTA == 1)
      allocate(full_P(1,1)); deallocate(full_P)
      ! This is a waste of CPU time for pasta calculations. 
      ! This return is ugly and will require more elegant inclusion later on.
      return
#endif

      allocate(full_P(nwt,nwt))
      if(allocated(canenergies) .and. (.not.allocated(P_CAN))) allocate(P_CAN(nwt))

      full_P  = spwf_parities(HFPsi, fullmatrices)
      if(.not. diagsphamil) then
        full_P  = matmul(full_P, HFtransfo)
        full_P  = matmul(transpose(HFtransfo), full_P)
      endif

      do i=1,nwt
        P_HF(i) = full_P(i,i)
      enddo

      if(allocated(canenergies)) then
        full_P = spwf_parities(denpsi,.false.)
        do i=1,nwt
          P_CAN(i) = full_P(i,i)
        enddo
      endif
    
      deallocate(full_P)
     
  end subroutine update_spwf_symmetries

  subroutine update_spwf_angmom(fullmatrices)
    !---------------------------------------------------------------------------
    ! Calculate all relevant single-particle matrix elements of 
    ! 
    !  (a) Jx, Jy, Jz
    !  (b) JxT, JyT, JzT => Real (JTR) and imaginary (JTI) parts
    !  (c) Jx^2, Jy^2, Jz^2
    !  (d) JJ 
    !  (e) Sx, Sy, Sz
    !  (f) STx, STy, STz
    !  
    ! where JJ is a simple number, such that J*(J+1) = Jx^2 + Jy^2 + Jz^2.
    !
    ! These things are calculated for 
    !  (1) the spwfs in memory    => direct integration over the box 
    !  (2) the Hartree-fock basis => matrix transformation
    !                                (only if fullmatrices == .true.)
    !  (3) the canonical basis    => direct integration over the box
    !
    ! Right now, all these things are calculated in CR8-like geometry.
    !
    ! Input: 
    !    fullmatrices : if .true., the full matrix elements are calculated 
    !                   for the set of spwfs in storage. If .false., only
    !                   diagonal matrix elements are calculated. 
    !                   Note: in the diagonal basis, we always only calculate
    !                         diagonal matrix elements. 
    !
    !---------------------------------------------------------------------------
    logical, intent(in)        :: fullmatrices
    logical                    :: diag
    integer                    :: k, wave
    real(KIND=dp), allocatable :: temp(:,:)

#if(PASTA == 1)
    ! This is a waste of CPU time for pasta calculations. 
    ! This return is ugly and will require more elegant inclusion later on.
    allocate(temp(1,1)); deallocate(temp)
    return
#endif

    call start_timer(T_spwfangmom)

    allocate(temp(nwt,nwt))
    if(.not.allocated(spwf_J)) then
      allocate(spwf_J(3,nwt,nwt))   ; spwf_J  = 0.0
      allocate(spwf_JTR(3,nwt,nwt)) ; spwf_JTR= 0.0
      allocate(spwf_JTI(3,nwt,nwt)) ; spwf_JTI= 0.0
      allocate(spwf_J2(3,nwt,nwt))  ; spwf_J2 = 0.0
      allocate(spwf_JJ(nwt))        ; spwf_JJ = 0.0

      allocate(spwf_spin(3,nwt,nwt)); spwf_spin= 0.0
      allocate(spwf_STR(3,nwt,nwt)) ; spwf_STR = 0.0
      allocate(spwf_STI(3,nwt,nwt)) ; spwf_STI = 0.0
    endif

    if(.not.allocated(HF_J)) then
     allocate(HF_J(3,nwt))   ;  HF_J = 0.0
     allocate(HF_JTR(3,nwt)) ;  HF_JTR = 0.0
     allocate(HF_JTI(3,nwt)) ;  HF_JTI = 0.0
     allocate(HF_J2(3,nwt))  ;  HF_J2 = 0.0
     allocate(HF_JJ(nwt))    ;  HF_JJ = 0.0

     allocate(HF_spin(3,nwt)); HF_spin= 0.0
     allocate(HF_STR (3,nwt)); HF_STR = 0.0
     allocate(HF_STI (3,nwt)); HF_STI = 0.0
    endif

    if(.not.allocated(can_J)) then
      allocate(can_J(3,nwt))   ; can_J  = 0.0
      allocate(can_JTR(3,nwt)) ; can_JTR= 0.0
      allocate(can_JTI(3,nwt)) ; can_JTI= 0.0
      allocate(can_J2(3,nwt))  ; can_J2 = 0.0
      allocate(can_JJ(nwt))    ; can_JJ = 0.0

      allocate(can_spin(3,nwt)); can_spin = 0.0
      allocate(can_STR(3,nwt)) ; can_STR  = 0.0
      allocate(can_STI(3,nwt)) ; can_STI  = 0.0
    endif

    diag = (.not. fullmatrices)
    ! Operators for which we need no derivatives
    call ME_function(spwf_STR (1,:,:),spin_xt_real,+1,diag,'HF')
    call ME_function(spwf_STI (2,:,:),spin_yt_imag,+1,diag,'HF')
    call ME_function(spwf_spin(3,:,:),spin_z_real, +1,diag,'HF')
    ! Operators for which we need one set of derivatives
    call ME_function_deriv1(spwf_J   (3,:,:),angmom_z_real, +1,diag,'HF')
    call ME_function_deriv1(spwf_JTR (1,:,:),angmom_xt_real,+1,diag,'HF')
    call ME_function_deriv1(spwf_JTI (2,:,:),angmom_yt_imag,+1,diag,'HF')
    ! Operators for which we need two sets of derivatives
    call ME_function_deriv2(spwf_J2  (1,:,:),angmom_x_quad,+1,diag,'HF')
    call ME_function_deriv2(spwf_J2  (2,:,:),angmom_y_quad,+1,diag,'HF')
    call ME_function_deriv2(spwf_J2  (3,:,:),angmom_z_quad,+1,diag,'HF')

    do wave=1,nwt
       spwf_JJ(wave) = (-1. + sqrt(1. + 4*sum(spwf_J2(:,wave,wave))))/2.
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Quantities in the canonical basis 
    if(allocated(CANPSI)) then
      ! Note: we NEVER need the full matrix of angular momenta in the canonical
      ! basis To save on memory, we pass through intermediate arrays.

      can_spin = 0.0d0; can_J   = 0.0d0 ; can_J2 = 0.0d0
      can_STR  = 0.0d0; can_JTR = 0.0d0 
      can_STI  = 0.0d0; can_JTI = 0.0d0
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Operators for which we need no derivatives
      call ME_function(temp,spin_xt_real,+1,.true.,'CAN')
      can_STR(1,:)  = diag_of_mat(temp)
      call ME_function(temp,spin_yt_imag,+1,.true.,'CAN')
      can_STI(2,:)  = diag_of_mat(temp)
      call ME_function(temp,spin_z_real, +1,.true.,'CAN')
      can_spin(3,:) = diag_of_mat(temp)

      ! Operators for which we need one set of derivatives
      call ME_function_deriv1(temp,angmom_z_real, +1,.true.,'CAN')
      can_J   (3,:) = diag_of_mat(temp)
      call ME_function_deriv1(temp,angmom_xt_real,+1,.true.,'CAN')
      can_JTR (1,:) = diag_of_mat(temp)
      call ME_function_deriv1(temp,angmom_yt_imag,+1,.true.,'CAN')
      can_JTI (2,:) = diag_of_mat(temp)

      ! Operators for which we need two sets of derivatives
      call ME_function_deriv2(temp,angmom_x_quad,+1,.true.,'CAN')
      can_J2  (1,:) = diag_of_mat(temp)
      call ME_function_deriv2(temp,angmom_y_quad,+1,.true.,'CAN')
      can_J2  (2,:) = diag_of_mat(temp)
      call ME_function_deriv2(temp,angmom_z_quad,+1,.true.,'CAN')
      can_J2  (3,:) = diag_of_mat(temp)
    endif
    can_JJ = (-1. + sqrt(1. + 4*sum(can_J2,1)))/2.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Quantities in the HF basis 
    if(fullmatrices .and. (.not.diagsphamil)) then
     do k=1,3
        HF_J   (k,:) = transform_mat_diag(spwf_J   (k,:,:), HFtransfo)
        HF_J2  (k,:) = transform_mat_diag(spwf_J2  (k,:,:), HFtransfo)
        HF_JTR (k,:) = transform_mat_diag(spwf_JTR (k,:,:), HFtransfo)
        HF_JTI (k,:) = transform_mat_diag(spwf_JTI (k,:,:), HFtransfo)
        HF_JTI (k,:) = transform_mat_diag(spwf_JTI (k,:,:), HFtransfo)
        HF_spin(k,:) = transform_mat_diag(spwf_spin(k,:,:), HFtransfo)
        HF_STR (k,:) = transform_mat_diag(spwf_STR (k,:,:), HFtransfo)
        HF_STI (k,:) = transform_mat_diag(spwf_STI (k,:,:), HFtransfo)
     enddo
   else
     do k=1,3
       HF_J   (k,:) = diag_of_mat(spwf_J   (k,:,:))
       HF_J2  (k,:) = diag_of_mat(spwf_J2  (k,:,:))
       HF_JTR (k,:) = diag_of_mat(spwf_JTR (k,:,:))
       HF_JTI (k,:) = diag_of_mat(spwf_JTI (k,:,:))
       HF_spin(k,:) = diag_of_mat(spwf_spin(k,:,:))
       HF_STR (k,:) = diag_of_mat(spwf_STR (k,:,:))
       HF_STI (k,:) = diag_of_mat(spwf_STI (k,:,:))
     enddo
   endif
   do wave=1,nwt
     HF_JJ(wave) = (-1. + sqrt(1. + 4*sum(HF_J2(:,wave))))/2.0d0
   enddo
   call stop_timer(T_spwfangmom)

  end subroutine update_spwf_angmom

  function diag_of_mat(A) result(diag)
    !---------------------------------------------------------------------------
    ! Simple function assigning the diagonal matrix elements of a matrix into
    ! a vector.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input: 
    !   A   :  real matrix whose diagonal matrix elements are extracted
    ! Output:
    !   diag:  real vector such that diag(i) = A(i,i)
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: A(nwt,nwt)
    real(KIND=dp)             :: diag(nwt)
    integer                   :: i

    do i=1,nwt
      diag(i) = A(i,i)
    enddo
  end function diag_of_mat
  
  pure function angmom_x_real(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     Re < wf2 | j_x | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_x = 1/2*( 0  1 ) + i z \partial_y - i y\partial_z
    !                 ( 1  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i  
    real(KIND=dp)             :: angmom

    angmom = 0    

    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(  wf2(i,1) * wf1(i,3)                       &
      &                            + wf2(i,2) * wf1(i,4)                       & 
      &                            + wf2(i,3) * wf1(i,1)                       &
      &                            + wf2(i,4) * wf1(i,2))                        
      ! Orbital part
      angmom = angmom &
      &           + wf2(i,2) * meshgrid_shifted(i,3) * dwf1(i,2,1)             &
      &           - wf2(i,2) * meshgrid_shifted(i,2) * dwf1(i,3,1)             &
      !
      &           - wf2(i,1) * meshgrid_shifted(i,3) * dwf1(i,2,2)             &
      &           + wf2(i,1) * meshgrid_shifted(i,2) * dwf1(i,3,2)             &
      !
      &           + wf2(i,4) * meshgrid_shifted(i,3) * dwf1(i,2,3)             &
      &           - wf2(i,4) * meshgrid_shifted(i,2) * dwf1(i,3,3)             &
      !
      &           - wf2(i,3) * meshgrid_shifted(i,3) * dwf1(i,2,4)             &
      &           + wf2(i,3) * meshgrid_shifted(i,2) * dwf1(i,3,4)
    enddo
    angmom = angmom * dv

  end function angmom_x_real

  pure function spin_x_real(wf2, wf1) result(sx)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     Re < wf2 | s_x | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       S_x = 1/2*( 0  1 ) 
    !                 ( 1  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4)
    integer                   :: i  
    real(KIND=dp)             :: sx

    sx = 0    
    do i=1,mv
      ! Spin part
      sx = sx               + 0.5*(  wf2(i,1) * wf1(i,3)                       &
      &                            + wf2(i,2) * wf1(i,4)                       & 
      &                            + wf2(i,3) * wf1(i,1)                       &
      &                            + wf2(i,4) * wf1(i,2))                        
    enddo
    sx = sx * dv

  end function spin_x_real

  pure function angmom_x_imag(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !    Im < wf2 | j_x | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_x = 1/2*( 0  1 ) + i z \partial_y - i y\partial_z
    !                 ( 1  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom

    angmom = 0    

    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(  wf2(i,1) * wf1(i,4)                       &
      &                            - wf2(i,2) * wf1(i,3)                       & 
      &                            + wf2(i,3) * wf1(i,1)                       &
      &                            - wf2(i,4) * wf1(i,2))                        
      ! Orbital part
      angmom = angmom &
      &           + wf2(i,1) * meshgrid_shifted(i,3) * dwf1(i,2,1)             &
      &           - wf2(i,1) * meshgrid_shifted(i,2) * dwf1(i,3,1)             &
      !
      &           + wf2(i,2) * meshgrid_shifted(i,3) * dwf1(i,2,2)             &
      &           - wf2(i,2) * meshgrid_shifted(i,2) * dwf1(i,3,2)             &
      !
      &           + wf2(i,3) * meshgrid_shifted(i,3) * dwf1(i,2,3)             &
      &           - wf2(i,3) * meshgrid_shifted(i,2) * dwf1(i,3,3)             &
      !
      &           + wf2(i,4) * meshgrid_shifted(i,3) * dwf1(i,2,4)             &
      &           - wf2(i,4) * meshgrid_shifted(i,2) * dwf1(i,3,4)
    enddo
    angmom = angmom * dv

  end function angmom_x_imag
  
  pure function spin_x_imag(wf2, wf1) result(sx)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !    Im < wf2 | s_x | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       s_x = 1/2*( 0  1 )
    !                 ( 1  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4)
    integer                   :: i
    real(KIND=dp)             :: sx

    sx = 0    
    do i=1,mv
      ! Spin part
      sx = sx               + 0.5*(  wf2(i,1) * wf1(i,4)                       &
      &                            - wf2(i,2) * wf1(i,3)                       & 
      &                            + wf2(i,3) * wf1(i,1)                       &
      &                            - wf2(i,4) * wf1(i,2))                        
    enddo
    sx = sx * dv

  end function spin_x_imag

  pure function angmom_xt_real(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j_x T | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_x = 1/2*( 0  1 ) + i z \partial_y - i y\partial_z
    !                 ( 1  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom

    angmom = 0    
    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(- wf2(i,1) * wf1(i,1)                       &
      &                            + wf2(i,2) * wf1(i,2)                       & 
      &                            + wf2(i,3) * wf1(i,3)                       &
      &                            - wf2(i,4) * wf1(i,4))                        
      ! Orbital part
      angmom = angmom &
      &           + wf2(i,2) * meshgrid_shifted(i,3) * dwf1(i,2,3)             &
      &           - wf2(i,2) * meshgrid_shifted(i,2) * dwf1(i,3,3)             &
      !
      &           + wf2(i,1) * meshgrid_shifted(i,3) * dwf1(i,2,4)             &
      &           - wf2(i,1) * meshgrid_shifted(i,2) * dwf1(i,3,4)             &
      !
      &           - wf2(i,4) * meshgrid_shifted(i,3) * dwf1(i,2,1)             &
      &           + wf2(i,4) * meshgrid_shifted(i,2) * dwf1(i,3,1)             &
      !
      &           - wf2(i,3) * meshgrid_shifted(i,3) * dwf1(i,2,2)             &
      &           + wf2(i,3) * meshgrid_shifted(i,2) * dwf1(i,3,2)
    enddo
    angmom = angmom * dv
  end function angmom_xt_real
  
  pure function spin_xt_real(wf2, wf1) result(sx)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | s_x T | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       s_x = 1/2*( 0  1 ) 
    !                 ( 1  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4)
    integer                   :: i
    real(KIND=dp)             :: sx

    sx = 0    
    do i=1,mv
      ! Spin part
      sx = sx               + 0.5*(- wf2(i,1) * wf1(i,1)                       &
      &                            + wf2(i,2) * wf1(i,2)                       & 
      &                            + wf2(i,3) * wf1(i,3)                       &
      &                            - wf2(i,4) * wf1(i,4))                        
    enddo
    sx = sx * dv
  end function spin_xt_real

  pure function angmom_x_quad(wf2, dwf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j^\dagger_x j_x | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_x = 1/2*( 0  1 ) + i z \partial_y - i y\partial_z
    !                 ( 1  0 )
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in) :: wf1(mv,4),wf2(mv,4),dwf1(mv,3,4),dwf2(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom, l1,l2,l3,l4, r1,r2,r3,r4

    angmom = 0    
    do i=1,mv

       ! Action of J_x to the right
       r1 = - meshgrid_shifted(i,3)* dwf1(i,2,2) &
       &    + meshgrid_shifted(i,2)* dwf1(i,3,2) & 
       &    + 0.5_dp       * wf1(i,3)
       
       r2 = + meshgrid_shifted(i,3)* dwf1(i,2,1) &
       &    - meshgrid_shifted(i,2)* dwf1(i,3,1) & 
       &    + 0.5_dp       * wf1(i,4)
       
       r3 = - meshgrid_shifted(i,3)* dwf1(i,2,4) &
       &    + meshgrid_shifted(i,2)* dwf1(i,3,4) & 
       &    + 0.5_dp       * wf1(i,1)
       
       r4 = + meshgrid_shifted(i,3)* dwf1(i,2,3) &
       &    - meshgrid_shifted(i,2)* dwf1(i,3,3) & 
       &    + 0.5_dp       * wf1(i,2)
            
       ! Action of J_x to the left
       l1 = - meshgrid_shifted(i,3)* dwf2(i,2,2) &
       &    + meshgrid_shifted(i,2)* dwf2(i,3,2) & 
       &    + 0.5_dp       * wf2(i,3)
       
       l2 = + meshgrid_shifted(i,3)* dwf2(i,2,1) &
       &    - meshgrid_shifted(i,2)* dwf2(i,3,1) & 
       &    + 0.5_dp       * wf2(i,4)
       
       l3 = - meshgrid_shifted(i,3)* dwf2(i,2,4) &
       &    + meshgrid_shifted(i,2)* dwf2(i,3,4) & 
       &    + 0.5_dp       * wf2(i,1)
       
       l4 = + meshgrid_shifted(i,3)* dwf2(i,2,3) &
       &    - meshgrid_shifted(i,2)* dwf2(i,3,3) & 
       &    + 0.5_dp       * wf2(i,2)
            
        angmom = angmom + l1*r1 + l2*r2 + l3*r3 + l4*r4
    enddo
    angmom = angmom * dv
  end function angmom_x_quad

  pure function angmom_y_real(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     Re < wf2 | j_y | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 0 -i ) + i x \partial_z - i z\partial_x
    !                 ( i  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i 
    real(KIND=dp)             :: angmom

    angmom = 0
    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(  wf2(i,1) * wf1(i,4)                       &
      &                            - wf2(i,2) * wf1(i,3)                       & 
      &                            - wf2(i,3) * wf1(i,2)                       &
      &                            + wf2(i,4) * wf1(i,1))                        
      ! Orbital part
      angmom = angmom                                                          &
      &           + wf2(i,2) * meshgrid_shifted(i,1) * dwf1(i,3,1)             &
      &           - wf2(i,2) * meshgrid_shifted(i,3) * dwf1(i,1,1)             &
      !
      &           - wf2(i,1) * meshgrid_shifted(i,1) * dwf1(i,3,2)             &
      &           + wf2(i,1) * meshgrid_shifted(i,3) * dwf1(i,1,2)             &
      !
      &           + wf2(i,4) * meshgrid_shifted(i,1) * dwf1(i,3,3)             &
      &           - wf2(i,4) * meshgrid_shifted(i,3) * dwf1(i,1,3)             &
      !
      &           - wf2(i,3) * meshgrid_shifted(i,1) * dwf1(i,3,4)             &
      &           + wf2(i,3) * meshgrid_shifted(i,3) * dwf1(i,1,4) 
     !--------------------------------------------------------------------------
    enddo
    angmom = angmom * dv

  end function angmom_y_real
  
  pure function spin_y_real(wf2, wf1) result(sy)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     Re < wf2 | s_y | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       s_z = 1/2*( 0 -i ) 
    !                 ( i  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4)
    integer                   :: i 
    real(KIND=dp)             :: sy

    sy = 0
    do i=1,mv
      ! Spin part
      sy = sy               + 0.5*(  wf2(i,1) * wf1(i,4)                       &
      &                            - wf2(i,2) * wf1(i,3)                       & 
      &                            - wf2(i,3) * wf1(i,2)                       &
      &                            + wf2(i,4) * wf1(i,1))                        
    enddo
    sy = sy * dv

  end function spin_y_real

  pure function angmom_y_imag(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     Im < wf2 | j_y | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_y = 1/2*( 0 -i ) + i x \partial_z - i z\partial_x
    !                 ( i  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom

    angmom = 0
    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(- wf2(i,1) * wf1(i,3)                       &
      &                            - wf2(i,2) * wf1(i,4)                       & 
      &                            + wf2(i,3) * wf1(i,1)                       &
      &                            + wf2(i,4) * wf1(i,2))                        
      ! Orbital part
      angmom = angmom &
      &           - wf2(i,1) * meshgrid_shifted(i,3) * dwf1(i,1,1)             &
      &           + wf2(i,1) * meshgrid_shifted(i,1) * dwf1(i,3,1)             &
      !
      &           - wf2(i,2) * meshgrid_shifted(i,3) * dwf1(i,1,2)             &
      &           + wf2(i,2) * meshgrid_shifted(i,1) * dwf1(i,3,2)             &
      !
      &           - wf2(i,3) * meshgrid_shifted(i,3) * dwf1(i,1,3)             &
      &           + wf2(i,3) * meshgrid_shifted(i,1) * dwf1(i,3,3)             &
      !
      &           - wf2(i,4) * meshgrid_shifted(i,3) * dwf1(i,1,4)             &
      &           + wf2(i,4) * meshgrid_shifted(i,1) * dwf1(i,3,4)
    enddo
    angmom = angmom * dv

  end function angmom_y_imag
  
  pure function spin_y_imag(wf2, wf1) result(sy)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     Im < wf2 | s_y | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       s_y = 1/2*( 0 -i ) 
    !                 ( i  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4)
    integer                   :: i
    real(KIND=dp)             :: sy

    sy = 0
    do i=1,mv
      ! Spin part
      sy = sy               + 0.5*(- wf2(i,1) * wf1(i,3)                       &
      &                            - wf2(i,2) * wf1(i,4)                       & 
      &                            + wf2(i,3) * wf1(i,1)                       &
      &                            + wf2(i,4) * wf1(i,2))                        
    enddo
    sy = sy * dv

  end function spin_y_imag

  pure function angmom_yt_imag(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j_yT | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_y = 1/2*( 0 -i ) + i x \partial_z - i z\partial_x
    !                 ( i  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom

    angmom = 0
    do i=1,mv
      ! Spin part
      angmom = angmom       + 0.5*(+ wf2(i,1) * wf1(i,1)                       &
      &                            - wf2(i,2) * wf1(i,2)                       & 
      &                            + wf2(i,3) * wf1(i,3)                       &
      &                            - wf2(i,4) * wf1(i,4))                        
      ! Orbital part
      angmom = angmom                                                          &
      &           - wf2(i,1) * meshgrid_shifted(i,3) * dwf1(i,1,3)             &
      &           + wf2(i,1) * meshgrid_shifted(i,1) * dwf1(i,3,3)             &
      !
      &           + wf2(i,2) * meshgrid_shifted(i,3) * dwf1(i,1,4)             &
      &           - wf2(i,2) * meshgrid_shifted(i,1) * dwf1(i,3,4)             &
      !
      &           + wf2(i,3) * meshgrid_shifted(i,3) * dwf1(i,1,1)             &
      &           - wf2(i,3) * meshgrid_shifted(i,1) * dwf1(i,3,1)             &
      !
      &           - wf2(i,4) * meshgrid_shifted(i,3) * dwf1(i,1,2)             &
      &           + wf2(i,4) * meshgrid_shifted(i,1) * dwf1(i,3,2) 
     !-------------------------------------------------------------------------
    enddo
    angmom = angmom * dv

  end function angmom_yt_imag
  
  pure function spin_yt_imag(wf2, wf1) result(sy)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | s_yT | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       s_y = 1/2*( 0 -i )
    !                 ( i  0 )
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4)
    integer                   :: i
    real(KIND=dp)             :: sy

    sy = 0
    do i=1,mv
      ! Spin part
      sy = sy               + 0.5*(+ wf2(i,1) * wf1(i,1)                       &
      &                            - wf2(i,2) * wf1(i,2)                       & 
      &                            + wf2(i,3) * wf1(i,3)                       &
      &                            - wf2(i,4) * wf1(i,4))                        
     !-------------------------------------------------------------------------
    enddo
    sy = sy * dv

  end function spin_yt_imag

  pure function angmom_y_quad(wf2, dwf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j^\dagger_y j_y | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_y = 1/2*( 0 -i ) + i x \partial_z - i z\partial_x
    !                 ( i  0 )    
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4),dwf2(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom, l1,l2,l3,l4, r1,r2,r3,r4

    angmom = 0    
    do i=1,mv
       ! Action of J_y to the right
       r1 = - meshgrid_shifted(i,1)* dwf1(i,3,2) &
       &    + meshgrid_shifted(i,3)* dwf1(i,1,2) & 
       &    + 0.5_dp       * wf1(i,4)
       
       r2 = + meshgrid_shifted(i,1)* dwf1(i,3,1) &
       &    - meshgrid_shifted(i,3)* dwf1(i,1,1) & 
       &    - 0.5_dp       * wf1(i,3)
       
       r3 = - meshgrid_shifted(i,1)* dwf1(i,3,4) &
       &    + meshgrid_shifted(i,3)* dwf1(i,1,4) & 
       &    - 0.5_dp       * wf1(i,2)
       
       r4 = + meshgrid_shifted(i,1)* dwf1(i,3,3) &
       &    - meshgrid_shifted(i,3)* dwf1(i,1,3) & 
       &    + 0.5_dp       * wf1(i,1)
            
       ! Action of J_y to the left
       l1 = - meshgrid_shifted(i,1)* dwf2(i,3,2) &
       &    + meshgrid_shifted(i,3)* dwf2(i,1,2) & 
       &    + 0.5_dp       * wf2(i,4)
       
       l2 = + meshgrid_shifted(i,1)* dwf2(i,3,1) &
       &    - meshgrid_shifted(i,3)* dwf2(i,1,1) & 
       &    - 0.5_dp       * wf2(i,3)
       
       l3 = - meshgrid_shifted(i,1)* dwf2(i,3,4) &
       &    + meshgrid_shifted(i,3)* dwf2(i,1,4) & 
       &    - 0.5_dp       * wf2(i,2)
       
       l4 = + meshgrid_shifted(i,1)* dwf2(i,3,3) &
       &    - meshgrid_shifted(i,3)* dwf2(i,1,3) & 
       &    + 0.5_dp       * wf2(i,1)
       
        angmom = angmom + l1*r1 + l2*r2 + l3*r3 + l4*r4
    enddo
    angmom = angmom * dv
  end function angmom_y_quad

  pure function angmom_z_real(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j_z | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 1  0 ) + i y \partial_x - i x\partial_y
    !                 ( 0 -1 ) 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    real(KIND=dp)             :: angmom
    integer                   :: i

    angmom = 0

    do i=1,nx*ny*nz
      !-------------------------------------------------------------------------
      ! Real part  
      ! Spin part
      angmom = angmom    + 0.5*(         wf2(i,1) * wf1(i,1)                   &
      &                                + wf2(i,2) * wf1(i,2)                   & 
      &                                - wf2(i,3) * wf1(i,3)                   &
      &                                - wf2(i,4) * wf1(i,4))                        
      ! Orbital part
      angmom = angmom              +                                           &
      &                      wf2(i,2) * meshgrid_shifted(i,2) * dwf1(i,1,1)    &
      &                    - wf2(i,2) * meshgrid_shifted(i,1) * dwf1(i,2,1)    &
      !                 
      &                    - wf2(i,1) * meshgrid_shifted(i,2) * dwf1(i,1,2)    &
      &                    + wf2(i,1) * meshgrid_shifted(i,1) * dwf1(i,2,2)    &
      !
      &                    + wf2(i,4) * meshgrid_shifted(i,2) * dwf1(i,1,3)    &
      &                    - wf2(i,4) * meshgrid_shifted(i,1) * dwf1(i,2,3)    &
      !                 
      &                    - wf2(i,3) * meshgrid_shifted(i,2) * dwf1(i,1,4)    &
      &                    + wf2(i,3) * meshgrid_shifted(i,1) * dwf1(i,2,4)
    enddo

    angmom = angmom * dv

  end function angmom_z_real
  
  pure function spin_z_real(wf2, wf1) result(sz)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | s_z | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       s_z = 1/2*( 1  0 )
    !                 ( 0 -1 ) 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4)
    real(KIND=dp)             :: sz
    integer                   :: i

    sz = 0
    do i=1,nx*ny*nz
      !-------------------------------------------------------------------------
      ! Real part  
      ! Spin part
      sz = sz            + 0.5*(         wf2(i,1) * wf1(i,1)                   &
      &                                + wf2(i,2) * wf1(i,2)                   & 
      &                                - wf2(i,3) * wf1(i,3)                   &
      &                                - wf2(i,4) * wf1(i,4))                        
    enddo
    sz = sz * dv

  end function spin_z_real

  pure function angmom_z_imag(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !    Im < wf2 | j_z | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 1  0 ) + i y \partial_x - i x\partial_y
    !                 ( 0 -1 ) 
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom

    angmom = 0    

    do i=1,mv
      !-------------------------------------------------------------------------
      ! Real part  

      ! Spin part
      angmom = angmom       + 0.5*(- wf2(i,2) * wf1(i,1)                       &
      &                            + wf2(i,1) * wf1(i,2)                       & 
      &                            + wf2(i,4) * wf1(i,3)                       &
      &                            - wf2(i,3) * wf1(i,4))                        
      ! Orbital part
      angmom = angmom &
      &           + wf2(i,1) * meshgrid_shifted(i,2) * dwf1(i,1,1)             &
      &           - wf2(i,1) * meshgrid_shifted(i,1) * dwf1(i,2,1)             &
      !
      &           + wf2(i,2) * meshgrid_shifted(i,2) * dwf1(i,1,2)             &
      &           - wf2(i,2) * meshgrid_shifted(i,1) * dwf1(i,2,2)             &
      !
      &           + wf2(i,3) * meshgrid_shifted(i,2) * dwf1(i,1,3)             &
      &           - wf2(i,3) * meshgrid_shifted(i,1) * dwf1(i,2,3)             &
      !
      &           + wf2(i,4) * meshgrid_shifted(i,2) * dwf1(i,1,4)             &
      &           - wf2(i,4) * meshgrid_shifted(i,1) * dwf1(i,2,4)
    enddo
    angmom = angmom * dv

  end function angmom_z_imag
  
  pure function spin_z_imag(wf2, wf1) result(sz)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !    Im < wf2 | s_z | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       s_z = 1/2*( 1  0 ) 
    !                 ( 0 -1 ) 
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4)
    integer                   :: i
    real(KIND=dp)             :: sz

    sz = 0    
    do i=1,mv
      !-------------------------------------------------------------------------
      ! Spin part
      sz = sz               + 0.5*(- wf2(i,2) * wf1(i,1)                       &
      &                            + wf2(i,1) * wf1(i,2)                       & 
      &                            + wf2(i,4) * wf1(i,3)                       &
      &                            - wf2(i,3) * wf1(i,4))                        
    enddo
    sz = sz * dv

  end function spin_z_imag

  pure function angmom_zt_real(wf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j_z T | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 1  0 ) + i y \partial_x - i x\partial_y
    !                 ( 0 -1 ) 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4)
    integer                   :: i 
    real(KIND=dp)             :: angmom

    angmom = 0
    do i=1,mv
      !-------------------------------------------------------------------------
      ! Real part  

      ! Spin part
      angmom = angmom    + 0.5*(         wf2(i,1) * wf1(i,3)                   &
      &                                - wf2(i,2) * wf1(i,4)                   & 
      &                                + wf2(i,3) * wf1(i,1)                   &
      &                                - wf2(i,4) * wf1(i,2))                        
      ! Orbital part
      angmom = angmom                                                          &
      &                    + wf2(i,2) * meshgrid_shifted(i,2) * dwf1(i,1,3)    &
      &                    - wf2(i,2) * meshgrid_shifted(i,1) * dwf1(i,2,3)    &
      !                 
      &                    + wf2(i,1) * meshgrid_shifted(i,2) * dwf1(i,1,4)    &
      &                    - wf2(i,1) * meshgrid_shifted(i,1) * dwf1(i,2,4)    &
      !
      &                    - wf2(i,4) * meshgrid_shifted(i,2) * dwf1(i,1,1)    &
      &                    + wf2(i,4) * meshgrid_shifted(i,1) * dwf1(i,2,1)    &
      !                 
      &                    - wf2(i,3) * meshgrid_shifted(i,2) * dwf1(i,1,2)    &
      &                    + wf2(i,3) * meshgrid_shifted(i,1) * dwf1(i,2,2)
    enddo
    angmom = dv * angmom
  end function angmom_zt_real
  
  pure function angmom_z_quad(wf2, dwf2, wf1, dwf1) result(angmom)
    !---------------------------------------------------------------------------
    ! calculate the matrix element
    !
    !     < wf2 | j^\dagger_z j_z | wf1 >
    !
    ! Note that this routine assumes that all checks regarding symmetries have 
    ! been performed.
    !
    !       J_z = 1/2*( 1  0 ) + i y \partial_x - i x\partial_y
    !                 ( 0 -1 ) 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: wf1(mv,4), wf2(mv,4), dwf1(mv,3,4),dwf2(mv,3,4)
    integer                   :: i
    real(KIND=dp)             :: angmom, l1,l2,l3,l4, r1,r2,r3,r4

    angmom = 0    
    do i=1,mv
       ! Action of J_z to the right
       r1 = - meshgrid_shifted(i,2)* dwf1(i,1,2) &
       &    + meshgrid_shifted(i,1)* dwf1(i,2,2) & 
       &    + 0.5_dp       * wf1(i,1)
       
       r2 = + meshgrid_shifted(i,2)* dwf1(i,1,1) &
       &    - meshgrid_shifted(i,1)* dwf1(i,2,1) & 
       &    + 0.5_dp       * wf1(i,2)
       
       r3 = - meshgrid_shifted(i,2)* dwf1(i,1,4) &
       &    + meshgrid_shifted(i,1)* dwf1(i,2,4) & 
       &    - 0.5_dp       * wf1(i,3)
       
       r4 = + meshgrid_shifted(i,2)* dwf1(i,1,3) &
       &    - meshgrid_shifted(i,1)* dwf1(i,2,3) & 
       &    - 0.5_dp       * wf1(i,4)
            
       ! Action of J_z to the left
       l1 = - meshgrid_shifted(i,2)* dwf2(i,1,2) &
       &    + meshgrid_shifted(i,1)* dwf2(i,2,2) & 
       &    + 0.5_dp       * wf2(i,1)
       
       l2 = + meshgrid_shifted(i,2)* dwf2(i,1,1) &
       &    - meshgrid_shifted(i,1)* dwf2(i,2,1) & 
       &    + 0.5_dp       * wf2(i,2)
       
       l3 = - meshgrid_shifted(i,2)* dwf2(i,1,4) &
       &    + meshgrid_shifted(i,1)* dwf2(i,2,4) & 
       &    - 0.5_dp       * wf2(i,3)
       
       l4 = + meshgrid_shifted(i,2)* dwf2(i,1,3) &
       &    - meshgrid_shifted(i,1)* dwf2(i,2,3) & 
       &    - 0.5_dp       * wf2(i,4)
       
       angmom = angmom + l1*r1 + l2*r2 + l3*r3 + l4*r4
    enddo
    angmom = angmom * dv
  end function angmom_z_quad
  
  function AngMomOperator(psi, dpsi, direction) result(Jpsi)
    !---------------------------------------------------------------------------
    ! Calculate the action of the angular momentum operator on a spwf.
    !
    !       J_x = 1/2*( 0  1 ) + i z \partial_y - i y\partial_z
    !                 ( 1  0 )
    !
    !       J_y = 1/2*( 0 -i ) + i x \partial_z - i z\partial_x
    !                 ( i  0 )
    !
    !       J_z = 1/2*( 1  0 ) + i y \partial_x - i x\partial_y
    !                 ( 0 -1 )
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     psi   : spwf values on the mesh (= spinor)
    !    dpsi   : derivative of the spwf on the mesh ( = 3 spinors)
    !  direction: Cartesian direction, indicates which angular momentum to 
    !             calculate; 1/2/3 = x/y/z
    ! Output:
    !   Jpsi    : the action of the angular momentum on the spwf ( = spinor)
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: psi(mv,4), dpsi(mv,3,4)
    integer, intent(in)       :: direction
    real(KIND=dp), allocatable :: Jpsi(:,:)
     
    JPsi = Orbital(dpsi, direction) + Pauli(psi,direction)   
 
  end function AngMomOperator
  
  pure function Orbital(dpsi, direction) result(Lpsi)
    !---------------------------------------------------------------------------
    ! Calculate the action of the orbital momentum operator on a spwf.
    !
    !       L_x = i z \partial_y - i y\partial_z
    !       L_y = i x \partial_z - i z\partial_x
    !       L_z = i y \partial_x - i x\partial_y
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     psi   : spwf values on the mesh (= spinor)
    !    dpsi   : derivative of the spwf on the mesh ( = 3 spinors)
    !  direction: Cartesian direction, indicates which angular momentum to 
    !             calculate; 1/2/3 = x/y/z
    ! Output:
    !   Lpsi    : the action of the orbital momentum on the spwf ( = spinor)
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: dpsi(mv,3,4)
    integer, intent(in)        :: direction
    real(KIND=dp), allocatable :: Lpsi(:,:), temp(:,:)
    integer                    :: k
    
    allocate(LPsi(mv,4)) ; LPsi = 0.0d0
    allocate(temp(mv,4)) ; temp = 0.0d0
    
    ! Calculate the spatial operation
    select case(Direction)
    case(1)
      do k=1,4
        temp(:,k) =  meshgrid_shifted(:,3) *dpsi(:,2,k)                        &
        &         -  meshgrid_shifted(:,2) *dpsi(:,3,k)
      enddo
    case(2)
      do k=1,4
        temp(:,k) =  meshgrid_shifted(:,1) *dpsi(:,3,k)                        &
        &         -  meshgrid_shifted(:,3) *dpsi(:,1,k)
      enddo
    case(3)
      do k=1,4
        temp(:,k) =  meshgrid_shifted(:,2) *dpsi(:,1,k)                        &
        &         -  meshgrid_shifted(:,1) *dpsi(:,2,k)
      enddo
    end select
    
    ! And multiply by 'i'
    LPsi(:,1) = - temp(:,2)
    LPsi(:,2) = + temp(:,1)
    LPsi(:,3) = - temp(:,4)
    LPsi(:,4) = + temp(:,3)

  end function Orbital
  
  pure function Pauli(psi, direction) result(Spsi)
    !---------------------------------------------------------------------------
    ! Calculate the action of the spin operator on a spwf.
    !
    !       S_x = 1/2*( 0  1 ) 
    !                 ( 1  0 )
    !
    !       S_y = 1/2*( 0 -i ) 
    !                 ( i  0 )
    !
    !       S_z = 1/2*( 1  0 ) 
    !                 ( 0 -1 )
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     psi   : spwf values on the mesh (= spinor)
    !  direction: Cartesian direction, indicates which angular momentum to 
    !             calculate; 1/2/3 = x/y/z
    ! Output:
    !   Spsi    : the action of the spin operator on the spwf ( = spinor)
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: psi(mv,4)
    integer, intent(in)        :: direction
    real(KIND=dp), allocatable :: Spsi(:,:)
    integer                    :: i
    
    allocate(SPsi(mv,4)) ; SPsi = 0.0d0
    
    if(Direction.eq.1) then
        !\sigma_x = ( 0  1 )
        !           ( 1  0 )
        do i=1,mv
          SPsi(i,1) = Psi(i,3)
          SPsi(i,2) = Psi(i,4)
          SPsi(i,3) = Psi(i,1)
          SPsi(i,4) = Psi(i,2)   
        enddo      
    elseif(Direction.eq.2) then
        !\sigma_y = ( 0 -i )
        !           ( i  0 )
        do i=1,mv          
          SPsi(i,1) =   Psi(i,4)
          SPsi(i,2) = - Psi(i,3)
          SPsi(i,3) = - Psi(i,2)
          SPsi(i,4) =   Psi(i,1)
        enddo       
    elseif(Direction.eq.3) then
        !\sigma_z = ( 1  0 )
        !           ( 0 -1 )
        do i=1,mv      
          SPsi(i,1) =   Psi(i,1)
          SPsi(i,2) =   Psi(i,2)
          SPsi(i,3) = - Psi(i,3)
          SPsi(i,4) = - Psi(i,4)
        enddo
    endif
  
  end function Pauli
  
  pure function ImagMultiplySpinor(Psi, Phi) result(ImPsiPhi)
    !-------------------------------------------------------------------------
    ! Computes the imaginary part of 
    !       Psi^{dagger} Phi
    !-------------------------------------------------------------------------
    real(KIND=dp),intent(in) :: Psi(nx*ny*nz,4), Phi(nx*ny*nz,4)
    real(KIND=dp)            :: ImPsiPhi(nx*ny*nz)
    integer                  :: i

    do i=1,nx*ny*nz
      ImPsiPhi(i) = Psi(i,1) * Phi(i,2) &
      &           - Psi(i,2) * Phi(i,1) &
      &           + Psi(i,3) * Phi(i,4) &
      &           - Psi(i,4) * Phi(i,3)
    enddo

  end function ImagMultiplySpinor

!===============================================================================
! Routines for calculating diverse properties of the spwfs
!===============================================================================

  subroutine update_spwf_properties( fullmatrices )
      !-------------------------------------------------------------------------
      ! Wrapper function to update all spwf information that needs to be
      ! recalculated. This is not hidden inside some other routine, simply 
      ! because the timing of this call is important: it needs to be AFTER
      ! the construction of the canonical basis.
      ! 
      ! Input:
      !    fullmatrices : if .true., force calculation in HF and canonical basis
      !                   even if diagsphamil = .false.
      !
      ! All of these calculations can be trivially executed in any basis which
      ! is explicitly stored. For the canonical basis hence, this is trivial in
      ! every runmode of the calculation. For the HFbasis, this is only trivial
      ! if diagsphamil = .true.. If diagsphamil is .false., then we can still
      ! calculate everything using the HF-transformation and a full set of 
      ! matrix elements. Since the latter are expensive to calculate, and 
      ! expectation values of operators in the HF-basis are not so relevant 
      ! to a HFB calculation (except for printing) this routine offers the 
      ! option to skip the expensive calculation by setting fullmatrices=.false.
      ! Ofcourse, this means that HFbasis values should not be trusted....
      !-------------------------------------------------------------------------

      logical, intent(in) :: fullmatrices

      call update_spwf_symmetries(fullmatrices) ! <symmetry operators>
      call update_spwf_angmom(fullmatrices)     ! angular momentum
      call update_spwf_r2(fullmatrices)         ! <r^2> 
  end subroutine update_spwf_properties
  
  subroutine update_spwf_r2(fullmatrices)
      !-------------------------------------------------------------------------
      ! Calculate the single-particle expectation <r^2> for every spwf in 
      ! the Hartree-Fock and canonical basis.
      !
      ! Input:
      !   fullmatrices : if .true., force calculation in the HF and canonical
      !                  basis even if diagsphamil = .false.
      !-------------------------------------------------------------------------
      logical, intent(in) :: fullmatrices

      real(KIND=dp), pointer     :: rme(:,:)
      real(KIND=dp)              :: r2(mv)
!      integer                    :: B, N, si

#if(PASTA == 1)
      ! This is a waste of CPU time for pasta calculations. 
      ! This return is ugly and will require more elegant inclusion later on.
      return
#endif

      ! Value of r^2 = X^2 + Y^2 + Z^2 on the mesh
      r2 = sum(meshgrid,2)**2

      if(diagsphamil) then
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! This is easy: both HFBasis and canbasis are explicitly stored
        call ME_scalar(spwf_r2_HF, r2, .true., 'HF')
        if(allocated(canpsi)) then
          call ME_scalar(spwf_r2_can, r2, .true., 'CAN')
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else
       if(.not.allocated(spwf_r2_HF )) allocate(spwf_r2_HF(nwt,nwt)) ; spwf_r2_HF = 0.0d0
       if(.not.allocated(spwf_r2_can)) allocate(spwf_r2_can(nwt,nwt)) ; spwf_r2_HF = 0.0d0
!        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!        ! The canonical basis is still trivial, but is now stored in HFPSI
!        ! Note: it is safe to assume the calculation is a HFB one; this is the
!        !       only case when diagsphamil should be set to false.
!        call ME_scalar(spwf_r2_can, r2, .true., 'HF')

!        if(fullmatrices) then
!          ! For the HFbasis, things are more involved.....
!          ! a) calculate the entire matrix of r2
!          call ME_scalar(spwf_r2_HF, r2, .false., 'HF')

!          ! b) transform the matrix elements to the HF basis
!          si = 0
!          do B=1,8  
!            N = HFBlocks_global(B) ! <---- This loop is over global spwf indices
!            rme => spwf_r2_hf(si+1:si+N, si+1:si+N)
!            ! .... and then transform to the real Hartree-Fock basis
!            rme = matmul(transpose(HFtransfo(si+1:si+N, si+1:si+N)), rme)
!            rme = matmul(            rme,HFtransfo(si+1:si+N, si+1:si+N))
!            si = si + N
!          enddo
!        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      endif
      
  end subroutine update_spwf_r2
  
  function spwf_parities(basis, fullmatrices) result(P)
      !-------------------------------------------------------------------------
      ! Calculation of the single-particle matrix elements of parity P.
      !
      !         < i | P | j >
      ! 
      ! Input:
      !     basis        : set of single-particle wavefunctions to perform the 
      !                    calculation for.
      !     fullmatrices : whether to calculate all matrix elements (.true.)
      !                    or only the diagonal ones (.false.)
      ! Output:
      !     P    : set of parities
      !
      ! Currently, this routine is somewhat hardcoded for the symmetry options
      ! corresponding to EV8/CR8/EV4; in particular it assumes the conservation
      ! of z-signature. 
      ! 
      ! In time, Hephaestos should be able to deal more properly with this.
      !-------------------------------------------------------------------------
      real(KIND=dp), allocatable         :: P(:,:)
      real(KIND=dp), intent(in), target  :: basis(nx*ny*nz,4,nwt)
      logical, intent(in)                :: fullmatrices
      integer                            :: wave 
      real(KIND=dp)                      :: trash
      
$PBROKEN      integer :: B, N, i, j,k, si, wave2, startind, endind
$PBROKEN      real(KIND=dp), pointer             :: left1(:,:,:), right1(:,:,:)
$PBROKEN      real(KIND=dp), pointer             :: left2(:,:,:), right2(:,:,:)
$PBROKEN      real(KIND=dp), pointer             :: left3(:,:,:), right3(:,:,:)
$PBROKEN      real(KIND=dp), pointer             :: left4(:,:,:), right4(:,:,:)
        
      ! A statement to stop the compiler complaining about unused variables
      if(fullmatrices) trash = basis(1,1,1)

#if(PASTA == 1)
    ! This is a waste of CPU time for pasta calculations. 
    ! This return is ugly and will require more elegant inclusion later on.
    return
#endif

      allocate(P(nwt,nwt))

$PCON      P = 0 
$PCON      do wave=1,nwt
$PCON         if    (wave .le. sum(HFBlocks_global(1:2))) then
$PCON               P(wave,wave) = +1
$PCON         elseif(wave .le. sum(HFBlocks_global(1:4))) then
$PCON               P(wave,wave) = -1
$PCON         elseif(wave .le. sum(HFBlocks_global(1:6))) then
$PCON               P(wave,wave) = +1
$PCON         else  
$PCON               P(wave,wave) = -1
$PCON         endif
$PCON      enddo 
     
$PBROKEN   si = 0
$PBROKEN   do B=1,8
$PBROKEN      N = HFBlocks(B) ; if(N.eq.0) cycle
$PBROKEN        
$PBROKEN      do wave = si+1, si+N
$PBROKEN        ! This is ugly, because the FORTRAN standard does not allow
$PBROKEN        ! for sufficiently general pointer remapping....
$PBROKEN        left1(1:nx,1:ny,1:nz) => Basis(1:nx*ny*nz,1,wave)
$PBROKEN        left2(1:nx,1:ny,1:nz) => Basis(1:nx*ny*nz,2,wave)
$PBROKEN        left3(1:nx,1:ny,1:nz) => Basis(1:nx*ny*nz,3,wave)
$PBROKEN        left4(1:nx,1:ny,1:nz) => Basis(1:nx*ny*nz,4,wave)
$PBROKEN
$PBROKEN        startind   = si + wave
$PBROKEN        if(fullmatrices) then
$PBROKEN          endind   = si + N
$PBROKEN        else
$PBROKEN          endind   = wave
$PBROKEN        endif
$PBROKEN        do wave2 = wave,endind
$PBROKEN
$PBROKEN         right1(1:nx,1:ny,1:nz) => Basis(1:nx*ny*nz,1,wave2)
$PBROKEN         right2(1:nx,1:ny,1:nz) => Basis(1:nx*ny*nz,2,wave2)
$PBROKEN         right3(1:nx,1:ny,1:nz) => Basis(1:nx*ny*nz,3,wave2)
$PBROKEN         right4(1:nx,1:ny,1:nz) => Basis(1:nx*ny*nz,4,wave2)
$PBROKEN
$PBROKEN         P(wave,wave2) = 0
$PBROKEN         do k=1,nz
$PBROKEN          do j=1,ny
$PBROKEN            do i=1,nx
$PBROKEN             P(wave,wave2)=P(wave,wave2)+left1(i,j,k)*right1(i,j,nz-k+1)
$PBROKEN             P(wave,wave2)=P(wave,wave2)+left2(i,j,k)*right2(i,j,nz-k+1)
$PBROKEN             P(wave,wave2)=P(wave,wave2)-left3(i,j,k)*right3(i,j,nz-k+1)
$PBROKEN             P(wave,wave2)=P(wave,wave2)-left4(i,j,k)*right4(i,j,nz-k+1)
$PBROKEN            enddo
$PBROKEN          enddo
$PBROKEN         enddo
$PBROKEN         if(mod(B,2) .eq. 0) P(wave,wave2) = - P(wave,wave2)
$PBROKEN         P(wave ,wave2) = P(wave,wave2) * dv
$PBROKEN         P(wave2,wave ) = P(wave,wave2) 
$PBROKEN        enddo
$PBROKEN      enddo
$PBROKEN      si = si + N
$PBROKEN   enddo
  
  end function spwf_parities

!===============================================================================
! Routines useful to simplify the parallelization of the calculation of 
! complete-matrices of the expectation values of single-particle operators in
! the case of MPI calculations.
!===============================================================================

subroutine ME_scalar(ME, oper, diag, basis)
  !-----------------------------------------------------------------------------
  ! Evaluate the matrix elements of a simple position-dependent operator 
  ! for all spwfs in memory. This routine is valid for operators whose only 
  ! non-zero matrix elements are to be found within symmetry blocks.
  ! - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -
  ! Input: 
  !     oper : value of an operator on all mesh points
  !     diag : if .true., only compute the diagonal matrix elements
  ! Output:
  !     ME   : matrix elements such that
  !               M(i,j) = <i|oper|j> = int d^3 psi_i^*(r) oper(r) psi_j(r)
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Important note: in the case of MPI calculations, this implementation is
  !                 naive and requires quite a bunch of communications between
  !                 ranks. For now, the load balancing is made simply on 
  !                 protons vs. neutrons.
  !-----------------------------------------------------------------------------
  
  real(KIND=dp), intent(out), allocatable :: ME(:,:)
  real(KIND=dp), intent(in)   :: oper(mv)
  character(len=*), intent(in):: basis
  logical, intent(in)         :: diag

  integer                :: B, si, N
  integer                :: wave, wave_global, wave2, wave2_global
  integer                :: designated_rank(8), calc_rank, ranki, rankj
  real(KIND=dp), pointer :: psi(:,:,:)
  real(KIND=dp)          :: psi_i(mv,4), psi_j(mv,4)
#if(USE_MPI>0)
  integer :: mpi_err
#endif

  if(.not.allocated(ME)) allocate(ME(nwt,nwt))
  ME = 0.0d0 ! clearly zero everything for the allreduce call later

  if(diag) then
    ! We calculate only diagonal matrix elements; we thus need no additional
    ! communications as each rank can just do its own calculation.
    if(to_upper(adjustl(basis)) .eq. 'HF') then
      psi => HFpsi
    elseif(to_upper(adjustl(basis)) .eq. 'CAN') then
      psi => canpsi
    elseif(to_upper(adjustl(basis)) .eq. 'DEN') then
      psi => denpsi
    endif

    do wave=1,nwt_local
      wave_global  = spwf_map(wave)

      ME(wave_global,wave_global) = dv * sum(oper*sum(psi(:,:,wave)**2,2))
    enddo
  else
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! We have to calculate the full matrix
    do B=1,8
      ! First, we designate a rank to do the calculations
      !    = the first rank storing spwfs in a symmetry block
      designated_rank(B) = rank_map(sum(HFBlocks_global(1:B-1))+1)
    enddo

    si = 0
    do B=1,8
      N = HFBlocks_global(B)  ! <= this loops over global spwfs
      if(N.eq.0) cycle
      calc_rank = designated_rank(B) !  This is the rank doing the integrations

      do wave_global=si+1,si+N
        ! Find the rank storing psi_i and its local index
        ranki = rank_map(wave_global)
        wave  = spwf_inverse(wave_global)
        ! transfer psi_i to the calculating rank
#if(USE_MPI>0)
        call Transfer_psi(psi_i, wave, basis,ranki, calc_rank)
#else
        call Transfer_psi(psi_i, wave, basis)
#endif
        do wave2_global=wave_global,si+N
          ! Find the rank storing psi_j and its local index
          rankj = rank_map(wave2_global)
          wave2 = spwf_inverse(wave2_global)
          ! transfer psi_j to the calculating rank
#if(USE_MPI>0)
          call Transfer_psi(psi_j, wave2, basis, rankj, calc_rank)
#else
          call Transfer_psi(psi_j, wave2, basis)
#endif
          if(MPI_RANK.eq.calc_rank) then
            ! Perform the integration.....
            ME(wave_global,wave2_global) = dv*sum(oper*sum(psi_i*psi_j,2))
            ! this type of matrix elements are always symmetric ....
            ME(wave2_global, wave_global) = ME(wave_global , wave2_global)
          endif
        enddo
      enddo
      si = si + N
    enddo
  endif
#if(USE_MPI>0) 
  ! Transferring results to all ranks
  call MPI_ALLREDUCE(MPI_IN_PLACE, ME, nwt**2, MPI_REAL8, MPI_SUM, &
  &                                                      MPI_COMM_WORLD,mpi_err)
#endif 

end subroutine ME_scalar

subroutine ME_function(ME, f, sym, diag, basis)
  !-----------------------------------------------------------------------------
  ! Evaluate the matrix elements of an operator on the spwfs using a function  
  ! defined for its evaluation between any two. I.e. given 
  !
  !      f(psi_i, psi_j) = < psi_i | O | psi_j >
  !
  ! this function evaluates either the full matrix M_ij = f(psi_i, psi_j) or
  ! just its diagonal elements. 
  !
  ! Assumptions:
  ! 1) A symmetry with sign sym, meaning that 
  !      f(psi_j, psi_i) = sym * f(psi_i , psi_j)
  ! 2) ALL matrix elements are real, hence f should return real values.
  ! 3) No derivatives are involved in f, i.e. that procedure takes as arguemts
  !    only the values of the spwfs; f(psi_i, psi_j)
  ! 4) the matrix elements are only calculated within symmetry blocks. 
  ! - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -
  ! Input: 
  !     f    : function procedure, with the interface of spin_z_real
  !     sym  : 
  !     diag : if .true., only compute the diagonal matrix elements
  ! Output:
  !     ME   : matrix elements such that
  !               M(i,j) = <i|oper|j> = int d^3 psi_i^*(r) oper(r) psi_j(r)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Important note
  !    in the case of MPI calculations, this implementation is naive and 
  !    requires quite a bunch of communications between ranks. For now, the 
  !    load balancing is made simply on protons vs. neutrons.
  !-----------------------------------------------------------------------------
  
  real(KIND=dp), intent(out)  :: ME(:,:)
  procedure (spin_z_real)     :: f
  integer, intent(in)         :: sym
  character(len=*), intent(in):: basis
  logical, intent(in)         :: diag

  integer                :: B, si, N
  integer                :: wave, wave_global, wave2, wave2_global
  integer                :: designated_rank(8), calc_rank, ranki, rankj
  real(KIND=dp), pointer :: psi(:,:,:)
  real(KIND=dp)          :: psi_i(mv,4), psi_j(mv,4)
#if(USE_MPI>0)
  integer :: mpi_err
#endif

  ME = 0.0d0 ! clearly zero everything for the allreduce call later

  if(diag) then
    ! We calculate only diagonal matrix elements; we thus need no additional
    ! communications as each rank can just do its own calculation.
    if(to_upper(adjustl(basis)) .eq. 'HF') then
      psi => HFpsi
    elseif(to_upper(adjustl(basis)) .eq. 'CAN') then
      psi => canpsi
    elseif(to_upper(adjustl(basis)) .eq. 'DEN') then
      psi => denpsi
    endif

    do wave=1,nwt_local
      wave_global  = spwf_map(wave)

      ME(wave_global,wave_global) = f(psi(:,:,wave),psi(:,:,wave))
    enddo
  else
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! We have to calculate the full matrix
    designated_rank = -1 
    do B=1,8
      ! First, we designate a rank to do the calculations
      !    = the first rank storing spwfs in a symmetry block
      N = HFBlocks_global(B) ;  if(N.eq.0) cycle
      designated_rank(B) = rank_map(sum(HFBlocks_global(1:B-1))+1)
    enddo

    si = 0
    do B=1,8
      N = HFBlocks_global(B)  ! <= this loops over global spwfs
      if(N.eq.0) cycle
      calc_rank = designated_rank(B) !  This is the rank doing the integrations

      do wave_global=si+1,si+N
        ! Find the rank storing psi_i and its local index
        ranki = rank_map(wave_global)
        wave  = spwf_inverse(wave_global)
        ! transfer psi_i to the calculating rank
#if(USE_MPI>0)
        call Transfer_psi(psi_i, wave, basis,ranki, calc_rank)
#else
        call Transfer_psi(psi_i, wave, basis)
#endif
        do wave2_global=wave_global,si+N
          ! Find the rank storing psi_j and its local index
          rankj = rank_map(wave2_global)
          wave2 = spwf_inverse(wave2_global)
          ! transfer psi_j to the calculating rank
#if(USE_MPI>0)
          call Transfer_psi(psi_j, wave2, basis,rankj, calc_rank)
#else
          call Transfer_psi(psi_j, wave2, basis)
#endif

          if(MPI_RANK.eq.calc_rank) then
            ! Perform the integration.....
            ME(wave_global,wave2_global) = f(psi_i, psi_j)
            ! .... and use the symmetry
            if(wave_global .ne. wave2_global) then
              ME(wave2_global, wave_global) = sym* ME(wave_global, wave2_global)
            endif
          endif
        enddo
      enddo
      si = si + N
    enddo
  endif

#if(USE_MPI>0) 
  ! Transferring results to all ranks
  call MPI_ALLREDUCE(MPI_IN_PLACE, ME, nwt**2, MPI_REAL8, MPI_SUM, &
  &                                                      MPI_COMM_WORLD,mpi_err)
#endif 

end subroutine ME_function

subroutine ME_function_deriv1(ME, f, sym, diag, basis)
  !-----------------------------------------------------------------------------
  ! Evaluate the matrix elements of an operator on the spwfs using a function  
  ! defined for its evaluation between any two. i.e. given 
  !
  !      f(psi_i, psi_j) = < psi_i | O | psi_j >
  !
  ! this function evaluates either the full matrix M_ij = f(psi_i, psi_j) or
  ! just its diagonal elements. 
  !
  ! Assumptions:
  ! 1) A symmetry with sign sym, meaning that 
  !      f(psi_j, psi_i) = sym * f(psi_i , psi_j)
  ! 2) ALL matrix elements are real, hence f should return real values.
  ! 3) One derivative is involved in f, i.e. that procedure takes as arguments
  !       f(psi_i, psi_j, der_j)
  ! 4) the matrix elements are only calculated within symmetry blocks. 
  ! - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -
  ! Input: 
  !     f    : function procedure, with the interface of spin_z_real
  !     sym  : 
  !     diag : if .true., only compute the diagonal matrix elements
  ! Output:
  !     ME   : matrix elements such that
  !               M(i,j) = <i|oper|j> = int d^3 psi_i^*(r) oper(r) psi_j(r)
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Important note
  !    in the case of MPI calculations, this implementation is naive and 
  !    requires quite a bunch of communications between ranks. For now, the 
  !    load balancing is made simply on protons vs. neutrons.
  !-----------------------------------------------------------------------------
  
  real(KIND=dp), intent(out)  :: ME(:,:)
  procedure (angmom_z_real)   :: f
  integer, intent(in)         :: sym
  character(len=*), intent(in):: basis
  logical, intent(in)         :: diag

  integer                :: B, si, N
  integer                :: wave, wave_global, wave2, wave2_global
  integer                :: designated_rank(8), calc_rank, ranki, rankj
  real(KIND=dp), pointer :: psi(:,:,:), dpsi(:,:,:,:)
  real(KIND=dp)          :: psi_i(mv,4), psi_j(mv,4), der_j(mv,3,4)
#if(USE_MPI>0)
  integer :: mpi_err
#endif

  ME = 0.0d0 ! clearly zero everything for the allreduce call later

  if(diag) then
    ! We calculate only diagonal matrix elements; we thus need no additional
    ! communications as each rank can just do its own calculation.
    if(to_upper(adjustl(basis)) .eq. 'HF') then
      psi  => HFpsi
      dpsi => HFdpsi 
    elseif(to_upper(adjustl(basis)) .eq. 'CAN') then
      psi  => canpsi
      dpsi => candpsi 
    elseif(to_upper(adjustl(basis)) .eq. 'DEN') then
      psi  => denpsi
      dpsi => dendpsi 
    endif

    do wave=1,nwt_local
      wave_global  = spwf_map(wave)

      ME(wave_global,wave_global) = &
      &                          f(psi(:,:,wave),psi(:,:,wave),dpsi(:,:,:,wave))
    enddo
  else
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! We have to calculate the full matrix
    designated_rank = -1 
    do B=1,8
      ! First, we designate a rank to do the calculations
      !    = the first rank storing spwfs in a symmetry block
      N = HFBlocks_global(B) ;  if(N.eq.0) cycle
      designated_rank(B) = rank_map(sum(HFBlocks_global(1:B-1))+1)
    enddo

    si = 0
    do B=1,8
      N = HFBlocks_global(B)  ! <= this loops over global spwfs
      if(N.eq.0) cycle
      calc_rank = designated_rank(B) !  This is the rank doing the integrations

      do wave_global=si+1,si+N
        ! Find the rank storing psi_i and its local index
        ranki = rank_map(wave_global)
        wave  = spwf_inverse(wave_global)
        ! transfer psi_i to the calculating rank
#if(USE_MPI>0)
        call Transfer_psi(psi_i, wave, basis,ranki, calc_rank)
#else
        call Transfer_psi(psi_i, wave, basis)
#endif
        do wave2_global=wave_global,si+N
          ! Find the rank storing psi_j and its local index
          rankj = rank_map(wave2_global)
          wave2 = spwf_inverse(wave2_global)

          ! transfer psi_j and its derivative to the calculating rank
#if(USE_MPI>0)
          call Transfer_psi(psi_j,wave2,basis,rankj,calc_rank)
          call Transfer_derpsi_complete(der_j, wave2, basis,rankj,calc_rank)
#else
          call Transfer_psi(psi_j,wave2,basis)
          call Transfer_derpsi_complete(der_j, wave2, basis)
#endif

          if(MPI_RANK.eq.calc_rank) then
            ! Perform the integration.....
            ME(wave_global,wave2_global) = f(psi_i, psi_j, der_j)
            ! .... and use the symmetry
            if(wave_global .ne. wave2_global) then
              ME(wave2_global, wave_global) = sym* ME(wave_global, wave2_global)
            endif
          endif
        enddo
      enddo
      si = si + N
    enddo
  endif

#if(USE_MPI>0) 
  ! Transferring results to all ranks
  call MPI_ALLREDUCE(MPI_IN_PLACE, ME, nwt**2, MPI_REAL8, MPI_SUM, &
  &                                                      MPI_COMM_WORLD,mpi_err)
#endif
end subroutine ME_function_deriv1

subroutine ME_function_deriv2(ME, f, sym, diag, basis)
  !-----------------------------------------------------------------------------
  ! Evaluate the matrix elements of an operator on the spwfs using a function  
  ! defined for its evaluation between any two. i.e. given 
  !
  !      f(psi_i, psi_j) = < psi_i | O | psi_j >
  !
  ! this function evaluates either the full matrix M_ij = f(psi_i, psi_j) or
  ! just its diagonal elements. 
  !
  ! Assumptions:
  ! 1) A symmetry with sign sym, meaning that 
  !      f(psi_j, psi_i) = sym * f(psi_i , psi_j)
  ! 2) ALL matrix elements are real, hence f should return real values.
  ! 3) Two derivatives are involved in f, i.e. that procedure takes as arguments
  !       f(psi_i, der_i, psi_j, der_j)
  ! 4) the matrix elements are only calculated within symmetry blocks. 
  ! - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -
  ! Input: 
  !     f    : function procedure, with the interface of spin_z_real
  !     sym  : 
  !     diag : if .true., only compute the diagonal matrix elements
  ! Output:
  !     ME   : matrix elements such that
  !               M(i,j) = <i|oper|j> = int d^3 psi_i^*(r) oper(r) psi_j(r)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Important note
  !    in the case of MPI calculations, this implementation is naive and 
  !    requires quite a bunch of communications between ranks. For now, the 
  !    load balancing is made simply on protons vs. neutrons.
  !-----------------------------------------------------------------------------
  
  real(KIND=dp), intent(out)  :: ME(:,:)
  procedure (angmom_z_quad)   :: f
  integer, intent(in)         :: sym
  character(len=*), intent(in):: basis
  logical, intent(in)         :: diag

  integer                :: B, si, N
  integer                :: wave, wave_global, wave2, wave2_global
  integer                :: designated_rank(8), calc_rank, ranki, rankj
  real(KIND=dp), pointer :: psi(:,:,:), dpsi(:,:,:,:)
  real(KIND=dp)          :: psi_i(mv,4), psi_j(mv,4)
  real(KIND=dp)          :: der_i(mv,3,4), der_j(mv,3,4)
#if(USE_MPI>0)
  integer :: mpi_err
#endif

  ME = 0.0d0 ! clearly zero everything for the allreduce call later

  if(diag) then
    ! We calculate only diagonal matrix elements; we thus need no additional
    ! communications as each rank can just do its own calculation.
    if(to_upper(adjustl(basis)) .eq. 'HF') then
      psi  => HFpsi
      dpsi => HFdpsi 
    elseif(to_upper(adjustl(basis)) .eq. 'CAN') then
      psi  => canpsi
      dpsi => candpsi 
    elseif(to_upper(adjustl(basis)) .eq. 'DEN') then
      psi => denpsi
      dpsi => dendpsi 
    endif

    do wave=1,nwt_local
      wave_global  = spwf_map(wave)

      ME(wave_global,wave_global) = &
      &                          f(psi(:,:,wave),dpsi(:,:,:,wave), &
      &                            psi(:,:,wave),dpsi(:,:,:,wave))
    enddo
  else
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! We have to calculate the full matrix
    designated_rank = -1 
    do B=1,8
      ! First, we designate a rank to do the calculations
      !    = the first rank storing spwfs in a symmetry block
      N = HFBlocks_global(B) ;  if(N.eq.0) cycle
      designated_rank(B) = rank_map(sum(HFBlocks_global(1:B-1))+1)
    enddo

    si = 0
    do B=1,8
      N = HFBlocks_global(B)  ! <= this loops over global spwfs
      if(N.eq.0) cycle
      calc_rank = designated_rank(B) !  This is the rank doing the integrations

      do wave_global=si+1,si+N
        ! Find the rank storing psi_i and its local index
        ranki = rank_map(wave_global)
        wave  = spwf_inverse(wave_global)
        ! transfer psi_i and its derivative to the calculating rank
#if(USE_MPI>0)
        call Transfer_psi(psi_i,wave, basis,ranki, calc_rank)
        call Transfer_derpsi_complete(der_i, wave, basis,ranki, calc_rank)
#else
        call Transfer_psi(psi_i,wave, basis)
        call Transfer_derpsi_complete(der_i, wave, basis)
#endif
        do wave2_global=wave_global,si+N
          ! Find the rank storing psi_j and its local index
          rankj = rank_map(wave2_global)
          wave2 = spwf_inverse(wave2_global)
          
          ! transfer psi_j and its derivative to the calculating rank
#if(USE_MPI>0)
          call Transfer_psi(psi_j,wave2,basis,rankj,calc_rank)
          call Transfer_derpsi_complete(der_j, wave2, basis,rankj, calc_rank)
#else
          call Transfer_psi(psi_j,wave2,basis)
          call Transfer_derpsi_complete(der_j,wave2, basis)
#endif
          if(MPI_RANK.eq.calc_rank) then
            ! Perform the integration.....
            ME(wave_global,wave2_global) = f(psi_i, der_i, psi_j, der_j)
            ! .... and use the symmetry
            if(wave_global .ne. wave2_global) then
              ME(wave2_global, wave_global) = sym* ME(wave_global, wave2_global)
            endif
          endif
        enddo
      enddo
      si = si + N
    enddo
  endif

#if(USE_MPI>0) 
  ! Transferring results to all ranks
  call MPI_ALLREDUCE(MPI_IN_PLACE, ME, nwt**2, MPI_REAL8, MPI_SUM, &
  &                                                      MPI_COMM_WORLD,mpi_err)
#endif
end subroutine ME_function_deriv2

subroutine Transfer_psi(psi,wave, basis &
#if(USE_MPI>0)
&                       , send_rank, calc_rank)
#else
&                       )
#endif
    !---------------------------------------------------------------------------
    ! Transfer the wavefunction  
    !        denpsi(:,:,wave)   
    ! from the sending MPI_rank to a rank fit for calculations. 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input : 
    !   send_rank : MPI rank storing the requested wavefunction
    !   calc_rank : MPI rank supposed to be doing calculations with the 
    !               requested wavefunction.
    !   wave      : LOCAL index of the requested spwf on the send_rank
    !   basis     : spwf in which basis to send; 'HF', 'CAN' or 'DEN'
    ! Output:
    !  psi        : the requested spwf, but only on CALC_RANK. For all 
    !               other MPI ranks, the result will be unallocated. 
    !---------------------------------------------------------------------------
    real(KIND=dp),  intent(out)  :: psi(mv,4)
    integer, intent(in)          :: wave
    real(KIND=dp), pointer       :: psis(:,:,:)
    character(len=*), intent(in) :: basis
#if(USE_MPI>0)
    integer, intent(in)          :: send_rank, calc_rank
    integer                      :: mpi_err
#endif    

    if(to_upper(adjustl(basis))     .eq. 'HF') then
      psis => HFpsi
    elseif(to_upper(adjustl(basis)) .eq. 'CAN') then
      psis => canpsi
    elseif(to_upper(adjustl(basis)) .eq. 'DEN') then
      psis => denpsi
    endif

#if(USE_MPI>0)
    if((MPI_RANK.eq. calc_rank) .AND. (send_rank.eq.calc_rank)) then
        ! nothing to send or receive
        psi   = psis(:,:,wave)
    elseif(MPI_RANK.eq.calc_rank) then
        ! calc_rank receives
        call MPI_RECV(            psi, 4*mv, MPI_REAL8, send_rank, 2, &
        &                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    else if(MPI_RANK .eq. send_rank)  then
        ! send_rank sends the wavefunction
        call MPI_SEND(psis(:,:,wave), 4*mv, MPI_REAL8, calc_rank, 2,&
        &                                           MPI_COMM_WORLD, mpi_err)
    endif
#else 
    psi   = psis(:,:,wave)
#endif

end subroutine Transfer_psi

subroutine Transfer_derpsi(derpsi,wave,direction, basis, TR &
#if(USE_MPI>0)
&                       , send_rank, calc_rank)
#else
&                       )
#endif
    !---------------------------------------------------------------------------
    ! Transfer the (possibly time-reversed) derivative of a wavefunction  
    !        dendpsi(:,:,direction,wave)   
    ! from the sending MPI_rank to a rank fit for calculations. 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input : 
    !   send_rank : MPI rank storing the requested wavefunction
    !   calc_rank : MPI rank supposed to be doing calculations with the 
    !               requested wavefunction.
    !   wave      : LOCAL index of the requested spwf on the send_rank
    !   direction : direction of the derivative (x/y/z = 1/2/3) of the 
    !               requested spwf
    !   basis     : character, 
    !   TR        : logical, whether or not to apply a time-reversal operation
    !               before returning.
    ! Output:
    !  derpsi     : the requested derivative of an spwf, but only on CALC_RANK.
    !               For all other MPI ranks, the array is not changed.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(out)   :: derpsi(mv,4)
    integer, intent(in)          :: wave
    integer, intent(in)          :: direction
    logical, intent(in)          :: TR
    real(KIND=dp), pointer       :: psis(:,:,:,:)
    character(len=*), intent(in) :: basis
#if(USE_MPI>0)
    integer, intent(in)          :: send_rank, calc_rank
    integer                      :: mpi_err
#endif    

    if(to_upper(adjustl(basis))     .eq. 'HF') then
      psis => HFdpsi
    elseif(to_upper(adjustl(basis)) .eq. 'CAN') then
      psis => candpsi
    elseif(to_upper(adjustl(basis)) .eq. 'DEN') then
      psis => dendpsi
    endif

#if(USE_MPI>0)
    if((MPI_RANK.eq. calc_rank) .AND. (send_rank.eq.calc_rank)) then
        ! nothing to send or receive
        derpsi   = psis(:,direction,:,wave)
        if(TR)   derpsi = TimeReverse(derpsi)
    elseif(MPI_RANK.eq.calc_rank) then
        ! calc_rank receives
        call MPI_RECV(                     derpsi, 4*mv, MPI_REAL8,send_rank,2,&
        &                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
        ! and time-reverses if needed
        if(TR)   derpsi = TimeReverse(derpsi)
    else if(MPI_RANK .eq. send_rank)  then
        ! ranki sends the wavefunction
        call MPI_SEND(psis(:,direction,:,wave), 4*mv, MPI_REAL8,calc_rank,2,&
        &                                           MPI_COMM_WORLD, mpi_err)
    endif
#else 
    derpsi   = psis(:,direction,:,wave)
    if(TR)   derpsi = TimeReverse(derpsi)
#endif

end subroutine Transfer_derpsi

subroutine Transfer_derpsi_complete(derpsi, wave, basis &
#if(USE_MPI>0)
&                       , send_rank, calc_rank)
#else
&                       )
#endif
    !---------------------------------------------------------------------------
    ! Transfer the complere gradient of a wavefunction  
    !        HFdpsi/candpsi/dendpsi(:,:,:,wave)   
    ! from the sending MPI_rank to a rank fit for calculations. 
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input : 
    !   send_rank : MPI rank storing the requested wavefunction
    !   calc_rank : MPI rank supposed to be doing calculations with the 
    !               requested wavefunction.
    !   wave      : LOCAL index of the requested spwf on the send_rank
    !   basis     : spwf in which basis to send; 'HF', 'CAN' or 'DEN'
    ! Output:
    !  derpsi     : the requested derivative of an spwf, but only on CALC_RANK.
    !               For all other MPI ranks, the array is not changed.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(out)   :: derpsi(mv,3,4)
    integer, intent(in)          ::  wave
    character(len=*), intent(in) :: basis
    real(KIND=dp), pointer       :: psis(:,:,:,:)
#if(USE_MPI>0)
    integer, intent(in)          :: send_rank, calc_rank
    integer                      :: mpi_err
#endif

    if(to_upper(adjustl(basis))     .eq. 'HF') then
      psis => HFdpsi
    elseif(to_upper(adjustl(basis)) .eq. 'CAN') then
      psis => candpsi
    elseif(to_upper(adjustl(basis)) .eq. 'DEN') then
      psis => dendpsi
    endif

#if(USE_MPI>0)
    if((MPI_RANK.eq. calc_rank) .AND. (send_rank.eq.calc_rank)) then
        ! nothing to send or receive
        derpsi   = psis(:,:,:,wave)
    elseif(MPI_RANK.eq.calc_rank) then
        ! calc_rank receives
        call MPI_RECV(derpsi          , 12*mv, MPI_REAL8,send_rank,2,&
        &                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    else if(MPI_RANK .eq. send_rank)  then
        ! ranki sends the wavefunction
        call MPI_SEND(psis(:,:,:,wave), 12*mv, MPI_REAL8,calc_rank,2,&
        &                                           MPI_COMM_WORLD, mpi_err)
    endif
#else 
    derpsi   = psis(:,:,:,wave)
#endif

end subroutine Transfer_derpsi_complete

!-------------------------------------------------------------------------------
! A natural place for this routine would be in basis_transform.f90, but it is 
! needed in this module....
!-------------------------------------------------------------------------------
 
function transform_mat_diag(M, transfo) result(Mc)
  !-----------------------------------------------------------------------------
  ! Identical to transform_mat, but only calculate the diagonal matrix elements.
  ! 
  ! Input:
  !     M    : matrix to transform
  !  transfo : unitary transformation C to employ 
  !            (in the conventions of this module)
  !
  ! Output:
  !     the diagonal elements of Mc = C^T M C
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in)    :: M(nwt,nwt)
  real(KIND=dp), intent(in)    :: transfo(nwt,nwt)
  real(KIND=dp)                :: Mc(nwt)
  integer                      :: B, N, si,i,j,l

  si = 0
  Mc = 0.0d0
  do B=1,8
    N = HFBlocks_global(B)  ;  if(N .eq. 0) cycle 

    do i=si+1,si+N
      do j=si+1,si+N
        do l=si+1,si+N
          !                     c^T           M         C
           Mc(i)  = Mc(i) + transfo(l,i) * M (l,j) * HFtransfo(j,i)
        enddo
      enddo
    enddo

    si = si +  N
  enddo  

 end function transform_mat_diag 
 
 function memory_wavefunctions(spwf_number) result(storage)
  !-----------------------------------------------------------------------------
  ! Estimate the total storage requirements for a given number of spwfs.
  !
  ! Input:
  !   spwf_number : number of spwfs stored
  ! Output:
  !   storage: total number of real numbers involved in storing the spwfs
  !-----------------------------------------------------------------------------
  integer, intent(in)    :: spwf_number
  integer(kind=LargeInt) :: storage

  ! storage for the HFPSI array 
  storage = mv * 4 * spwf_number
  ! factors two account for
  !   (*) additional storage of momentum updates for heavy-ball machinery
  !   (*) additional storage for h | psi > in evolution
  storage = 4 * storage
  
  ! 2 bonus wavefunctions in estimation of iterative parameters
  storage = 4*mv*2 + 3*mv*4*2 + 6*mv*2 

  ! storage for the first order derivatives
  storage = storage + 3 * mv * 4 * spwf_number
  ! storage for the second order derivatives
  storage = storage + 6 * mv * 4 * spwf_number

 end function memory_wavefunctions
 
 function memory_wavefunctions_2D(spwf_number, row, nrow, col, ncol) result(mem)
   ! TODO: document  
   integer,intent(in)     :: row, nrow, col, ncol, spwf_number
   integer(kind=LargeInt) :: mem
#if(USE_MPI>0)
   integer(kind=LargeInt) :: xs,ys
   integer, external :: numroc
 
   xs = NUMROC(          4*mv,block_factor_row,row,0,nrow)
   ys = NUMROC(spwf_number,block_factor_col,col,0,ncol)
  
   mem = 3*xs*ys ! factor 3 for mom_2D and hpsi_2D
#else
   mem = 0
#endif
   end function memory_wavefunctions_2D

  subroutine clean_wavefunctions()

    if(allocated(HFPsi))    deallocate(HFPsi)
    if(allocated(HFdPsi))   deallocate(HFdPsi)
    if(allocated(HFddPsi))  deallocate(HFddPsi)
    if(allocated(HFdddPsi)) deallocate(HFdddPsi)

    if(allocated(CANPsi))    deallocate(CANPsi)
    if(allocated(CANdPsi))   deallocate(CANdPsi)
    if(allocated(CANddPsi))  deallocate(CANddPsi)
    if(allocated(CANdddPsi)) deallocate(CANdddPsi)

    if(allocated(spenergies))  deallocate(spenergies)
    if(allocated(dispersions)) deallocate(dispersions)
    if(allocated(canenergies)) deallocate(canenergies)

    if(allocated(sx)) deallocate(sx)
    if(allocated(sy)) deallocate(sy)
    if(allocated(sz)) deallocate(sz)

  end subroutine clean_wavefunctions

end module wavefunctions
