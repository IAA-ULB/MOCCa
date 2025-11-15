!!! ============================================================================
!!! modified nil8.f90 for building a python module ni8_f90 that exposes
!!! `subroutine nilsson`.
!!! triple ! indicate comments relating to the modifications

!!! As we rely on numpy.f2py for building Python modules from Fortran code,
!!! and we want numpy arrays as extent(inout) arguments the memory management
!!! of these arguments must be taken out of the subroutine. Here is the full 
!!! list of arguments and their declarations. The `allocatable` arguments pose  
!!! problems because they are allocated on the fly inside `subroutine nilsson`
!!! or in between the two calls
!!!
!!!  -> real(kind=8),  allocatable, intent(inout):: wfs(:,:,:)
!!!         `wfs` corresponds to `HFPsi.data`. it is allocated between the two
!!!         calls.
!!!  -> integer,       allocatable, intent(inout):: kparz(:)
!!!         allocated inside the first call as `kparz(nwt)`, as `nwt` equals
!!!         `nwp` + `nwn` this can be allocated before the first call.
!!!  -> real(kind=8),  allocatable, intent(inout):: esp1(:)
!!!         same as `kparz`.
!!!     integer                   , intent(in)   :: meven
!!!     integer                   , intent(in)   :: modd
!!!     integer                   , intent(in)   :: nwt
!!!     integer                   , intent(in)   :: nwp
!!!     integer                   , intent(in)   :: nwn
!!!     integer                   , intent(in)   :: nwp
!!!     integer                   , intent(in)   :: nwn
!!!     integer                   , intent(in)   :: mx
!!!     integer                   , intent(in)   :: my
!!!     integer                   , intent(in)   :: mz
!!!     real(kind=8)             , intent(in)   :: dx
!!!     real(kind=8)             , intent(in)   :: osc_freq(3)
!!!  -> integer, allocatable      , intent(in)   :: spwf_map(:)
!!!         allocated between the two calls.
!!!
!!! ============================================================================

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
     !         ATTENTION: this set of spwfs will not be orthonormal at all.
     !    spe: a guess (trivial in this routine) of the single-particle energies
     !         of these random states. 
     !    par: the parity quantum numbers of the spwfs
     !--------------------------------------------------------------------------
      
     !!! Argument declaration: order according to the arguents in the call
     !!! replace kind=8 with kind=8
     !!! remove allocatable qualifier
     real(kind=8), intent(inout) :: psi(:,:,:)
     integer     , intent(inout) :: par(:)
     real(kind=8), intent(inout) :: spe(:)
     integer     , intent(in)    :: nshells_even, nshells_odd
     integer     , intent(in)    :: nw,nwn,nwp,neut,prot,mx,my,mz
     real(kind=8), intent(in)    :: dx, osc_freq(3)
     integer     , intent(in)    :: map(:)

     !!! Internal variable declaration

     real(kind=8), allocatable   :: trash(:,:), scalar_trash
     integer                     :: s, minindex, k
     integer, allocatable        :: seed(:)
     !!! copied from
     character(len=20)           :: random_numbers = 'FAST'

     integer :: MPI_RANK
     MPI_rank = 1

     allocate(trash(mx*my*mz,4))

     !---------------------------------------------------------------------------
     ! Use these dummy arguments to stop compiler complaints
!     scalar_trash = mx*my*mz*dx*product(osc_freq)*nshells_even*nshells_odd
!     scalar_trash = neut * prot

     !!!<Memory management of arguments is done in the Python caller
     !!! if(allocated(par)) deallocate(par)
     !!! if(allocated(spe)) deallocate(spe)
     !!! allocate(par(nw),spe(nw))
     !!!>

     ! We just take spwfs of positive and negative in 50/50 proportion
     par(            1:nwn/2    ) = +1
     par(nwn/2      +1:nwn      ) = -1
     par(nwn        +1:nwn+nwp/2) = +1
     par(nwn  +nwp/2+1:nw       ) = -1
     
     spe = 100
     
     !!!<Instead of testing allocation we test for non-zero size
     !!! if(allocated(map)) then
     if (size(map).gt.0) then
     !!!>

        call random_seed(size=s)
        allocate(seed(s))

        if(adjustl(random_numbers) .eq. 'FAST') then
        ! Let all MPI ranks generate different sets of random numbers
            seed = 987654321 + MPI_RANK * 123456789 ! Seed value needs to depend on
                                                    ! MPI_RANK; if not, we all ranks
                                                    ! will generate the same sequence
                                                    ! and we will run in trouble with
                                                    ! orthonormalisation
            call random_seed(put=seed)
        else
            ! Let all MPI ranks generate the same random numbers, but take a subset
            seed = 987654321                        ! Seed value is fixed
            call random_seed(put=seed)

            ! Determine the offset for this particular MPI rank
            minindex = minval(map)
            do k=1,minindex-1
              call random_number(trash)
            enddo
        endif
        ! Generate random numbers as (non-orthogonal) wavefunctions
        call random_number(psi)
     endif
  end subroutine randomspwfs
