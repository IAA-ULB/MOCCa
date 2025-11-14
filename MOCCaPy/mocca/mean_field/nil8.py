### This is a python translation of nil8.f90
### We are staying as close as possible to the original.
### Triple # indicate new comment lines.
### Single # indicate original lines commented out.
### the many python loops make it terribly slow.

### import numpy as np ### np is used as a variable by nilsson()
import numpy

# module nil8
#  !==============================================================================
#  !  #######   ##   #    # #####   ##   #      #    #  ####
#  !     #     #  #  ##   #   #    #  #  #      #    # #
#  !     #    #    # # #  #   #   #    # #      #    #  ####
#  !     #    ###### #  # #   #   ###### #      #    #      #
#  !     #    #    # #   ##   #   #    # #      #    # #    #
#  !     #    #    # #    #   #   #    # ######  ####   ####
#  !
#  !  Copyright W. Ryssens & M. Bender
#  !
#  !==============================================================================
#  !
#  ! Module able to generate new wave-functions from a Nilsson model
#  ! Hamiltonian.  Code was taken from nil8 (v1.0.0), but I do not
#  ! guarantee it actually diagonalizes the correct nilsson Hamiltonian.
#  !
#  ! In any case, more info can (in principle) be found in
#  !
#  !==============================================================================
 
# use compilation
### `compilation` defines sp (single precision = numpy.float32) and dp (double
### precision = float or numpy.float64)
 
# implicit none
  
# contains

# subroutine nilsson (wfs,kparz,esp1,meven,modd,nwt,nwp,nwn,npp,npn,mx,my,mz,   &
#  &                   dx,osc_freq, spwf_map)
def nilsson(wfs, kparz, esp1, meven, modd,
            nwt, nwp, nwn, npp, npn, mx, my, mz, dx, osc_freq, spwf_map):
    """
    !---------------------------------------------------------------------------
    ! Subroutine taken from nil8.1.0.0.f, written by
    !         Bonche, Flocard and Heenen
    ! in the depths of time.
    !
    ! The wave-functions are generated on a (mx,my,mz) mesh with spacing dx.
    ! They have
    !   * good isospin
    !   * good Rz signature quantum number (all +i)
    !   * good parity quantum number (value in kparz on exit)
    !
    ! Note that they are thus generated on an EV8-like mesh and further actions
    ! should be taken if the user wants other quantum numbers.
    !
    ! How this routine works:
    !   Step 1) Enumerate all of the eigenstates of the harmonic oscillator
    !           that exist within the first [meven] and [modd] shells.
    !           They are enumerated as (nx,ny,nz), where these three numbers
    !           are the oscillator quanta.
    !   Step 2) The corresponding 1D Hermite functions are constructed on the
    !           mesh in a rather ad-hoc way. Rather than use the recursive
    !           relations, Hn(x) is simply multiplied by x and afterwards
    !           orthonormalized.
    ! |-Step 3) The matrix elements of the Nilsson Hamiltonian in the basis of
    ! |          these states are calculated.
    ! | Step 4) The hamiltonian is diagonalized.
    ! | Step 5) The lowest nwn/p states are actually constructed from the
    ! |          hermite functions constructed before.
    ! |- Isospin loop.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   wfs:
    !       allocatable array containing the constructed wave-functions on exit
    !       if spwf_map is allocated on entry. not referenced if spwf_map is
    !       not allocated.
    !   kparz:
    !       allocatable array containing the parities on exit
    !       will be re-allocated to match the rest of the input
    !   esp1  :
    !       single-particle energies on exit
    !       will be re-allocated to match the rest of the input
    !   meven, modd:
    !       the number of oscillator shells with even/odd parity
    !   nwt, nwn, nwp:
    !       the number of total/neutron/proton wave-functions on exit
    !   npp, npn
    !       the number of protons/neutrons, needed for the hamiltonian
    !       parameters
    !   mx,my,mz
    !       number of mesh points in the x/y/z direction (on the EV8 mesh )
    !   dx
    !       mesh spacing on the EV8 mesh
    !   homegax,homegay,homegaz
    !       harmonic oscillator parameters
    !   spwf_map
    !       allocatable, integer
    !       This controls the construction of spwfs. If allocated, it will
    !       construct the 'spwf_map(i)'-th spwf in the Nilsson spectrum
    !       and store it in the i-th entry in the wfs array. Obviously
    !       len(spwf_map) <= nwt.
    !---------------------------------------------------------------------------
    !implicit real*8 (a-h,o-z)
    !
    !  101 format (/,' neutron levels kappa=',e10.3,' mu=',e9.2,   &
    !     &          ' al0n=al0*(1+',f5.2,'*(n-z)/a)',/,           &
    !     &          ' (n0,nor,energy/(hbar*omega0),parity)',/,' ')
    !  102 format (/,' proton  levels kappa=',e10.3,' mu=',e9.2,   &
    !     &          ' al0p=al0*(',f6.3,'-',f5.2,'*(n-z)/a)',/,    &
    !     &          ' (n0,nor,energy/(hbar*omega0),parity)',/,' ')
    !  103 format (' (',2i4,f8.3,i3,') (',2i4,f8.3,i3,') (',2i4,f8.3,i3,')')
    """

    # external :: DSYEV
    ### LAPACK:DSYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices.

    # integer              , intent(in)        :: meven, modd,mx,my,mz,nwt,nwp,nwn
    # integer              , intent(in)        :: npp, npn
    # integer, allocatable, intent(inout)      :: kparz(:)
    # real(KIND=dp), allocatable, intent(inout):: wfs(:,:,:), esp1(:)
    # real(KIND=dp), intent(in)                :: osc_freq(3), dx
    # integer, allocatable, intent(in)         :: spwf_map(:)
    #
    # real(KIND=dp)              :: hox, hoy, hoz
    # real(KIND=dp), allocatable :: h(:,:), s(:,:), d(:), wd(:), e(:)
    # real(KIND=dp), allocatable :: he(:,:,:) , a(:), work(:)
    # real(KIND=dp)              :: psi(mx,my,mz,4), hbm(2), cf(2)
    # real(KIND=dp)              :: xho, x, x1, x2, x3, x4, hex, hey, hez, xph
    # real(KIND=dp)              :: ho0, ax, ay, az, y,an, am, xis
    #
    # integer                    :: npar(2,2), nvv, nz2, nz1, nx1, nx2, ny1, ny2
    # integer                    :: nwave, nodd, nnn2, nnn1, nn2, nn1, nn, nmax
    # integer                    :: nij,i,i1,ia,ii,it, iwave,ix, nb, n, kk, iy, iz
    # integer                    :: j,ja, k, nw, neven, ni, ni1, np, nvec, ind
    # integer                    :: mblc, mq, mqa, ms, nblc, ndd, ndim, ifail
    # integer                    :: lwork, store_counter
    # integer, allocatable       :: nsi(:,:),ns(:), nx(:), ny(:), nz(:), irep(:)
    # integer, allocatable       :: nor(:), npa(:), ntrs(:)
    #
    # real*8, parameter :: hhbar=6.58218d0, xxmn =1.044673d0
    hhbar, xxmn = 6.58218, 1.044673
    # real*8, parameter :: ca   =0.986d0  , cb   =0.14
    ca, cb = 0.986, 0.14
    # real*8, parameter :: xk(4)  = (/0.08d0,0.08d0,   0.0637d0,0.0637d0/)
    xk = numpy.array([0.08, 0.08, 0.0637, 0.0637], dtype=float)
    # real*8, parameter :: xmu(4) = (/0.08d0,0.08d0,   0.0637d0,0.0637d0/)
    xmu = numpy.array([0.08, 0.08, 0.0637, 0.0637], dtype=float)

  # hox = osc_freq(1) ; hoy = osc_freq(2) ; hoz = osc_freq(3)
    hox = osc_freq[0] ; hoy = osc_freq[1] ; hoz = osc_freq[2]

    mblc =  meven+modd+1
  # ms   = (mblc*(mblc**2-1))/6
    ms   = (mblc*(mblc**2-1))//6
  # ndim = (mblc*(mblc-1))/2
    ndim = (mblc*(mblc-1))//2
  # mqa  = (ms*(3*mblc**2-2))/10
    mqa  = (ms*(3*mblc**2-2))//10
    mq   = my * mx * mz * 4

    # if(ms .lt. nwn .or. ms .lt. nwp) then
    #   ! The basis does not allow for this many states
    #   print *, ' Harmonic oscillator basis in nilsson is not sufficient.'
    #   stop
    # endif
    if ms < nwn or \
       ms < nwp:
        raise ValueError("Harmonic oscillator basis in nilsson is not sufficient.")

    # hbm(1)  = hhbar*hhbar/xxmn
    # hbm(2)  = hhbar*hhbar/xxmn
    hbm =numpy.array( [hhbar*hhbar/xxmn,hhbar*hhbar/xxmn], dtype=float)

    # allocate(h(ndim,ndim),s(ndim,ndim),d(ndim),wd(ndim))
    ndim_x_ndim = (ndim,ndim)
    h  = numpy.zeros(ndim_x_ndim, dtype=float)
    s  = numpy.zeros(ndim_x_ndim, dtype=float)
    d  = numpy.zeros(ndim       , dtype=float)
    wd = numpy.empty(ndim       , dtype=float)

    # allocate(nsi(mblc+1,4),ns(mblc+1))
    nsi = numpy.zeros((mblc+1, 4), dtype=int)
    ns  = numpy.zeros( mblc+1    , dtype=int)

    # allocate(nx(ms),ny(ms), nz(ms), e(ms), nor(ms),npa(ms))
    nx  = numpy.zeros(ms, dtype=int)
    ny  = numpy.zeros(ms, dtype=int)
    nz  = numpy.zeros(ms, dtype=int)
    e   = numpy.zeros(ms, dtype=float)
    nor = numpy.zeros(ms, dtype=int)
    npa = numpy.zeros(ms, dtype=int)

    # allocate(he(mblc,max(mx,my,mz),3), a(mqa))
    he = numpy.zeros((mblc, max(mx, my, mz), 3), dtype=float, order='F')
    a  = numpy.zeros(mqa, dtype=float)

    # allocate(irep(mblc+1), ntrs(ms))
    irep = numpy.zeros(mblc + 1, dtype=int)
    ntrs = numpy.zeros(ms, dtype=int)

  # if(allocated(kparz)) deallocate(kparz)
  # if(allocated(esp1))  deallocate(esp1)
  # allocate(kparz(nwt),esp1(nwt))
    ### those were allocated already in HFPsi.__init__

  # irep = 0 ; ntrs = 0                                   ### allocated with numpy.zeros.
  # h = 0.0d0 ; s = 0.0d0 ; d = 0.0d0                     ### allocated with numpy.zeros.
  # nsi = 0 ; ns = 0                                      ### allocated with numpy.zeros.
  # nx = 0 ; ny = 0 ; nz = 0 ; e = 0.0d0; nor =0 ; npa =0 ### allocated with numpy.zeros.
  # he = 0.0d0 ; kparz=0; a= 0.0d0                        ### allocated with numpy.zeros.
    kparz = 0

    # ! In order for the compiler not to complain about non-initialised stuff.
    nvv = 0
    # !c......................... mz must be larger or equal than both mx and my
    #
    # !c     neven   nodd   nvec+   nvec-   nblc    ms   mblc   ndim     mqa
    # !c       1       0      1       0       1      1     2       1       1
    # !c       1       1      0       3       2      4     3       3      10
    # !c       2       1      6       0       3     10     4       6      46
    # !c       2       2      0      10       4     20     5      10     146
    # !c       3       2     15       0       5     35     6      15     371
    # !c       3       3      0      21       6     56     7      21     812
    # !c       4       3     28       0       7     84     8      28    1592
    # !c       4       4      0      36       8    120     9      36    2892
    # !c       5       4     45       0       9    165    10      45    4917
    # !c       5       5      0      55      10    220    11      55    7942
    # !c       6       5     66       0      11    286    12      66   12298
    # !c       6       6      0      78      12    364    13      78   18382
    # !c       7       6     91       0      13    455    14      91   26663
    # !c       7       7      0     105      14    560    15     105   37688

    # !c............. this subroutine determines the starting point by selecting
    # !c              the nwn(nwp) lowest eigenstates of the nilsson hamiltonian
    # !c              h = sum(1 to 3) of hoi*(ni+1/2)-xk*ho0*(2*l*s+xmu*l**2)
    # !c              the quantities hox(hoy,hoz) are in fact m*omegax(y,z)/hbar
    # !c              are given in data
    # !c              (G.Gustafson, I.L.Lamm, B.Nilsson and S.G.Nilsson,
    # !c               Ark Fys 36(1966)613)
    # !c              the operator l is the stretched angular momentum.
    # !c              the operator l**2 is corrected so that its trace over each
    # !c               major shell is 0.0d0.(see also copybook n0 17)


    # !c................................... neven = number of even parity shells
    # !c                                    nodd  = number of odd  parity shells

    neven = meven
    nodd  = modd
    nmax  = max(neven,nodd)
      
    # !c..................................................... ordering the basis
    # nvec  = 0
    # nblc  = 0
    nvec, nblc = -1, -1 ### because of zero based indexing in Python
    # ns(1) = 0
    ns[0] = 0 ### mind the change from 1-based (fortran) indexing to 0-based (python)

  # !c..................................................... loop on the blocks
  # !c                                                             even parity
  # do ni=1,neven
    for ni in range(1,neven+1):
        n    = 2*(ni-1)
        ### nblc was initialized to -1 above
        nblc = nblc + 1
      # nsi(nblc,1) = ((n+2)*(n+4))/8
        nsi[nblc,0] = ((n+2)*(n+4))//8
      # nsi(nblc,2) = ((n+2)*n)/8
        nsi[nblc,1] = ((n+2)*n)//8
      # nsi(nblc,3) = nsi(nblc,2)
        nsi[nblc,2] = nsi[nblc,1]
      # nsi(nblc,4) = nsi(nblc,2)
        nsi[nblc,3] = nsi[nblc,1]

      # !c............................................................ sub-block 1
      # do i=1,ni
        for i in range(1,ni+1):
            # nij = ni - i + 1
            nij = ni - i + 1
            # do j=1,nij
            for j in range(1,nij+1):
                ### nvec was set to -1 above
                nvec = nvec + 1
              # nx(nvec) = 2*(i-1)
                nx[nvec] = 2*(i-1)
              # ny(nvec) = 2*(j-1)
                ny[nvec] = 2*(j-1)
              # nz(nvec) = n - nx(nvec) - ny(nvec)
                nz[nvec] = n - nx[nvec] - ny[nvec]
          # enddo
      # enddo

    # !c............................................................ sub-block 2
        ni1 = ni - 1
      # if (ni1.ne.0) then
        if ni1 != 0:
          # do i=1,ni1
            for i in range(1,ni1+1):
                nij = ni - i
              # do j=1,nij
                for j in range(1,nij+1):
                    nvec = nvec + 1
                  # nx(nvec) = 2*(i-1) + 1
                    nx[nvec] = 2*(i-1) + 1
                  # ny(nvec) = 2*(j-1) + 1
                    ny[nvec] = 2*(j-1) + 1
                  # nz(nvec) = n - nx(nvec) - ny(nvec)
                    nz[nvec] = n - nx[nvec] - ny[nvec]
              # enddo
          # enddo

    # !c............................................................ sub-block 3
          # do i=1,ni1
            for i in range(1,ni1+1):
                nij = ni - i
              # do j=1,nij
                for j in range(1,nij+1):
                    nvec = nvec + 1
                  # nx(nvec) = 2*(i-1) + 1
                    nx[nvec] = 2*(i-1) + 1
                  # ny(nvec) = 2*(j-1)
                    ny[nvec] = 2*(j-1)
                  # nz(nvec) = n - nx(nvec) - ny(nvec)
                    nz[nvec] = n - nx[nvec] - ny[nvec]
              # enddo
          # enddo

    # !c............................................................ sub-block 4
          # do i=1,ni1
            for i in range(1,ni1+1):
                nij = ni - i
    #           do j=1,nij
                for j in range(1,nij+1):
                    nvec = nvec + 1
                  # nx(nvec) = 2*(i-1)
                    nx[nvec] = 2*(i-1)
                  # ny(nvec) = 2*(j-1) + 1
                    ny[nvec] = 2*(j-1) + 1
                  # nz(nvec) = n - nx(nvec) - ny(nvec)
                    nz[nvec] = n - nx[nvec] - ny[nvec]
              # enddo
          # enddo
      # endif
      # ns(nblc+1) = ns(nblc) + ((n+1)*(n+2))/2
        ns[nblc+1] = ns[nblc] + ((n+1)*(n+2))/2
  # enddo
  # !c............................................................. odd parity
  # do ni=1,nodd
    for ni in range(1,nodd+1):
        n = 2*ni - 1
        nblc = nblc + 1
      # nsi(nblc,1) = ((n+1)*(n+3))/8
        nsi[nblc,0] = ((n+1)*(n+3))/8
      # nsi(nblc,2) = ((n-1)*(n+1))/8
        nsi[nblc,1] = ((n-1)*(n+1))/8
      # nsi(nblc,3) = nsi(nblc,1)
        nsi[nblc,2] = nsi[nblc,0]
      # nsi(nblc,4) = nsi(nblc,1)
        nsi[nblc,3] = nsi[nblc,0]

    # !c............................................................ sub-block 1
      # do i=1,ni
        for i in range(1,ni+1):
            nij = ni - i + 1
          # do j=1,nij
            for j in range(1,nij+1):
                nvec = nvec + 1
              # nx(nvec) = 2*(i-1)
                nx[nvec] = 2*(i-1)
              # ny(nvec) = 2*(j-1)
                ny[nvec] = 2*(j-1)
              # nz(nvec) = n - nx(nvec) - ny(nvec)
                nz[nvec] = n - nx[nvec] - ny[nvec]
          # enddo
      # enddo

    # !c............................................................ sub-block 2
        ni1 = ni - 1
      # if (ni1.ne.0) then
        if ni1 != 0:
          # do i=1,ni1
            for i in range(1,ni1+1):
                nij = ni - i
              # do j=1,nij
                for j in range(1,nij+1):
                    nvec = nvec + 1
                  # nx(nvec) = 2*i-1
                    nx[nvec] = 2*i-1
                  # ny(nvec) = 2*j-1
                    ny[nvec] = 2*j-1
                  # nz(nvec) = n - nx(nvec) - ny(nvec)
                    nz[nvec] = n - nx[nvec] - ny[nvec]
              # enddo
          # enddo
      # endif

    # !c............................................................ sub-block 3
      # do i=1,ni
        for i in range(1,ni+1):
            nij = ni - i + 1
          # do j=1,nij
            for j in range(1,nij+1):
                nvec = nvec + 1
              # nx(nvec) = 2*i-1
                nx[nvec] = 2*i-1
              # ny(nvec) = 2*(j-1)
                ny[nvec] = 2*(j-1)
              # nz(nvec) = n - nx(nvec) - ny(nvec)
                nz[nvec] = n - nx[nvec] - ny[nvec]
          # enddo
      # enddo

    # !c............................................................ sub-block 4
      # do i=1,ni
        for i in range(1,ni+1):
            nij=ni - i + 1
          # do j=1,nij
            for j in range(1,nij+1):
                nvec = nvec + 1
              # nx(nvec) = 2*(i-1)
                nx[nvec] = 2*(i-1)
              # ny(nvec) = 2*j-1
                ny[nvec] = 2*j-1
              # nz(nvec) = n - nx(nvec) - ny(nvec)
                nz[nvec] = n - nx[nvec] - ny[nvec]
          # enddo
      # enddo
      # ns(nblc+1) = ns(nblc) + ((n+1)*(n+2))/2
        ns[nblc+1] = ns[nblc] + ((n+1)*(n+2))/2
  # enddo

  # !c............................. one dimensional oscillator wave-functions
  # xis   =     (npn-npp)
  # xis   = xis/(npn+npp)
    xis   = (npn-npp)/(npn+npp)
  # cf(2) = ca   -cb*xis
  # cf(1) = 1.0d0+cb*xis
    cf = numpy.array([1.0 + cb*xis, ca - cb*xis])

  # !c.................................................... loop on the isospin
    nwave         = 0
    store_counter = 1
    # do it=1,2
    for it in range(2):
        nn = max(mx, my, mz)
      # do ndd=1,3
        for ndd in range(3):
            # if (ndd.eq.1) xho = sqrt(hox)*dx
            # if (ndd.eq.2) xho = sqrt(hoy)*dx
            # if (ndd.eq.3) xho = sqrt(hoz)*dx
            # xho = xho * sqrt(cf(it))
            xho = numpy.sqrt(hox * cf[it]) if ndd == 0 else \
                  numpy.sqrt(hoy * cf[it]) if ndd == 1 else \
                  numpy.sqrt(hoz * cf[it]) * dx
          # do j=1,nmax
            for j in range(nmax):
              # nn1 = 2*j - 1
                nn1 = 2*j
              # nn2 = 2*j
                nn2 = 2*j + 1
              # if (j.eq.1) then
                if j == 0:
                  # do i=1,nn
                    for i in range(nn):
                      # x = xho*(i-0.5d0)
                        x = xho*(i+0.5) ### i starts at 0, henc +0.5
                      # he(1,i,ndd) = exp(-x*x/2.0d0)
                        he[0,i,ndd] = numpy.exp(-0.5*x*x)
                      # he(2,i,ndd) = x*he(1,i,ndd)
                        he[1,i,ndd] = x*he[0,i,ndd]
                  # enddo
                else:
                  # do i=1,nn
                    for i in range(nn):
                      # x = (xho*(i-0.5d0))**2
                        x = (xho*(i+0.5))**2 ### i starts at 0, hence + 0.5
                      # he(nn1,i,ndd) = he(nn1-2,i,ndd)*x
                        he[nn1,i,ndd] = he[nn1-2,i,ndd]*x
                      # he(nn2,i,ndd) = he(nn2-2,i,ndd)*x
                        he[nn2,i,ndd] = he[nn2-2,i,ndd]*x
                  # enddo
                  # do k=1,j-1
                    for k in range(j-1):
                      # nnn1 = 2*k-1
                        nnn1 = 2*k
                      # nnn2 = 2*k
                        nnn2 = 2*k+1
                      # x1   = 0.0d0
                      # x2   = 0.0d0
                        x1, x2 = 0.0, 0.0
                      # do i=1,nn
                        for i in range(nn):
                          # x1 = x1 + he(nn1,i,ndd)*he(nnn1,i,ndd)
                            x1 = x1 + he[nn1,i,ndd]*he[nnn1,i,ndd]
                          # x2 = x2 + he(nn2,i,ndd)*he(nnn2,i,ndd)
                            x2 = x2 + he[nn2,i,ndd]*he[nnn2,i,ndd]
                      # enddo
                      # x1 = x1*dx*2.0d0
                        x1 *= dx*2.0
                      # x2 = x2*dx*2.0d0
                        x2 *= dx*2.0
                      # do i=1,nn
                        for i in range(nn):
                          # he(nn1,i,ndd) = he(nn1,i,ndd) - x1*he(nnn1,i,ndd)
                            he[nn1,i,ndd] = he[nn1,i,ndd] - x1*he[nnn1,i,ndd]
                          # he(nn2,i,ndd) = he(nn2,i,ndd) - x2*he(nnn2,i,ndd)
                            he[nn2,i,ndd] = he[nn2,i,ndd] - x2*he[nnn2,i,ndd]
                      # enddo
                  # enddo
              # endif
              # x1 = 0.0d0
              # x2 = 0.0d0
                x1, x2 = 0., 0.
              # do i=1,nn
                for i in range(nn):
                  # x1 = x1 + he(nn1,i,ndd)**2
                    x1 = x1 + he[nn1,i,ndd]**2
                  # x2 = x2 + he(nn2,i,ndd)**2
                    x2 = x2 + he[nn2,i,ndd]**2
              # enddo
              # x1 = sqrt(0.5d0/(dx*x1))
                x1 = numpy.sqrt(0.5/(dx*x1))
              # x2 = sqrt(0.5d0/(dx*x2))
                x2 = numpy.sqrt(0.5/(dx*x2))
              # do i=1,nn
                for i in range(nn):
                  # he(nn1,i,ndd) = x1*he(nn1,i,ndd)
                    he[nn1,i,ndd] = x1*he[nn1,i,ndd]
                  # he(nn2,i,ndd) = x2*he(nn2,i,ndd)
                    he[nn2,i,ndd] = x2*he[nn2,i,ndd]
              # enddo
          # enddo
      # enddo

      # !c........................ building and diagonalization of the hamiltonian
      # !c                         storage of the eigenvectors and the eigenvalues
      # ho0 = (hox*hoy*hoz)**(1.0d0/3.0d0)
        ho0 = (hox*hoy*hoz)**(1.0 / 3.0)
        ax  = hox/ho0
        ay  = hoy/ho0
        az  = hoz/ho0
      # ho0 = ho0*hbm(it)
        ho0 = ho0*hbm[it]
      # ho0 = ho0*cf(it)
        ho0 = ho0*cf[it]
      # nw  = nwn
      # if (it.eq.2) nw = nwp
        nw = nwn if (it == 1) else \
             nwp
      # np  = npn
      # if (it.eq.2) np = npp
        np = npn if (it == 1) else \
             npp

      # !c................................................................ nucleus
      # x = xk(it)
        x = xk[it]
      # y = xmu(it)
        y = xmu[it]
      # if (numpy.gt.50) x = xk(it+2)
      # if (numpy.gt.50) y = xmu(it+2)
        if np > 50:
            x = xk[it+2]
            y = xmu[it+2]

      # !c..................................................... loop on the blocks
      # this is terrible code
      # ia = 0
        ia = -1 ### an inline incremented index
      # do 16 ni=1,nblc
      ### ni is an index and must therefor be 0-based
        for ni in range(nblc):
          # n  = ns(ni+1)-ns(ni)
            n  = ns[ni+1]-ns[ni]
          # nn = ns(ni)
            nn = ns[ni]
            # !c............................................... loop on the first vector
          # do 17 i=1,n
          ### i is an index and must therefor be 0-based
            for i in range(n):
              # nx1 = nx(nn+i)
                nx1 = nx[nn+i]
              # ny1 = ny(nn+i)
                ny1 = ny[nn+i]
              # nz1 = nz(nn+i)
                nz1 = nz[nn+i]
                nb  = nx1 + ny1 + nz1
              # x1  = 1 - 2*mod(nb-nz1,2)
                x1  = 1 - 2*  ((nb-nz1) % 2)
              # x2  = 1 - 2*mod(ny1,2)
                x2  = 1 - 2*   (ny1 % 2)
          # !c.............................................. loop on the second vector
              # do j=i,n
              ### j is an index and must therefor be 0-based
                for j in range(i,n):
                  # nx2 = nx(nn+j)
                    nx2 = nx[nn+j]
                  # ny2 = ny(nn+j)
                    ny2 = ny[nn+j]
                  # nz2 = nz(nn+j)
                    nz2 = nz[nn+j]
          # !c..................................... computation of the matrix elements
                  # if (j.ne.i) go to 19
                    if i == j:
                      # h(i,j) = ax*(nx1+0.5d0) + ay*(ny1+0.5d0) + az*(nz1+0.5d0) &
                      # &        -x*y*((nb*(nb+1))/2.0d0 -nx1**2 -ny1**2 -nz1**2)
                        h[i,j] = ax*(nx1+0.5) + ay*(ny1+0.5) + az*(nz1+0.) - x*y*( nb*(nb+1)*0.5 - nx1**2 -ny1**2 -nz1**2 )
                      # go to 18
                    else:
                      # 19 if (nx1.ne.nx2) go to 20
                        if nx1 == nx2:
                          # h(i,j) = 0.0d0
                            h[i,j] = 0.0
                            an = ny2
                            am = ny1
                          # if (ny2.eq.ny1+1) h(i,j) = x*sqrt(an*nz1)
                          # if (ny2.eq.ny1-1) h(i,j) = x*sqrt(am*nz2)
                          # if (ny2.eq.ny1+2) h(i,j) =-x*y*sqrt(an*(ny2-1)*nz1*(nz1-1))
                          # if (ny2.eq.ny1-2) h(i,j) =-x*y*sqrt(am*(ny1-1)*nz2*(nz2-1))
                            if ny2 == ny1 + 1:
                                h[i,j] = x * numpy.sqrt(an*nz1)
                            elif ny2 == ny1 - 1:
                                h[i,j] = x * numpy.sqrt(am*nz2)
                            elif ny2 == ny1 + 2:
                                h[i, j] = -x*y*numpy.sqrt(an*(ny2-1)*nz1*(nz1-1))
                            elif ny2 == ny1 - 2:
                                h[i, j] = -x*y*numpy.sqrt(am*(ny1-1)*nz2*(nz2-1))
                          # go to 18
                        else:
                       # 20 if (nz1.ne.nz2) go to 21
                            if nz1 == nz2:
                              # h(i,j) = 0.0d0
                                h[i,j] = 0.
                                am = nx1
                                an = nx2
                              # if (nx2.eq.nx1+1) h(i,j) =-x*x1*sqrt(an*ny1)
                              # if (nx2.eq.nx1-1) h(i,j) =-x*x1*sqrt(am*ny2)
                              # if (nx2.eq.nx1+2) h(i,j) =-x*y*sqrt(an*(nx2-1)*ny1*(ny1-1))
                              # if (nx2.eq.nx1-2) h(i,j) =-x*y*sqrt(am*(nx1-1)*ny2*(ny2-1))
                                if nx2 == nx1+1:
                                    h[i,j] = -x*x1*numpy.sqrt(an*ny1)
                                if nx2 == nx1-1:
                                    h[i,j] = -x*x1*numpy.sqrt(am*ny2)
                                if nx2 == nx1+2:
                                    h[i,j] = -x*y*numpy.sqrt(an*(nx2-1)*ny1*(ny1-1))
                                if nx2 == nx1-2:
                                    h[i,j] = -x*y*numpy.sqrt(am*(nx1-1)*ny2*(ny2-1))
                                # go to 18
                            else:
                           # 21 h(i,j) = 0.0d0
                                h[i,j] = 0.0
                              # if (ny1.ne.ny2) go to 18
                                if ny1 == ny2:
                                  # x3 = 1.0d0-2.0d0*mod(nx1,2)
                                    x3 = 1.0 - 2.0*(nx1 % 2)
                                    an = nx2
                                    am = nx1
                                  # x4 = 1.0d0-2.0d0*mod(nx2,2)
                                    x4 = 1.0 - 2.0 * (nx2 % 2)
                                  # if (nz2.eq.nz1+1) h(i,j) =-x*x2*x3*sqrt(am*nz2)
                                  # if (nz2.eq.nz1-1) h(i,j) =-x*x2*x4*sqrt(an*nz1)
                                  # if (nz2.eq.nz1+2) h(i,j) = x*y*sqrt(nz2*(nz2-1)*am*(nx1-1))
                                  # if (nz2.eq.nz1-2) h(i,j) = x*y*sqrt(nz1*(nz1-1)*an*(nx2-1))
                                    if nz2 == nz1+1:
                                        h[i,j] =-x*x2*x3 * numpy.sqrt(am*nz2)
                                    if nz2 == nz1-1:
                                        h[i,j] =-x*x2*x4 * numpy.sqrt(an*nz1)
                                    if nz2 == nz1+2:
                                        h[i,j] = x*y * numpy.sqrt(nz2*(nz2-1)*am*(nx1-1))
                                    if nz2 == nz1-2:
                                        h[i,j] = x*y * numpy.sqrt(nz1*(nz1-1)*an*(nx2-1))

                   # 18 h(j, i) = h(i, j)
                        h[j, i] = h[i, j]
                  # enddo
          # 17 continue
          # ! call diagon (h,ndim,n,s,d,wd, ifail)

          # ! Diagonalization in the subblock
          # # ! Inquire about the optimal size of work
          #   allocate(work(1)) ; lwork = -1
          #   call DSYEV( 'V', 'U', n, h(1:n,1:n), n, d, work, lwork, ifail)
          # # ! Change to the optimal value
          #   lwork = int(work(1)) ;  deallocate(work) ; allocate(work(lwork))
          # # ! Do the diagonalization
          #   call DSYEV( 'V', 'U', n, h(1:n,1:n), n, d(1:n), work, lwork, ifail)
          #   deallocate(work)
            eigenvalues, eigenvectors = numpy.linalg.eig(h[:n,:n])
          #   s(1:n,1:n)= h(1:n,1:n)
            s = h
          # !c.......................storage and shift of the single particle energies
          # irep(ni) = ia
            irep[ni] = ia
          # do i=1,n
            for i in range(n):
              # do j=1,n
                for j in range(n):
                    ia    = ia + 1
                  # a(ia) = s(i,j)
                    a[ia] = s[i,j]
              # enddo
              # e(nn+i) = d(i)*ho0 - 50.0
                e[nn+i] = d[i]*ho0 - 50.0
          # enddo
      # 16 continue

        # !if (it.eq.1) print 101,x,y,cb
        # !if (it.eq.2) print 102,x,y,ca,cb
        # !c.............. ordering the eigenvalues according to increasing energies
        # !c              nor(i) gives the original position of the i th s.p. energy
      # do i=1,nvec
        for i in range(nvec):
          # nor(i) = i
            nor[i] = i
      # enddo
        i1 = nvec - 1
      # do i=1,i1
        for i in range(i1):
            ii = i + 1
          # do j=ii,nvec
            for j in range(ii,nvec):
              # if (e(j).lt.e(i)) then
                if e[j] < e[i]:
                  # x = e(i)
                    x = e[i]
                  # e(i) = e(j)
                    e[i] = e[j]
                  # e(j) = x
                    e[j] = x
                  # k      = nor(i)
                    k = nor[i]
                  # nor(i) = nor(j)
                    nor[i] = nor[j]
                  # nor(j) = k
                    nor[j] = k
              # endif
          # enddo
      # enddo
      # do i=1,nvec
        for i in range(nvec):
          # j = nor(i)
            j = nor[i]
          # k = nx(j) + ny(j) + nz(j)
            k = nx[j] + ny[j] + nz[j]
          # npa(i) = 1 - 2*mod(k,2)
            npa[i] = 1 - 2*(k % 2)
      # enddo
      #   !print 103,(i,nor(i),e(i),npa(i),i=1,nvec)
      #   !c..................................... selection of the nw wave-functions
      #   !c       different filling is obtained by previous change of the array nor
        j = 0
      # do i=1,nw
        for i in range(nw):
          # if (npa(i).ne.-1) then
            if npa[i] != -1:
                j = j + 1
              # ntrs(j) = i
                ntrs[j] = i
          # endif
      # enddo
      # npar(1,it) = j
        npar[0,it] = j
      # do i=1,nw
        for i in range(nw):
          # if (npa(i).ne.1) then
            if npa[i] != -1:
                j = j + 1
              # ntrs(j) = i
                ntrs[j] = i
          # endif
      # enddo

      # npar(2,it) = j - npar(1,it)
        npar[2,it] = j - npar[1,it]

      # do iwave=1,nw
        for iwave in range(nw):
            nwave = nwave + 1
          # if (iwave.le.npar(1,it)) go to 49
            if iwave > npar[1,it]:
              # kparz(nwave) =-1
                kparz[nwave] = -1
              # i = iwave - npar(1,it)
                i = iwave - npar[1,it]
              # go to 50
            else:
              # 49 kparz(nwave) =+1
                kparz[nwave] =-1
          # 50 i = ntrs(iwave)
            i = ntrs[iwave]
          # esp1(nwave) = e(i)
            esp1[nwave] = e[i]
          # j = nor(i)
            j = nor[i]
          # do nn=1,nblc
            for nn in range(nblc):
              # n = ns(nn+1) - ns(nn)
                n = ns[nn+1] - ns[nn]
              # ia = irep(nn)
                ia = irep[nn]
              # do i=1,n
                for i in range(n):
                  # do ja=1,n
                    for ja in range(n):
                        ia = ia + 1
                      # s(i,ja) = a(ia)
                        s[i,ja] = a[ia]
                  # enddo
              # enddo
                nvv = nn
              # if (j.le.ns(nn+1)) go to 45
                if j <= ns[nn+1]: break
          # enddo
          # 45 nn  = ns(nvv)
            nn = ns[nvv]

          # if(allocated(spwf_map)) then
            if spwf_map is not None:
              # ! Only construct the spwfs if the user asks for it
              # !         => allocated status of spwf_map

              # ! we have constructed all spwfs for this particular MPI rank
              # if(store_counter.gt.size(spwf_map)) cycle
                if store_counter > spwf_map.size: continue

              # ! This is an spwf we want to store
              # if(spwf_map(store_counter) .eq. nwave) then
                if spwf_map[store_counter] == nwave:
                  # ! construct the spwf in the array psi
                  # psi = 0.0d0
                    psi = .0

                    ny2 = 0
                    kk  = 0
                  # do k=1,4
                    for k in range(4): ### 0 based indexing
                      # TODO: check that nvv is 0-based
                      # if (nsi(nvv,k).eq.0) go to 46
                        if nsi[nvv,k] != 0:
                            nx2 = ny2 + 1
                          # ny2 = ny2 + nsi(nvv,k)
                            ny2 = ny2 + nsi[nvv,k]
                          # nz2 = 0
                          # if (k.eq.1.or.k.eq.3) nz2=3
                            nz = 0 if k % 2 == 0 else \
                                 3
                          # do i=nx2,ny2
                            for i in range(nx2,ny2):
                              # nx1 = nx(nn+i) + 1
                              # ny1 = ny(nn+i) + 1
                              # nz1 = nz(nn+i) + 1
                              # TODO: chech the removal of + 1
                                nx1 = nx[nn+i]
                                ny1 = ny[nn+i]
                                nz1 = nz[nn+i]
                              # xph = s(i,j-nn)
                                xph = s[i,j-nn]

                              # if (mod(ny1,4).eq.nz2) xph =-xph
                                if ny1 % 4 == nz2: xph *= -1.
                              # do ix=1,mx
                                for ix in range(mx):
                                  # hex = he(nx1,ix,1)
                                    hex = he[nx1, ix, 0]
                                  # do iy=1,my
                                    for iy in range(my):
                                      # hey = he(ny1,iy,2)
                                        hey = he[ny1,iy,1]
                                      # do iz=1,mz
                                        for iz in range(mz):
                                          # hez = he(nz1,iz,3)
                                            hez = he[nz1,iz,2]
                                          # psi(ix,iy,iz,k) = psi(ix,iy,iz,k) + xph*hex*hey*hez
                                            psi[ix,iy,iz,k] = psi[ix,iy,iz,k] + xph*hex*hey*hez
                                      # enddo
                                  # enddo
                              # enddo
                          # enddo

                      # 46 kk = kk + mz
                        kk = kk + mz
                  # enddo
                  # ! ... and copy it. ### another terrible idea!
                    ind = 0
                  # do k = 1,mz
                    for k in range(mz):
                      # do j = 1,my
                        for j in range(my):
                          # do i = 1,mx
                            for i in range(mx):
                              # ind = ind + 1 ### zero based indexing
                              # wfs(ind,1,store_counter) = psi(i,j,k,1)
                                wfs[ind,0,store_counter] = psi[i,j,k,1]
                              # wfs(ind,2,store_counter) = psi(i,j,k,2)
                                wfs[ind,1,store_counter] = psi[i,j,k,2]
                              # wfs(ind,3,store_counter) = psi(i,j,k,3)
                                wfs[ind,2,store_counter] = psi[i,j,k,3]
                              # wfs(ind,4,store_counter) = psi(i,j,k,4)
                                wfs[ind,3,store_counter] = psi[i,j,k,4]
                          # enddo
                      # enddo
                  # enddo
                    store_counter = store_counter + 1
              # endif
          # endif
      # enddo
  # enddo

  # deallocate(h,s,d,wd)
  # deallocate(nsi,ns)
  # deallocate(nx,ny, nz, e, nor,npa)
  # deallocate(he, a)
  # deallocate(irep, ntrs)

# end subroutine nilsson
# end module nil8

