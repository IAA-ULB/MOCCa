# Things I struggled with

## Numroc

    na_rows = numroc(na, nblk, my_prow, 0, np_rows)
    na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

    na: dimension of the full matrix (in the direction considered:rows or columns)
    nblk: block size of distributed matrix (in the block-cyclic sense).
    my_prow: row index of the process in the process grid
    np_rows: number of processes in a process row of the process grid
    my_pcol: column index of the process in the process grid
    np_cols: number of processes in a process column of the process grid
    
## Setting up a blacs context

### without MPI

There is no explicit mpi_init() for the user. ! BLACS_PINFO does the mpi_init() (see https://www.netlib.org/scalapack/explore-html/d6/dd3/blacs__pinfo___8c_a14adb195fd2f7ffd95e208c3a5d332fb.html)

``` fortran    
    integer :: rank, nranks, ctxt_sys, ctxt_all, nprow, npcol, myrow, mycol
    call BLACS_PINFO( rank, nranks )
    call BLACS_GET( 0, 0, ctxt_sys )
    ctxt_all = ctxt_sys
    call BLACS_GRIDINIT( ctxt_all, 'C', 2, 2)
    call BLACS_GRIDINFO( ctxt_all, nprow,npcol,myrow,mycol)
```

### with MPI

Here MPI is initialized by the calling program 

``` fortran    
    integer :: info, mpi_comm_world, rank, nranks, ictxt, nprow, npcol, myrow, mycol
    call mpi_init(info)
    call mpi_comm_rank(mpi_comm_world,rank,info)
    call mpi_comm_size(mpi_comm_world,nranks,info)  
    ictxt = mpi_comm_world
    call blacs_get(-1, 0, ictxt)
    call blacs_gridinit(ictxt, ’C’, np_rows, np_cols)
    call blacs_gridinfo(ictxt, np_rows, np_cols, my_prow, my_pcol)
```
So, if we want to use PyScalapack with mpi4py, we must extend PyScalapack to use this scheme.

# Local and global indices

functions `indxl2g` and `indxg2l` deal with Fortran indices: `1<=iF<=n`. In PyScalapack you are typically dealing with Python indices: `0<=iP<n` and you must convert: `iF = iP +1`