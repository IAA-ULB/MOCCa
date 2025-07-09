# Things I struggled with

## Numroc

```fortran
    na_rows = numroc(na, nblk, my_prow, 0, np_rows)
    na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

    ! na: dimension of the full matrix (in the direction considered:rows or columns)
    ! nblk: block size of distributed matrix (in the block-cyclic sense).
    ! my_prow: row index of the process in the process grid
    ! np_rows: number of processes in a process row of the process grid
    ! my_pcol: column index of the process in the process grid
    ! np_cols: number of processes in a process column of the process grid
```
    
# Local and global indices

functions `indxl2g` and `indxg2l` deal with Fortran indices: `1<=iF<=n`. In PyScalapack you are typically dealing with Python indices: `0<=iP<n` and you must convert: `iF = iP +1`.

## Setting up a blacs context

Can we combine `MPI` and `BLACS`? In Fortran? In Python? Apparently, yes, but it took me about a week to find out how ...

### without MPI

There is no explicit `mpi_init()` for the user. `BLACS_PINFO()` does the mpi_init() (see [here](https://www.netlib.org/scalapack/explore-html/d6/dd3/blacs__pinfo___8c_a14adb195fd2f7ffd95e208c3a5d332fb.html)).

``` fortran    
    integer :: rank, nranks, ctxt_sys, ctxt_all, nprow, npcol, myrow, mycol
    call BLACS_PINFO( rank, nranks )
    call BLACS_GET( 0, 0, ctxt_sys )
    ctxt_all = ctxt_sys
    call BLACS_GRIDINIT( ctxt_all, 'C', 2, 2)
    call BLACS_GRIDINFO( ctxt_all, nprow,npcol,myrow,mycol)
```

### with MPI

Here MPI is explicitly initialized by the calling program, which thereby gains access to MPI_COMM_WORLD. This allows the calling program to do other MPI calls. 

``` fortran    
    integer :: info, mpi_comm_world, rank, nranks, ictxt, nprow, npcol, myrow, mycol
    call mpi_init(info)
    call mpi_comm_rank(mpi_comm_world,rank,info)
    call mpi_comm_size(mpi_comm_world,nranks,info)  
    call blacs_get(-1, 0, ictxt)
    call blacs_gridinit(ictxt, ’C’, np_rows, np_cols)
    call blacs_gridinfo(ictxt, np_rows, np_cols, my_prow, my_pcol)
```

`BLACS_GET` picks up the initialized MPI. 

I also found out that the BLACS context variables in fact hold the MPI communicator. And that `MPI_WORLD_COMM` is available via `BLACS_GET`, even when `mpi_init` is not called by the main program.

What happens if we call `MPI_INIT` ánd `BLACS_PINFO`? Apparently, this does not seem to be a problem. So, it looks like `BLACS_PINFO` picks up MPI_COMM_WORLD if `mpi_init` was already called, and initializes MPI itself. if `mpi_init` was not yet called.

This made me wonder why we had problems to combine `mpi4py` and `PyScalapack`. `from mpi4py import MPI` automatically initializes MPI and exposes `MPI.COMM_WORLD`, but this was not picked up by `BLACS`. Finally, it came to my mind that the `mpi4py` module we were loading (`mpi4py/4.0.1-gompi-2024a`) is linked with a different `libmpi.so` than the `iimkl` module.

    /apps/antwerpen/zen2/rocky8/OpenMPI/5.0.3-GCC-13.3.0/lib: libmpi.so (for mpi4py)
    /apps/antwerpen/zen2/rocky8/impi/2021.13.0-intel-compilers-2024.2.0/mpi/2021.13/lib: libmpi.so (for iimlk)

Obviously, `BLACS` cannot pick up the MPI initialization by `mpi4py` as this happen in a `.so` file to which it has no access.

Atfer Franky built `mpi4py/4.0.1-iimpi-2024a` which links mpi4py to the same `libmpi.so` file as `iimkl` the problem was gone.

