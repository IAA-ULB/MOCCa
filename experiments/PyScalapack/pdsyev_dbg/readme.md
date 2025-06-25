# Solve an eigenvalue/eigenvector problem with a real symmetric matrix.

we start small on an 8x8 matrix with 1x1 processes.

[pdsyev](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/p-syev.html) 

the reason why we took so long to get this working is that i thought that a ScaLAPACK array descriptor in PyScalapack tranlates to `*array.scalapack_params()`. That happened to be wrong! 

`array.scalapack_params()` returns a 4-tuple:

``` python
   def scalapack_params(self):
        """
        Get parameters use to pass to scalapack. Scalapack usually take matrix as parameter with the following format.
        The pointer to data, the row and column indices in the global matrix indicating the start point of all local matrix,
        and array desc.

        Returns
        -------
        args : tuple
            The arguments to be passed to scalapack function.
        """
        return (
            Scalapack.numpy_ptr(self.data),
            Scalapack.one,
            Scalapack.one,
            self,
        )
```

Hence, in the fortran function call:

    CALL PDSYEV ('V','L',n,A_sub,1,1,A_sub_dsc,eigenvalues,eigenvectors,1,1,A_sub_dsc, work,lwork,info)

the 4 parameters `A_sub,1,1,A_sub_dsc` are to be replaced with `*A_sub.scalapack_params()` in a PyScalasca call!

# Note!

The integers `1,1`, which ar fixed by `scalapack_params()`  represent:

    ia, ja
    (global) The row and column indices in the array A indicating the first row and the first column, respectively, of the submatrix of A `to copy`. 1 ≤ia≤total_rows_in_a - m +1, 1 ≤ja≤total_columns_in_a - n +1.

This is usually correct, unless A is part of a larger data array. 