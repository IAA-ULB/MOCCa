# extending PyScalapack with ELPA

How would that work? Elpa uses Fortran object orientation ... no clue how to wrap that with ctypes.

It gets hard from the beginning. adding the libelpa.so to the pyscalapack makes mpi crash:


``` Python
scalapack = pyscalapack(
    ...
    "/apps/antwerpen/skylake/rocky9/ELPA/2024.05.001-intel-2024a/lib/libelpa.so"
)  
```

    ===================================================================================
    =   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
    =   RANK 0 PID 3916298 RUNNING AT login9.breniac
    =   KILLED BY SIGNAL: 9 (Killed)
    ===================================================================================

It might be the case that the pyscalapack `lib*.so` files and `libelpa.so` contain their own blacs version... and maybe what else is duplicated. 

So let's try to use only `libelpa.so`. Weirdly enough, the program crashes with the same error, but yield the correct solution both with `numpy.linalg.eigh` and `pdsyev`. That probably means that `libelpa.so` has its own ScaLAPACK version. That problem is solved by calling `blacs_exit(0)` at the end of the program.

``` Python
if __name__ == "__main__":
    ...
    pyscalapack.blacs_exit(0)
```
It would be nice it that could be done automatically.

So, we can now start finding out how to wrap the elpa object

#### interesting links
- [fortran classes](https://cyber.dabamos.de/programming/modernfortran/object-oriented-programming.html)
- [numpy.f2py](https://numpy.org/doc/stable/f2py/index.html#f2py) not sure if it can handle fortran classes
- [interfacing python with fortran](https://www-uxsup.csx.cam.ac.uk/courses/moved.PythonFortran/f2py.pdf)
- [fmodpy](https://github.com/tchlux/fmodpy) can maybe handle fortran classes
- [Which package should I use to wrap Modern Fortran Code with Python?](https://scicomp.stackexchange.com/questions/2283/which-package-should-i-use-to-wrap-modern-fortran-code-with-python)
- [Fortran Wiki - Object-oriented programming](https://fortranwiki.org/fortran/show/Object-oriented+programming) 
- [Types from Fortran to Python via Opaque Pointers](https://rgoswami.me/posts/types-fortran-python-opaque/)
- [How to Call Fortran from Python](https://www.matecdev.com/posts/fortran-in-python.html)

[Fortran Wiki - Object-oriented programming](https://fortranwiki.org/fortran/show/Object-oriented+programming) has useful info about object-orientation in Modern Fortran. It seems as if these two ways exist:  

``` fortran
    object%member_function(args)
    member_function(object, args)
```

we can find the functions exported by a `.so` file as 

    nm libelpa.so

So maybe, it is possible to access theses functions with the pyscalapack ctypes interface.

After some trial and error in `experiments/elpa/sandbox` we learned this:

This is the way elpa exposes allocation and deallocation of Elpa objects. The allocator returns a `c_ptr` to the allocated object:

``` fortran
    ! in file elpa_impl.F90 line 347
    function elpa_impl_allocate_c(error) result(ptr) bind(C, name="elpa_allocate")
      integer(kind=c_int)        :: error
      type(c_ptr)                :: ptr
      type(elpa_impl_t), pointer :: obj

      obj => elpa_impl_allocate(error)
      ptr = c_loc(obj)
    end function
```

The deallocator accepts one, which it converts to a type(circle) pointer to free it: 

``` fortran
    ! in file elpa_impl.F90 line 382
    subroutine elpa_impl_deallocate_c1(handle) bind(C, name="elpa_deallocate1")
      type(c_ptr), value         :: handle
      type(elpa_impl_t), pointer :: self

      call c_f_pointer(handle, self)
      call self%destroy()
      deallocate(self)
    end subroutine
```
The manipulators (setters) of elpa_objects do the same:

``` fortran
    subroutine elpa_set_integer_c(handle, name_p, value, error) bind(C, name="elpa_set_integer")
      type(c_ptr), intent(in), value                :: handle
      type(elpa_impl_t), pointer                    :: self
      type(c_ptr), intent(in), value                :: name_p
      character(len=elpa_strlen_c(name_p)), pointer :: name
      integer(kind=c_int), intent(in), value        :: value
      integer(kind=c_int) , intent(in)              :: error

      call c_f_pointer(handle, self)
      call c_f_pointer(name_p, name)
      call elpa_set_integer(self, name, value, error)
    end subroutine
```

`libelpa.so` exports the following procedures dealing with elpa objects:

```
...
0000000000094280 T elpa_allocate
...
00000000000942b0 T elpa_deallocate
...
0000000000096fb0 T elpa_set_double
0000000000096b80 T elpa_set_float
00000000000964f0 T elpa_set_integer
...
```

At this point (fingers crossed) it seems that we do not need to build a fortran `elpa_object_interface.F90` so file to access and manipulate `elpa` objects.