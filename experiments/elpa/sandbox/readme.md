Folder `experiments/elpa/sandbox` contains a simple example to understand how we can have a Fortran class with dynamic memory allocation that is wrapped in Pythnn using ctypes. The ultimate goal is to build such a wrapper around Elpa.

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

Warning:

``` python
p_circle = circleso.alloc() # p_circle holds a raw pointer
print(f"{p_circle=}") # YIELDS A SEGMENTATION FAULT!
print(f"{ctypes.c_void_p(p_circle)=}") # fine!
```

On the command line we observed that de first test after rebuilding the so file is usually ok. Rerunning the test often yields a segmentation fault. In a new terminal the first test always worked, and in a batch job generally too 

- unless the so file was missing because another job was rebuilding it
- unless - I guess - because two jobs were running simultaneously
  
So it seems as if the problems are related to caching of so files and not to the code.