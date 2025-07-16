module elpa_object_interface

  use elpa 
  use iso_c_binding

  implicit none

  class(elpa_t), pointer           :: e
  integer  :: initialised = -1

  contains

  subroutine initialize() bind(C, name='elpa_initialize_')
  
    if (elpa_init(20240501) /= elpa_ok) then
      print *, "ELPA API version not supported"
      stop 1
    endif
    
    e => elpa_allocate()
    initialised = 0
  end subroutine

  subroutine finalize() bind(C, name='elpa_finalize_')

    if (initialised == 0) then
      call elpa_deallocate(e)
      call elpa_uninit()
    endif
  
  end subroutine

  subroutine set_integer(name_p, value, error) bind(C, name='set_integer_')

      type(c_ptr)        , intent(in), value   :: name_p
      integer(kind=c_int), intent(in), value   :: value
      integer            , intent(inout)       :: error ! on input the length of the string name_p, output 0 (succes) or error code.
      character(len=error)           , pointer :: name_f
      
      call c_f_pointer(name_p, name_f)
      ! write (*,*) 'set_integer ',len(name_f), '<<',name_f,'>>'
      call e%set(name_f, value, error)
      ! write (*,*) 'set_integer error=', error
      
  end subroutine

  subroutine set_double(name_p, value, error) bind(C, name='set_double_')

      type(c_ptr)           , intent(in), value   :: name_p
      integer(kind=c_double), intent(in), value   :: value
      integer               , intent(inout)       :: error ! on input the length of the string name_p, output 0 (succes) or error code.
      character(len=error)              , pointer :: name_f
      
      ! call c_f_pointer(name_p, name_f)
      write (*,*) 'set_integer ',len(name_f), '<<',name_f,'>>'
      ! call e%set(name_f, value, error)
      write (*,*) 'set_integer error=', error
      
  end subroutine

  subroutine set_float(name_p, value, error) bind(C, name='set_float_')

      type(c_ptr)          , intent(in), value   :: name_p
      integer(kind=c_float), intent(in), value   :: value
      integer              , intent(inout)       :: error ! on input the length of the string name_p, output 0 (succes) or error code.
      character(len=error)             , pointer :: name_f
      
      ! call c_f_pointer(name_p, name_f)
      write (*,*) 'set_integer ',len(name_f), '<<',name_f,'>>'
      ! call e%set(name_f, value, error)
      write (*,*) 'set_integer error=', error
      
  end subroutine
 

end module elpa_object_interface