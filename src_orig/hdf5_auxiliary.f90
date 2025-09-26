module HDF5_auxiliary 
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
 ! Module providing auxiliary routines for HDF5 I/O. 
 ! Notes: 
 !  - this file needs no preprocessing by Hephaestos
 !  - this file does not have a preprocessor directive because the Makefile
 !    ensures that it is only compiled when HDF5 is enabled.
 !
 ! Most of this module is the work of Nikolai Shchechilin.
 ! 
 !
 ! TODO: 
 ! - refactor the code to overload subroutines instead of having separate
 !   routines for different data types
 ! - look into the possibility of writing scalar attributes as 1d attributes
 !   of lentgh 0/1 to reduce code duplication
 ! - refactor the different ways of error reporting
 !==============================================================================

  use HDF5
  use geninfo, only : dp, stp

  implicit none

  ! The compression level of compression for hdf5
  !  - only active for potentials and densities, not for spwfs
  !  - level 6 seems to be optimal
  integer, parameter  :: comprlvl = 6

  ! NOTE: there is currently NO interface that overloads hdf5_write_dataset_1d
  !       because that would invalidate the use of hdf5_write_dataset_1d for 
  !       real arrays of arbitrary rank.
  !interface hdf5_write_dataset_1d
  !    module procedure hdf5_write_dataset_1d_real
  !    module procedure hdf5_write_dataset_1d_integer
  !end interface hdf5_write_dataset_1d

contains 
  
  subroutine hdf5_write_attr_char(id, name, attribute)
    !----------------------------------------------------------------------------
    ! This routine writes a character scalar attribute with to an open hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    ! attribute : string, value of the attribute
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    integer(hid_t), intent(in) :: id
    integer(HID_T)             :: space_id, attribute_id !identifiers
    integer(HID_T)             :: type_id !identifiers
    integer                    :: error
    integer(size_t)            :: alen
    character(len=*), intent(in) :: name, attribute
    Integer(hsize_t), dimension (1) :: dims

    dims(1)=1
    alen=len(attribute)
    if(alen.eq.0) call stp('ERROR: zero length attribute in hdf5_write_attr_char')
    !creating datatype
    call h5tcopy_f(h5t_native_character, type_id, error)
    call h5tset_size_f(type_id, alen, error)
    !create space
    call h5screate_f(h5s_scalar_f, space_id, error)
    !create attribute
    call h5acreate_f(id, name, type_id, space_id, attribute_id, error)
    !write attribute
    call h5awrite_f(attribute_id, type_id, attribute, dims, error)
    !close attribute
    call h5aclose_f(attribute_id, error)
    !close space
    call h5sclose_f(space_id, error)
    !close type
    call h5tclose_f(type_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'writing')
  end subroutine hdf5_write_attr_char

  subroutine hdf5_write_attr_char_1d(id, name, attribute, n, alen)
    !----------------------------------------------------------------------------
    ! This routine writes a character array attribute with to an open hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    ! attribute : string, value of the attribute
    ! n         : integer, length of the array
    ! alen      : integer, length of each string in the array
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    integer(hid_t), intent(in) :: id
    integer, intent(in)        :: n
    integer, intent(in)        :: alen
    character(len=*), intent(in)    :: name
    character(len=alen), intent(in) :: attribute(n)
    Integer(hsize_t), dimension (1) :: dims
    integer                        :: error
    integer(HID_T)             :: space_id, attribute_id, type_id, alen_hdf5 

    dims(1)=n
    alen_hdf5 = int(alen, kind=hid_t) ! type conversion for hdf5
    if(alen.eq.0) call stp('ERROR: zero length attribute in hdf5_write_attr_char_1d')
    !creating datatype
    call h5tcopy_f(h5t_native_character, type_id, error)
    call h5tset_size_f(type_id, alen_hdf5, error)
    !create space
    call h5screate_f(h5s_simple_f, space_id, error)
    !create attribute
    call h5acreate_f(id, name, type_id, space_id, attribute_id, error)
    !write attribute
    call h5awrite_f(attribute_id, type_id, attribute, dims, error)
    !close attribute
    call h5aclose_f(attribute_id, error)
    !close space
    call h5sclose_f(space_id, error)
    !close type
    call h5tclose_f(type_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'writing')
  end subroutine hdf5_write_attr_char_1d

  subroutine hdf5_write_attr_integer(id, name, attribute)
    !----------------------------------------------------------------------------
    ! This routine writes integer scalar attribute in an open HDF5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    ! attribute : integer, value of the attribute
    ! 
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    integer(hid_t), intent(in) :: id
    integer(hid_t)             :: space_id, attribute_id !identifiers
    integer, intent(in)        :: attribute
    integer                    :: error
    character(len=*), intent(in) :: name
    Integer(size_t), dimension (1)        :: dims=(/0/)

    !create space
    call h5screate_f(h5s_scalar_f, space_id, error)
    !create attribute
    call h5acreate_f(id, name, H5T_Native_Integer, space_id, attribute_id, error)
    !write attribute
    call h5awrite_f(attribute_id, H5T_Native_Integer, attribute, dims, error)
    !close attribute
    call h5aclose_f(attribute_id, error)
    !close space
    call h5sclose_f(space_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'writing')
  end subroutine hdf5_write_attr_integer
  
  subroutine hdf5_write_attr_integer_1d(id, name, attribute, n)
    !---------------------------------------------------------------------------- 
    ! This routine writes integer array attribute to an open hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    ! attribute : integer array, value of the attribute
    ! n         : integer, length of the array
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    integer(hid_t), intent(in) :: id
    integer(hid_t)             :: space_id, attribute_id !identifiers
    integer, intent(in)        :: n
    integer, intent(in)        :: attribute(n)
    integer                    :: error
    character(len=*), intent(in) :: name
    Integer(size_t), dimension (1) :: dims
    
    dims=(/n/)
    !create space
    call h5screate_simple_f(1, dims, space_id, error)
    !create attribute
    call h5acreate_f(id, name, H5T_NATIVE_INTEGER, space_id, attribute_id, error)
    !write attribute
    call h5awrite_f(attribute_id, H5T_NATIVE_INTEGER, attribute, dims, error)
    !close attribute
    call h5aclose_f(attribute_id, error)
    !close space
    call h5sclose_f(space_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'writing')
  end subroutine hdf5_write_attr_integer_1d

  subroutine hdf5_write_attr_double(id, name, attribute)
    !----------------------------------------------------------------------------
    ! This routine writes a double precision scalar attribute to an open hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    ! attribute : double precision, value of the attribute
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    integer(hid_t), intent(in) :: id
    integer(hid_t)             :: space_id, attribute_id !identifiers
    real(kind=dp), intent(in)  :: attribute
    integer                    :: error
    character(len=*), intent(in) :: name
    Integer(size_t), dimension (1)        :: dims=(/0/)

    !create space
    call h5screate_f(h5s_scalar_f, space_id, error)
    !create attribute
    call h5acreate_f(id, name, H5T_Native_Double, space_id, attribute_id, error)
    !write attribute
    call h5awrite_f(attribute_id, H5T_Native_Double, attribute, dims, error)
    !close attribute
    call h5aclose_f(attribute_id, error)
    !close space
    call h5sclose_f(space_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'writing')
  end subroutine hdf5_write_attr_double

  subroutine hdf5_write_attr_double_1d(id, name, attribute, n)
    !----------------------------------------------------------------------------
    ! This routine writes a double precision scalar attribute to an open hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    ! attribute : double precision, value of the attribute
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    integer(hid_t), intent(in) :: id
    integer(hid_t)             :: space_id, attribute_id !identifiers
    real(kind=dp), intent(in)  :: attribute(n)
    integer, intent(in)        :: n
    integer                    :: error
    character(len=*), intent(in) :: name
    Integer(size_t), dimension (1) :: dims

    dims=(/n/)
    !create space
    call h5screate_f(h5s_scalar_f, space_id, error)
    !create attribute
    call h5acreate_f(id, name, H5T_Native_Double, space_id, attribute_id, error)
    !write attribute
    call h5awrite_f(attribute_id, H5T_Native_Double, attribute, dims, error)
    !close attribute
    call h5aclose_f(attribute_id, error)
    !close space
    call h5sclose_f(space_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'writing')
  end subroutine hdf5_write_attr_double_1d

  subroutine hdf5_write_dataset_1d(id, name, dset, n, groupname)
    !----------------------------------------------------------------------------
    ! writes double precision dataset array of length n with some name in the hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group)
    ! name      : string, name of the dataset; will get trimmed for superfluous spaces
    ! dset      : double precision array, value of the dataset
    ! n         : integer, length of the array
    ! groupname : optional string, name of the group to write the dataset in 
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: groupname
    integer(hid_t), intent(in) :: id
    integer,        intent(in) :: n
    real(kind=dp),  intent(in) :: dset(n)
    integer(hid_t)             :: space_id, dset_id, plist_id !identifiers
    integer                    :: error
    integer(hsize_t), dimension(1) :: dims,data_dims

    dims=(/n/)
    data_dims=(/n/)
    ! Create dataspace for data_set 
    call h5screate_simple_f(1, dims, space_id, error)
    ! create property list
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    if(error.ne.0) call stp('ERROR: h5pcreate_f failed in hdf5_write_dataset_1d for name='//trim(name))
    ! create chunks with property list for compression, as of now size of chunk      
    ! is just equal to the size of array. Modify for MPI reading?
    call h5pset_chunk_f(plist_id, 1, dims, error)
    if(error.ne.0) call stp('ERROR: h5pset_chunk_f failed in hdf5_write_dataset_1d for name='//trim(name))
    ! shuffling for better compression?
    call h5pset_shuffle_f(plist_id, error)
    if(error.ne.0) call stp('ERROR: h5pset_shuffle_f failed in hdf5_write_dataset_1d for name='//trim(name))
    ! zlib compression with deflate
    call h5pset_deflate_f(plist_id, comprlvl, error)
    if(error.ne.0) call stp('ERROR: h5pset_deflate_f failed in hdf5_write_dataset_1d for name='//trim(name))
    ! Create dataset with default properties "dset_id" is returned
    if(present(groupname)) then
      call h5dcreate_f(id,groupname//'/'//trim(name),H5T_NATIVE_DOUBLE,space_id,dset_id,error,plist_id)
    else 
      call h5dcreate_f(id,trim(name),H5T_NATIVE_DOUBLE,space_id,dset_id,error,plist_id)    
    endif
    ! Write dataset 
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dset, data_dims, error)
    if(error.ne.0) call stp('ERROR: h5dwrite_f failed in hdf5_write_dataset_1d for name='//trim(name))
    ! Close access to dataset 
    call h5dclose_f(dset_id, error)
    ! Close access to data space 
    call h5sclose_f(space_id, error)
    ! close access to plist
    call h5pclose_f(plist_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'writing')
  end subroutine hdf5_write_dataset_1d

  subroutine hdf5_write_dataset_1d_integer(id, name, dset, n, groupname)
    !----------------------------------------------------------------------------
    ! writes double precision dataset array of length n with some name in the hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group)
    ! name      : string, name of the dataset; will get trimmed for superfluous spaces
    ! dset      : integer array, value of the dataset
    ! n         : integer, length of the array
    ! groupname : optional string, name of the group to write the dataset in 
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: groupname
    integer(hid_t), intent(in) :: id
    integer,        intent(in) :: n
    integer      ,  intent(in) :: dset(n)
    integer(hid_t)             :: space_id, dset_id
    integer                    :: error
    integer(hsize_t), dimension(1) :: dims,data_dims

    dims=(/n/)
    data_dims=(/n/)
    ! Create dataspace for data_set 
    call h5screate_simple_f(1, dims, space_id, error)
    if(present(groupname)) then
      call h5dcreate_f(id,groupname//'/'//trim(name),H5T_NATIVE_INTEGER,space_id,dset_id,error)
    else 
      call h5dcreate_f(id,trim(name),H5T_NATIVE_INTEGER,space_id,dset_id,error)    
    endif
    ! Write dataset 
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dset, data_dims, error)
    if(error.ne.0) call stp('ERROR: h5dwrite_f failed in hdf5_write_dataset_1d_integer for name='//trim(name))
    ! Close access to dataset 
    call h5dclose_f(dset_id, error)
    ! Close access to data space 
    call h5sclose_f(space_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'writing')
  end subroutine hdf5_write_dataset_1d_integer

  subroutine hdf5_write_dataset_2d(id, name, dset, n1, n2, groupname)
    !----------------------------------------------------------------------------
    ! writes double precision dataset array of length n with some name in the hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group)
    ! name      : string, name of the dataset; will get trimmed for superfluous spaces
    ! dset      : double precision array, value of the dataset
    ! n         : integer, length of the array
    ! groupname : optional string, name of the group to write the dataset in 
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: groupname
    integer(hid_t), intent(in) :: id
    integer,        intent(in) :: n1,n2
    real(kind=dp),  intent(in) :: dset(n1,n2)
    integer(hid_t)             :: space_id, dset_id, plist_id !identifiers
    integer                    :: error
    integer(hsize_t), dimension(2) :: dims,data_dims

    dims=(/n1, n2/)
    data_dims=(/n1, n2/)
    ! Create dataspace for data_set 
    call h5screate_simple_f(2, dims, space_id, error)
    if(error.ne.0) call stp('ERROR: h5screate_simple_f failed in hdf5_write_dataset_2d')
    ! create property list
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    if(error.ne.0) call stp('ERROR: h5pcreate_f failed in hdf5_write_dataset_2d')
    ! Create dataset with default properties "dset_id" is returned
    if(present(groupname)) then
      call h5dcreate_f(id,groupname//'/'//trim(name),H5T_NATIVE_DOUBLE,space_id,dset_id,error,plist_id)
    else 
      call h5dcreate_f(id,trim(name),H5T_NATIVE_DOUBLE,space_id,dset_id,error,plist_id)    
    endif
    if(error.ne.0) call stp('ERROR: h5dcreate_f failed in hdf5_write_dataset_2d')
    ! Write dataset 
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dset, data_dims, error)
    if(error.ne.0) call stp('ERROR: h5dwrite_f failed in hdf5_write_dataset_2d')
    ! Close access to dataset 
    call h5dclose_f(dset_id, error)
    ! Close access to data space 
    call h5sclose_f(space_id, error)
    ! close access to plist
    call h5pclose_f(plist_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'writing')
  end subroutine hdf5_write_dataset_2d

  subroutine hdf5_read_attr_char(id, name, attribute, n)
    !----------------------------------------------------------------------------
    ! reads character scalar attribute with some name from the hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    ! n         : integer, length of the character string
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    integer(hid_t), intent(in) :: id
    integer(HID_T)             :: attribute_id !identifiers
    integer(HID_T)             :: type_id !identifiers
    integer                    :: error
    integer(size_t)            :: n
    character(len=*), intent(in) :: name
    character(len=*), intent(inout) :: attribute
    Integer(hsize_t), dimension (1) :: dims

    dims(1)=n
    !open attribute
    call h5aopen_name_f(id, name, attribute_id, error)
    !get the type
    call h5aget_type_f(attribute_id, type_id, error)
    !read attribute
    call h5aread_f(attribute_id, type_id, attribute, dims, error)
    !close the attribute
    call h5aclose_f(attribute_id,error)

    if (error.ne.0) call report_hdf5_error(id, name, 'reading')
  end subroutine hdf5_read_attr_char
  
  subroutine hdf5_read_attr_integer(id, name, attribute)
    !----------------------------------------------------------------------------
    ! reads integer scalar attribute with some name from the hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    integer(hid_t), intent(in) :: id
    integer(HID_T)             :: attribute_id !identifiers
    integer                    :: error
    character(len=*), intent(in) :: name
    integer, intent(out)         :: attribute
    Integer(size_t), dimension (1) :: dims

    dims(1)=1
    !open attribute
    call h5aopen_name_f(id, name, attribute_id, error)
    !read attribute
    call h5aread_f(attribute_id, H5T_Native_Integer, attribute, dims, error)
    !close the attribute
    call h5aclose_f(attribute_id,error)

    if (error.ne.0) call report_hdf5_error(id, name, 'reading')
  end subroutine hdf5_read_attr_integer
  
  subroutine hdf5_read_attr_integer_1d(id, name, attribute, n)
    !----------------------------------------------------------------------------
    ! reads integer array attribute with some name from the hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    ! n         : integer, length of the array
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------

    integer(hid_t), intent(in) :: id
    integer, intent(in)        :: n
    integer(HID_T)             :: attribute_id !identifiers
    integer                    :: error
    character(len=*), intent(in) :: name
    integer, intent(out)         :: attribute(n)
    Integer(size_t), dimension (1) :: dims

    dims(1)=n
    !open attribute
    call h5aopen_name_f(id, name, attribute_id, error)
    !read attribute
    call h5aread_f(attribute_id, H5T_Native_Integer, attribute, dims, error)
    !close the attribute
    call h5aclose_f(attribute_id,error)

    if (error.ne.0) call report_hdf5_error(id, name, 'reading')
  end subroutine hdf5_read_attr_integer_1d

  subroutine hdf5_read_attr_double(id, name, attribute)
    !----------------------------------------------------------------------------
    ! reads double precision scalar attribute with some name from the hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    integer(hid_t), intent(in) :: id
    integer(HID_T)             :: attribute_id !identifiers
    integer                    :: error
    character(len=*), intent(in) :: name
    real(kind=dp), intent(out) :: attribute
    integer(size_t), dimension (1) :: dims

    dims(1)=1
    !open attribute
    call h5aopen_name_f(id, name, attribute_id, error)
    !read attribute
    call h5aread_f(attribute_id, H5T_Native_Double, attribute, dims, error)
    !close the attribute
    call h5aclose_f(attribute_id,error)

    if (error.ne.0) call report_hdf5_error(id, name, 'reading')
  end subroutine hdf5_read_attr_double

  subroutine hdf5_read_attr_double_1d(id, name, attribute, n)
    !----------------------------------------------------------------------------
    ! reads double precision scalar attribute with some name from the hdf5 file
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the attribute
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    integer(hid_t), intent(in)   :: id
    character(len=*), intent(in) :: name
    real(kind=dp), intent(out)   :: attribute(n)
    integer, intent(in)          :: n

    integer(HID_T)             :: attribute_id !identifiers
    integer                    :: error
    integer(size_t), dimension (1) :: dims

    dims(1)=n
    !open attribute
    call h5aopen_name_f(id, name, attribute_id, error)
    !read attribute
    call h5aread_f(attribute_id, H5T_Native_Double, attribute, dims, error)
    !close the attribute
    call h5aclose_f(attribute_id,error)

    if (error.ne.0) call report_hdf5_error(id, name, 'reading')
  end subroutine hdf5_read_attr_double_1d
    
  subroutine hdf5_read_dataset_1d(id, name, dset, n, groupname)
    !----------------------------------------------------------------------------
    ! reads double precision dataset array of length n with some name in the hdf5 file
    !
    ! Input:
    ! id   : hid_t,  identifier of the hdf5 object (file, group)
    ! name : string, name of the dataset
    ! dset : double precision array, value of the dataset
    ! n    : integer, length of the array
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: id
    integer, intent(in)          :: n
    character(len=*), intent(in), optional :: groupname
    real(kind=dp), intent(inout) :: dset(n)
    integer(hid_t)               :: dset_id !identifiers
    integer(size_t), dimension(1):: dims, data_dims
    integer                      :: error
    dims=(/n/)
    data_dims(1)=n
    ! open dataset, "dset_id" is returned
    if(present(groupname)) then
      call h5dopen_f(id,groupname//'/'//trim(name), dset_id, error)
    else 
      call h5dopen_f(id,trim(name), dset_id, error)
    endif
    ! read dataset 
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dset, data_dims, error)
    ! Close access to dataset 
    call h5dclose_f(dset_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'reading')

  end subroutine hdf5_read_dataset_1d

  subroutine hdf5_read_dataset_1d_integer(id, name, dset, n, groupname  )
    !----------------------------------------------------------------------------
    ! reads integer dataset array of length n with some name in the hdf5 file
    !
    ! Input:
    ! id   : hid_t,  identifier of the hdf5 object (file, group)
    ! name : string, name of the dataset
    ! dset : integer, value of the dataset
    ! n    : integer, length of the array
    ! groupname : optional string, name of the group to read the dataset from
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    character(len=*), intent(in)           :: name
    integer(hid_t), intent(in)             :: id
    integer, intent(in)                    :: n
    integer, intent(out)                   :: dset(n)
    character(len=*), intent(in), optional :: groupname
    integer(hid_t)               :: dset_id !identifiers
    integer(size_t), dimension(1):: dims, data_dims
    integer                      :: error
    dims=(/n/)
    data_dims(1)=n
    ! open dataset, "dset_id" is returned
    if(present(groupname)) then
      call h5dopen_f(id,groupname//'/'//trim(name), dset_id, error)
    else 
      call h5dopen_f(id,trim(name), dset_id, error)
    endif
    ! read dataset 
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, dset, data_dims, error)
    ! Close access to dataset 
    call h5dclose_f(dset_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'reading')

  end subroutine hdf5_read_dataset_1d_integer

  subroutine hdf5_read_dataset_2d(id, name, dset, n1,n2)
    !----------------------------------------------------------------------------
    ! reads double precision dataset rank 2 ALLOCATABLE array of dimension 
    ! (n1,n2) with some name in the hdf5 file
    !
    ! Input:
    ! id   : hid_t,  identifier of the hdf5 object (file, group)
    ! name : string, name of the dataset
    ! dset : double precision array, value of the dataset
    ! n1   : integer, length of the array
    ! n2   : integer, length of the array
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)              :: id
    integer, intent(in)                     :: n1,n2
    real(kind=dp), intent(out), allocatable :: dset(:,:)
    integer(hid_t)               :: dset_id !identifiers
    integer(size_t), dimension(2):: dims, data_dims
    integer                      :: error
    dims     =(/n1,n2/)
    data_dims=(/n1,n2/)

    allocate(dset(n1,n2)) ; dset = 0.0d0
    ! open dataset, "dset_id" is returned
    call h5dopen_f(id, name, dset_id, error)
    if(error.ne.0) call stp('ERROR: h5dopen_f failed in hdf5_read_dataset_2d')
    ! read dataset 
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dset, data_dims, error)
    if(error.ne.0) call stp('ERROR: h5dread_f failed in hdf5_read_dataset_2d')
    ! Close access to dataset 
    call h5dclose_f(dset_id, error)

    if (error.ne.0) call report_hdf5_error(id, name, 'reading')

  end subroutine hdf5_read_dataset_2d

  subroutine report_hdf5_error(id, name, operation)
    !----------------------------------------------------------------------------
    ! This routine reports an HDF5 error with some additional information
    !
    ! Input:
    ! id        : hid_t,  identifier of the hdf5 object (file, group, dataset)
    ! name      : string, name of the object
    ! operation : string, "reading" or "writing".
    !
    ! Output:
    ! none
    !----------------------------------------------------------------------------
    use geninfo, only : stp

    1 format ('ERROR: '//a7//' HDF5 object '//A20//' with identifier'//I5//'.')

    integer(hid_t), intent(in)    :: id
    integer                       :: error
    character(len=*), intent(in)  :: name
    character(len=7), intent(in)  :: operation
    CHARACTER(len=70)             :: msg

    write(msg,fmt=1) trim(operation), trim(name), id 
    call stp(msg)
  end subroutine report_hdf5_error
end module hdf5_auxiliary