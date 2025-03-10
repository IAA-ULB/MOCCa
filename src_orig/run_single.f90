program tantalus_single
  !-----------------------------------------------------------------------------
  ! This is a very simple driver module for the entire code; it
  ! parses any command line arguments (if there are any) and
  ! then passes control to Run_Tantalus.
  !
  !-----------------------------------------------------------------------------
  use Tantalus
  use FAM

  implicit none

  !-----------------------------------------------------------------------------
  ! Dealing with the optional command line arguments of the code
  ! You can invoke the code as
  !
  !  ./Tantalus.exe filename file_number
  !
  !  filename    = string indicating the name of the file that
  !                contains the input data
  !  file_number = integer indicating the channel to open the
  !                input file on
  !
  ! Both these arguments are optional:
  ! - if filename is absent, the code reads STDIN for input
  ! - if file_number is absent, the code will use a default value
  !------------------------------------------------------------------------------
  !
  integer           :: Narguments, status
  integer(dp)       :: file_number = 10
  character(len=32) :: filename = 'input.dat', numberstring
  logical           :: do_fam = .true.


  Narguments = COMMAND_ARGUMENT_COUNT()
  print *, 'Nargs', Narguments
  if(Narguments .eq. 0) then
    ! Run the code from STDIN
    ! call Run_Tantalus('Single-mode')

    if(do_fam) then
      ! Run a FAM QRPA calculation
      call Run_FAM()
    endif
  else
    ! Read filename
    call get_command_argument(1,filename,status=status)
    print *, 'filename', filename
    if(status.gt.0) then
      call stp('Unknown error when reading the first command line argument.')
    elseif(status.eq.-1) then
      call stp('Filename is too long.')
    endif

    if(Narguments.gt.1) then
      ! Override the filenumber if provided
      call get_command_argument(2,numberstring,status=status)
      if(status.ne.0) then
        call stp('Unknown error when reading the second command line argument.')
      endif
      READ(numberstring, "(i10)") file_number
      if(file_number .eq. 12) then
        print *, 'FILENUMBER 12 is reserved for parameterization reading.'
        print *, 'Please choose a different file_number.'
        call stp('')
      endif
    endif
    ! Run the code from input on file "filename"
    call Run_Tantalus('Single-mode', file_number, filename)
  endif
end program 
