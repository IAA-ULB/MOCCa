!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
!    Copyright (C) 2026 W. Ryssens and M. Bender
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Affero General Public License as published
!    by the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!===============================================================================
program mocca_single
  !-----------------------------------------------------------------------------
  ! This is a very simple driver module for the entire code; it  parses any
  ! command line arguments (if there are any) and then passes control to
  ! Run_MOCCa.
  !-----------------------------------------------------------------------------
  use MOCCa

  implicit none

  !-----------------------------------------------------------------------------
  ! Dealing with the optional command line arguments of the code
  ! You can invoke the code as
  !
  !  ./MOCCa.exe filename file_number
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
  integer           :: Narguments, status
  integer(dp)       :: file_number = 10
  character(len=32) :: filename = 'input.dat', numberstring

  Narguments = COMMAND_ARGUMENT_COUNT()
  if(Narguments .eq. 0) then
    ! Run the code from STDIN
    call Run_MOCCa()
  else
    ! Read filename
    call get_command_argument(1,filename,status=status)
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
    call Run_MOCCa(file_number, filename)
  endif
end program 
