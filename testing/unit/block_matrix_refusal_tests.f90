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
program block_matrix_refusal_tests
!==============================================================================
  !
  ! Refusal tests for the structure guards of the symmetry-block matrix
  !  module's rotation machinery (design section 6a). Every scenario hands
  !  the module an input it must REFUSE loudly (print + stop 1) - the same
  !  loud convention as the guards of init_sym_matrix, which cannot be
  !  caught inside a running process either.
  !
  ! The program executes exactly ONE scenario, selected by the
  !  command-line argument: it states what it expects the module to refuse,
  !  makes the offending call, and - if the module failed to stop -
  !  complains and exits with status 2. The run script drives the
  !  scenarios and requires exit status 1 from each.
  !
  ! The guards under test:
  !   1: init_from_dense refuses a C that mixes sectors
  !   2: init_from_dense refuses a non-orthonormal sector-diagonal C
  !   3: init_from_dense refuses a wrongly sized C
  !   4: set_sector refuses a non-orthonormal block
  !   5: set_sector refuses a wrongly sized block
  !   6: set_sector refuses a block for an empty sector
  !   7: set_sector refuses to act on an uninitialised rotation
  !   8: rotate refuses a rotation built from different HFBlocks
  !   9: full_trace refuses the TR-folded composite (D-7, 2026-11-20)
  !==============================================================================
  use block_matrix_module

  implicit none

  integer :: scenario, ios
  character(len=32) :: arg
  integer :: hf(8), pp(8), i
  real(8), allocatable :: C(:, :), Ublk(:, :)
  complex(8) :: tr
  type(sym_rotation) :: R
  type(sym_matrix) :: M

  scenario = 0
  if (command_argument_count() > 0) then
    call get_command_argument(1, arg)
    read (arg, *, iostat=ios) scenario
    if (ios .ne. 0) scenario = 0
  endif

  ! One basis for the scenarios: MOCCa blocks 1 and 3 populated
  !  (sectors 0 and 2, counts 4 and 3), the companion pattern pp counts
  !  the same total over sector 1 instead of 2.
  hf = [4, 0, 3, 0, 0, 0, 0, 0]
  pp = [4, 3, 0, 0, 0, 0, 0, 0]

  select case (scenario)
  case (1)
    print *, 'expect refusal: init_from_dense with sector-mixing C'
    allocate(C(7, 7))
    C = 0.0d0
    do i = 1, 7
      C(i, i) = 1.0d0
    enddo
    C(1, 7) = 0.5d0
    C(7, 1) = 0.5d0
    call R%init_from_dense(hf, C)
  case (2)
    print *, 'expect refusal: init_from_dense with non-orthonormal blocks'
    allocate(C(7, 7))
    C = 0.0d0
    do i = 1, 7
      C(i, i) = 2.0d0
    enddo
    call R%init_from_dense(hf, C)
  case (3)
    print *, 'expect refusal: init_from_dense with a wrongly sized C'
    allocate(C(5, 5))
    C = 0.0d0
    do i = 1, 5
      C(i, i) = 1.0d0
    enddo
    call R%init_from_dense(hf, C)
  case (4)
    print *, 'expect refusal: set_sector with a non-orthonormal block'
    call R%init(hf)
    allocate(Ublk(4, 4))
    Ublk = 0.0d0
    do i = 1, 4
      Ublk(i, i) = 2.0d0
    enddo
    call R%set_sector(0, Ublk)
  case (5)
    print *, 'expect refusal: set_sector with a wrongly sized block'
    call R%init(hf)
    allocate(Ublk(3, 3))
    Ublk = 0.0d0
    do i = 1, 3
      Ublk(i, i) = 1.0d0
    enddo
    call R%set_sector(0, Ublk)
  case (6)
    print *, 'expect refusal: set_sector for an empty sector'
    call R%init(hf)
    allocate(Ublk(2, 2))
    Ublk = 0.0d0
    do i = 1, 2
      Ublk(i, i) = 1.0d0
    enddo
    call R%set_sector(1, Ublk)
  case (7)
    print *, 'expect refusal: set_sector on an uninitialised rotation'
    allocate(Ublk(4, 4))
    Ublk = 0.0d0
    do i = 1, 4
      Ublk(i, i) = 1.0d0
    enddo
    call R%set_sector(0, Ublk)
  case (8)
    print *, 'expect refusal: rotate with mismatched HFBlocks'
    call R%init(hf)
    call M%init(pp, 0, HERMITIAN, STY_EVEN, TAU_EVEN)
    call M%rotate(R)
  case (9)
    print *, 'expect refusal: full_trace on the TR-folded composite'
    call M%init(hf, 0, HERMITIAN, STY_EVEN, TAU_EVEN, folded=.true.)
    tr = M%full_trace()
  case default
    print *, 'usage: block_matrix_refusal_tests <scenario 1..9>'
    stop 2
  end select

  print *, 'FAILURE: scenario', scenario, 'was not refused'
  stop 2
end program block_matrix_refusal_tests
