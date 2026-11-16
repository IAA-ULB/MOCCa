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
module block_matrix_module
 !==============================================================================
 !
 ! This module contains the machinery to efficiently represent matrices
 !  in the configuration space spanned by a finite set of single-particle
 !  wavefunctions on the MOCCa grid.
 !
 ! Minimal placeholder type, to be elaborated per the design summary
 ! 'mocca-symmetry-block-matrix-design-summary.md', section 2: the 3-bit
 ! symmetry word of an operator. bit 0: parity (0 = even, 1 = odd);
 ! bit 1: z-signature (0 = +1, 1 = -1); bit 2: isospin (0 = isoscalar,
 ! 1 = flips proton/neutron).

 type operator_symmetry
    integer :: word = 0         ! 3-bit sector word as per the design summary
 end type operator_symmetry

 type block_matrix

    integer :: conjugation_class           ! real or imaginary
    integer :: hermeticity                 ! hermitian or antihermitian
    type(operator_symmetry) :: op_symmetry ! quantum numbers - P, Rz, T
    real(8), allocatable :: blocks(:) ! Flattened storage for all non-zero blocks

    integer :: total_size            ! Total matrix size (e.g., 1000)
    integer :: block_counts(2, 2, 2) = 0  ! Number of basis states per (P, S_z, T_z)
    integer :: block_offsets(2, 2, 2) = 0 ! Starting index for each block in the global basis

    logical :: is_nonzero(2, 2, 2, 2, 2, 2) = .false. ! Which (P_i, S_i, T_i, P_j, S_j, T_j) blocks are non-zero?
    integer :: block_indices(2, 2, 2, 2, 2, 2) = 0 ! Mapping to flattened storage
!  contains
!    procedure :: init => init_block_matrix
!    procedure :: construct => construct_block_matrix
!    procedure :: set_block => set_block_matrix
!    procedure :: get_block => get_block_matrix
!    procedure :: matvec => block_matrix_matvec
  end type block_matrix

end module block_matrix_module
