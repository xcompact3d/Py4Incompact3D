!!! py4incompact3d.f90 provides light-weight interface to
!!!                    2decomp&fft/Xcompact3d for wrapping with f2py.
!!! Copyright (C) 2021  University of Edinburgh
!!!
!!! This program is free software: you can redistribute it and/or modify
!!! it under the terms of the GNU General Public License as published by
!!! the Free Software Foundation, either version 3 of the License, or
!!! (at your option) any later version.
!!!
!!! This program is distributed in the hope that it will be useful,
!!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!! GNU General Public License for more details.
!!!
!!! You should have received a copy of the GNU General Public License
!!! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module py4incompact3d
  ! A module to wrap the 2decomp&fft library such that it can be wrapped by f2py.
  
  use mpi
  use decomp_2d
  
  implicit none

  private
  public :: init_py4incompact3d

contains

  subroutine init_py4incompact3d(nx, ny, nz, p_row, p_col)
    ! Initialises the 2decomp&fft library for a domain discretised by :math:`n_x \times n_y \times
    ! n_z` mesh points with a :math:`p_{row} \times p_{col}` pencil decomposition.
    
    integer, intent(in) :: nx, ny, nz
    integer, intent(in) :: p_row, p_col

    integer :: ierr

    ! XXX: Need to initialise the nproc and nrank variables from 2decomp&fft.
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    
    call decomp_2d_init(nx, ny, nz, p_row, p_col)
    
  end subroutine init_py4incompact3d
  
end module py4incompact3d
