!!! py4incompact3d.f90 provides light-weight interface to
!!!                    2decomp&fft/Xcompact3d for wrapping with f2py.
!!! Copyright (C) 2021  University of Edinburgh
!!!
!!! Licensed under the Apache License, Version 2.0 (the "License");
!!! you may not use this file except in compliance with the License.
!!! You may obtain a copy of the License at
!!!
!!!     http://www.apache.org/licenses/LICENSE-2.0
!!!
!!! Unless required by applicable law or agreed to in writing, software
!!! distributed under the License is distributed on an "AS IS" BASIS,
!!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!! See the License for the specific language governing permissions and
!!! limitations under the License.

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
