! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
SUBROUTINE timer_mpi( name, opt )
  use mpi_siesta, timer_mpi_module => timer_mpi
  character(len=*), intent(in):: name
  integer,          intent(in):: opt

  call timer_mpi_module( name, opt )

END SUBROUTINE timer_mpi

