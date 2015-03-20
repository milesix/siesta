! ---
! Copyright (C) 1996-2014	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
MODULE sparse_matrices
  use precision

  implicit none

  ! Max. nonzero elements in a row of the hamiltonian matrix
  integer:: nh 

  ! Max. number of nonzero H matrix elements    
  integer :: maxnh = 10 



  integer, pointer :: listh(:), listhptr(:),  numh(:)

  real(dp), pointer :: Dold(:,:), Dscf(:,:), Eold(:,:), &
                       Escf(:,:), H(:,:)
  real(dp), pointer :: xijo(:,:)

  real(dp), pointer :: H0(:), S(:)

END MODULE sparse_matrices
