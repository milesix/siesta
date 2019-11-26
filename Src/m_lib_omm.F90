! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_lib_omm

use MatrixSwitch

use parallel,       only : BlockSize, Node, Nodes, ionode
use precision,      only : dp
use sys,            only : die
#ifdef MPI
use mpi_siesta,     only : mpi_comm_world
use parallelsubs,   only : set_BlockSizeDefault
#endif

implicit none

!**** PUBLIC ************************************!

public :: omm_min

!************************************************!

contains

!================================================!
! use the orbital minimization method (OMM) from libOMM to   !
! solve the eigenvalue problem (double precision !
! routine, for Gamma point-only calculations)    !
! not parallel yet
!================================================!
subroutine omm_min(CalcE,PreviousCallDiagon,iscf,istp,nbasis,nspin,h_dim,nhmax,numh,listhptr,listh,d_sparse,eta,qs,h_sparse,&
    s_sparse,t_sparse)
  

  implicit none

  !**** INPUT ***********************************!

  logical, intent(in) :: CalcE              ! Calculate the energy-density matrix from the existing coeffs.?
  logical, intent(in) :: PreviousCallDiagon ! Previous SCF iteration solved by diagonalization?

  integer, intent(in) :: iscf               ! SCF iteration num.
  integer, intent(in) :: istp               ! MD iteration num.
  integer, intent(in) :: nbasis             ! dimension of numh and listhptr
  integer, intent(in) :: nspin              ! num. of spins
  integer, intent(in) :: h_dim              ! num. of AOs (global)
  integer, intent(in) :: nhmax              ! first dimension of listh and sparse matrices
  integer, intent(in) :: numh(1:nbasis)     ! num. of nonzero elements of each row of sparse matrices
  integer, intent(in) :: listhptr(1:nbasis) ! pointer to start of row in listh
  integer, intent(in) :: listh(1:nhmax)     ! list of nonzero elements of each row of sparse matrices

  real(dp), intent(in) :: qs(1:2)                             ! num. of electrons per spin
  real(dp), intent(in) :: eta(1:2)                            ! chemical potential for Kim functional
  real(dp), intent(in), optional :: h_sparse(1:nhmax,1:nspin) ! hamiltonian matrix (sparse)
  real(dp), intent(in), optional :: t_sparse(1:nhmax)         ! kinetic energy matrix (sparse)
  real(dp), intent(in), optional :: s_sparse(1:nhmax)         ! overlap matrix (sparse)

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: d_sparse(1:nhmax,1:nspin) ! (energy-)density matrix (sparse)

  !**** LOCAL ***********************************!

  character(5) :: m_storage
  character(3) :: m_operation

  logical :: init_C

  integer :: mpi_rank
#ifdef MPI
  integer :: mpi_err, mpi_size, BlockSize_c_default
#endif

  integer :: i, j, io, jo, ind, k, l, N_occ
  real(dp) :: he, se, e_min
  type(matrix) :: H, S, D_min, T, C_min  
  integer, allocatable, save :: h_dim_l2g(:) ! local-to-global index transform for AOs
  real(dp), allocatable, save :: C_old(:,:)
  logical, save :: first_call = .true.

  !**********************************************!
  
  if (ionode) then
    print'(a)','OMM with libOMM'
  end if

  if (nspin == 1) then
    N_occ = nint(0.5_dp*qs(1))
  else
    N_occ = nint(qs(1))
  end if
  
#ifdef MPI
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

! calculate the ScaLAPACK blocking factor for distributing the WF coeffs. matrix
  call set_blocksizedefault(Nodes,N_occ,BlockSize_c_default)

  if (first_call) call ms_scalapack_setup(mpi_size,1,'c',BlockSize_c_default)

  m_storage ='pddbc'
  m_operation ='lap'
#else
  mpi_rank = 0
  m_storage = 'sdden'
  m_operation = 'ref'
#endif
  
  if(first_call) then    
    allocate(h_dim_l2g(1:nbasis))
    h_dim_l2g(:) = 0
    allocate(C_old(1:h_dim, 1:nbasis))
    C_old(:,:) = 0.0_dp

#ifdef MPI
    j = 0
    k = 0
    l = 0
    do i = 1, h_dim
      k = k+1
      if (j == Node) then
        l = l+1
        h_dim_l2g(l) = i
      end if
      if (k == BlockSize) then
        k = 0
        j = j+1
        if (j == Nodes) j = 0
      end if
    end do
#else
    do i = 1, nbasis
      h_dim_l2g(i) = i
    end do
#endif
  end if

  if (.not. H%is_initialized) call m_allocate(H,h_dim,h_dim,m_storage)
  if (.not. S%is_initialized) call m_allocate(S,h_dim,h_dim,m_storage)
  if (.not. D_min%is_initialized) call m_allocate(D_min,h_dim,h_dim,m_storage)
  if (.not. C_min%is_initialized) call m_allocate(C_min,N_occ,h_dim,m_storage)

  init_C = .true.  

  call m_set(H,m_operation,0.0_dp,0.0_dp)
  call m_set(S,m_operation,0.0_dp,0.0_dp)

  do io = 1, nbasis
    do j = 1, numh(io)
      ind = listhptr(io)+j
      jo = listh(ind)
      he = h_sparse(ind,1)
      se = s_sparse(ind)
      call m_set_element(H,io,jo,he,0.0_dp)
      call m_set_element(S,io,jo,se,0.0_dp)
    end do
  end do 
   
  call m_set(D_min,m_operation,0.0_dp,0.0_dp)  
  call m_set(C_min,m_operation,0.0_dp,0.0_dp)

  do io = 1, N_occ
    do jo = 1, nbasis
      ind = h_dim_l2g(jo)
      call m_set_element(C_min,io,ind,C_old(io,jo),0.0_dp)
    end do
  end do

  if(.not. calcE) then
    if(first_call) init_C = .false.
    call omm(h_dim,N_occ,H,S,.true.,e_min,D_min,.false.,0.0_dp,&
      C_min,init_C,T,0.0_dp,0,nspin,1,-1.0_dp,.true.,.true.,&
      m_storage,m_operation,mpi_rank)
    C_old(:,:) = 0.0_dp
    do io = 1, N_occ
      do jo = 1, nbasis
        ind = h_dim_l2g(jo)
        call m_get_element(C_min,io,ind,C_old(io,jo))
      end do
    end do
  else
    call omm(h_dim,N_occ,H,S,.true.,e_min,D_min,.false.,0.0_dp,&
      C_min,init_C,T,0.0_dp,0,nspin,1,-1.0_dp,.true.,.false.,&
      m_storage,m_operation,mpi_rank)
    call omm(h_dim,N_occ,H,S,.false.,e_min,D_min,.true.,0.0_dp,&
      C_min,init_C,T,0.0_dp,0,nspin,1,-1.0_dp,.true.,.true.,&
      m_storage,m_operation,mpi_rank)
  end if	  
  if (ionode) then
    print'(a, f13.7)','e_min = ', e_min
  end if
  d_sparse(:,:) = 0.0_dp
  do io = 1, nbasis
    do j = 1, numh(io)
      ind = listhptr(io) + j
      jo = listh(ind)
      call m_get_element(D_min,io,jo,d_sparse(ind,1))
      if(nspin == 1) d_sparse(ind,1) = 2.0_dp*d_sparse(ind,1)
      if(nspin == 2) d_sparse(ind,2) = d_sparse(ind,1)
    end do
  end do

  if(first_call) first_call = .false.
   
  call m_deallocate(D_min)
  call m_deallocate(S)
  call m_deallocate(H)
  call m_deallocate(C_min)

end subroutine omm_min

end module m_lib_omm
