! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_lib_omm

use MatrixSwitch
use omm_rand
use libomm

use atomlist,       only : qa, lasto
use fdf,            only : fdf_integer
use parallel,       only : BlockSize, Node, Nodes, ionode
use precision,      only : dp
use sys,            only : die
#ifdef MPI
use mpi_siesta
use parallelsubs,   only : set_BlockSizeDefault
#endif

implicit none

!**** PUBLIC ************************************!

public :: omm_min, omm_min_block

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

#ifdef MPI
  integer :: MPIerror
  integer, allocatable :: ind_o(:)
  integer :: jo1, kk, ib_r, ib_c, i_node, ind1 
#endif
  
  integer :: i, j, io, jo, ind, N_occ
  real(dp) :: he, se, e_min, de
  type(matrix), save :: H, S, D_min, T, C_min  
  logical, save :: first_call = .true.
  logical :: found
  !**********************************************!

  if (nspin == 1) then
    N_occ = nint(0.5_dp*qs(1))
  else
    N_occ = nint(qs(1))
  end if
  
#ifdef MPI
  if (first_call) call ms_scalapack_setup(mpi_comm_world,1,'c',BlockSize)
  m_storage ='pddbc'
  m_operation ='lap'
#else
  m_storage = 'sdden'
  m_operation = 'ref'
#endif
  
  if(first_call) then
    if(ionode) print'(a)','OMM with libOMM'
    if (.not. H%is_initialized) call m_allocate(H,h_dim,h_dim,m_storage)
    if (.not. S%is_initialized) call m_allocate(S,h_dim,h_dim,m_storage)
    if (.not. D_min%is_initialized) call m_allocate(D_min,h_dim,h_dim,m_storage)
    if (.not. C_min%is_initialized) call m_allocate(C_min,N_occ,h_dim,m_storage)
  end if

  call m_set(H,'a',0.0_dp,0.0_dp,m_operation)
  call m_set(S,'a',0.0_dp,0.0_dp,m_operation)

#ifdef MPI
  allocate(ind_o(1:nhmax))
  do io = 1, nbasis
    ind_o(listhptr(io) + 1) = listhptr(io) + 1
    do j = 2, numh(io)
      ind = listhptr(io) + j
      ind_o(ind) = ind
      jo = listh(ind)
      if(jo .lt. listh(ind_o(ind - 1))) then
        do kk = 1, j - 1
          ind1 = listhptr(io) + j - kk
          jo1 = listh(ind_o(ind1))
          if(jo1 .gt. jo) then
            ind_o(ind1 + 1) = ind_o(ind1)
            ind_o(ind1) = ind
          end if
        end do
      end if
    end do
  end do

  ind = 1
  do ib_c = 1, h_dim
    do ib_r = 1, h_dim
      he = 0.0_dp 
      se = 0.0_dp
      found = .false.
      i_node = int((ib_c - 1) / BlockSize)
      if(Node == i_node) then
        io = mod((ib_c - 1), BlockSize) + 1
        if((io .le. nbasis) .and. (ind .le. (listhptr(io) + numh(io))) .and. (ind .ge. (listhptr(io) + 1))) then
          jo = listh(ind_o(ind))
          if(jo == ib_r) then
            found = .true.
            he = h_sparse(ind_o(ind), 1)
            se = s_sparse(ind_o(ind))
            ind = ind + 1           
          end if
        end if
      end if
      call MPI_Bcast(found,1,MPI_Logical,i_node,MPI_Comm_World,MPIerror)
      if(found) then
        call MPI_Bcast(he,1,MPI_Double_Precision,i_node,MPI_Comm_World,MPIerror)
        call MPI_Bcast(se,1,MPI_Double_Precision,i_node,MPI_Comm_World,MPIerror)
        call m_set_element(H,ib_r,ib_c,he,0.0_dp)
        call m_set_element(S,ib_r,ib_c,se,0.0_dp)
      end if
    end do
  end do 
#else
  do io = 1, nbasis
    do j = 1, numh(io)
      ind = listhptr(io) + j
      jo = listh(ind)
      he = h_sparse(ind, 1)
      se = s_sparse(ind)
      call m_set_element(H,jo,io,he,0.0_dp, m_operation)
      call m_set_element(S,jo,io,se,0.0_dp, m_operation)
    end do
  end do
#endif
  
  call m_set(D_min,m_operation,0.0_dp,0.0_dp)  

  if(first_call) then
    call m_set(C_min,m_operation,0.0_dp,0.0_dp)
  end if


  if(.not. calcE) then
    init_C = .true.
    if(first_call) init_C = .false. 
    call omm(h_dim,N_occ,H,S,.true.,e_min,D_min,.false.,0.0_dp,&
      C_min,init_C,T,0.0_dp,0,nspin,1,-1.0_dp,.true.,.false.,&
      m_storage,m_operation)
    if(ionode) print'(a, f13.7)','e_min = ', e_min
  else
    call omm(h_dim,N_occ,H,S,.false.,e_min,D_min,.true.,0.0_dp,&
      C_min,init_C,T,0.0_dp,0,nspin,1,-1.0_dp,.true.,.true.,&
      m_storage,m_operation)
  end if

#ifdef MPI
  ind = 1
  do ib_c = 1, h_dim
    do ib_r = 1, h_dim
      found = .false.
      i_node = int((ib_c - 1) / BlockSize)
      if(Node == i_node) then
        io = mod((ib_c - 1), BlockSize) + 1
        if((io .le. nbasis) .and. (ind .le. (listhptr(io) + numh(io))) .and. (ind .ge. (listhptr(io) + 1))) then
          jo = listh(ind_o(ind))
          if(jo == ib_r) found = .true.
        end if
      end if 
      call MPI_Bcast(found,1,MPI_Logical,i_node,MPI_Comm_World,MPIerror)
      if(found) then
        call m_get_element(D_min,ib_r,ib_c,de)
        if(Node == i_node) then
          d_sparse(ind_o(ind), 1) = 2.0 * de 
          if(nspin == 2) d_sparse(ind_o(ind), 2) = d_sparse(ind_o(ind), 1)
          ind = ind + 1          
        end if
      end if
    end do
  end do
  deallocate(ind_o)
#else
  do io = 1, nbasis
    do j = 1, numh(io)
      ind = listhptr(io) + j
      jo = listh(ind)
      call m_get_element(D_min,jo,io,de)
      d_sparse(ind, 1) = 2.0 * de
      if(nspin == 2) d_sparse(ind, 2) = d_sparse(ind, 1)
    end do
  end do
#endif

  if(first_call) first_call = .false.

end subroutine omm_min

subroutine omm_min_block(CalcE,PreviousCallDiagon,iscf,istp,nbasis,nspin,h_dim,nhmax,numh,listhptr,listh,d_sparse,&
    eta,qs,h_sparse,s_sparse,t_sparse)
  
  use neighbour,      only : jan, mneighb
  use siesta_geom,    only : na_u, xa, ucell
  use siesta_options, only : rcoor
 
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

#ifdef MPI
  integer :: MPIerror, MPI_Size
#endif
  integer, dimension(2) :: dims
  integer, dimension(:), pointer, save :: row_blk_sizes, col_blk_sizes, row_blk_sizes1, col_blk_sizes1
  integer, save :: nblocks_r, nblocks_c, nblocks_r1, nblocks_c1
  integer :: ind_c, ind_r, ib_r, ib_c, ind_c2, ind1, jo1, ib, blk_c, blk_r, num_c, ia_off, ib_co, ind_co
  integer :: i, j, io, jo, ind, k, l, N_occ, kk, i_node, BlockSize_c
  integer :: blk_co, blk_ro, blk_f
  integer :: index, neib0, nna, ind_r1, ind_r2, ind_a1, ind_a2, ia, ib_r1, ib_r2, indexi, iorb, norb, ja
  integer, allocatable :: ind_o(:), ind_u(:), neib(:)
  real(dp) :: he, se, e_min
  real(dp), allocatable :: block_data(:,:), block_data_s(:,:)
  real(dp),pointer :: myblock(:,:)
  type(matrix), save :: H, S, D_min, T, C_min  
  logical, save :: first_call = .true.
  logical :: found
  logical, allocatable :: found_c(:)
  integer :: seed
  real(dp) :: rn(2), el
  real(dp) :: rmax,rr(3),rrmod, cgval

  !**********************************************!
  if (nspin == 1) then
    N_occ = nint(0.5_dp*qs(1))
  else
    N_occ = nint(qs(1))
  end if
   
#ifdef MPI
  dims(:) = 2
  if (first_call) then 
    call MPI_Comm_Size(MPI_Comm_World,MPI_Size,MPIerror)
!    call MPI_Dims_Create(MPI_Size, 2, dims, MPIerror)
    call ms_dbcsr_setup(MPI_Comm_World)
  end if
  m_storage ='pdcsr' 
  m_operation ='lap'
  BlockSize_c = BlockSize
  if (ionode) print'(a,i5)','BlockSize_c = ', BlockSize_c
#else
  dims(1) = 1
  dims(2) = 1
  m_storage = 'sdcsr'
  m_operation = 'ref'
  BlockSize_c = h_dim
  if(ionode) then
    write(6,*) 'omm_min: Block matrices are not supported in serial'
  endif
  call die()
#endif
    
  if(first_call) then
    if (ionode) print'(a)','OMM with libOMM and block matrices'
    
    ! Set the block sizes
    nblocks_r = ceiling(real(h_dim,dp)/dims(1))
    nblocks_c = ceiling(real(h_dim,dp)/dims(2))
    nblocks_r1 = ceiling(real(N_occ,dp)/dims(1))
    nblocks_c1 = ceiling(real(N_occ,dp)/dims(2))
    allocate(row_blk_sizes(1:nblocks_r))
    allocate(col_blk_sizes(1:nblocks_c))
    row_blk_sizes(:) = dims(1)
    col_blk_sizes(:) = dims(2)
    allocate(row_blk_sizes1(1:nblocks_r1))
    allocate(col_blk_sizes1(1:nblocks_c1))
    row_blk_sizes1(:) = dims(1)
    col_blk_sizes1(:) = dims(2)
  end if
      
  allocate(ind_o(1:nhmax))
  allocate(ind_u(1:nbasis))   
  do io = 1, nbasis 
    ind_u(io) = listhptr(io) + 1
    ind_o(listhptr(io) + 1) = listhptr(io) + 1
    do j = 2, numh(io)
      ind = listhptr(io) + j
      ind_o(ind) = ind
      jo = listh(ind)
      if(jo .lt. listh(ind_o(ind - 1))) then
        do kk = 1, j - 1
          ind1 = listhptr(io) + j - kk
          jo1 = listh(ind_o(ind1))
          if(jo1 .gt. jo) then
            ind_o(ind1 + 1) = ind_o(ind1)
            ind_o(ind1) = ind
          end if
        end do
      end if
    end do
  end do
   
  if(first_call) then
    if (.not. H%is_initialized) call m_allocate(H,row_blk_sizes,col_blk_sizes,m_storage)
    if (.not. S%is_initialized) call m_allocate(S,row_blk_sizes,col_blk_sizes,m_storage)
    if (.not. D_min%is_initialized) call m_allocate(D_min,row_blk_sizes,col_blk_sizes,m_storage)
    if (.not. C_min%is_initialized) call m_allocate(C_min,row_blk_sizes1,col_blk_sizes,m_storage)
  end if

  if(ionode) print'(a)',    'Allocation'
  
  ib_c = 1
  ib_r = 1 
  blk_co = col_blk_sizes(ib_c)
  blk_ro = row_blk_sizes(ib_r)
  blk_f = col_blk_sizes(ib_c)
  allocate(block_data(1:blk_ro,1:blk_co))
  allocate(block_data_s(1:blk_ro,1:blk_co))
  allocate(found_c(1:blk_f))
  num_c = 0
  do ib_c = 1, nblocks_c
    ind_c2 = col_blk_sizes(ib_c)
    if(ib_c == nblocks_c) ind_c2 = h_dim - num_c
    do ib_r = 1, nblocks_r  
      blk_c = col_blk_sizes(ib_c)
      blk_r = row_blk_sizes(ib_r)
      if(blk_f .ne. blk_c) then
        deallocate(found_c) 
        allocate(found_c(1:blk_c))
        blk_f = blk_c
      end if
      found_c(:) = .false.
      found = .false.
      do ind_c = 1, ind_c2
        i_node = int((num_c + ind_c - 1) / BlockSize_c)
        if(Node == i_node) then
          io = mod((num_c + ind_c - 1), BlockSize_c) + 1
          ind = ind_u(io)
          if(ind .le. (listhptr(io) + numh(io))) then
            jo = listh(ind_o(ind))
            call get_index(jo, row_blk_sizes, nblocks_r, ib, ind_r)
          else
            ib = 0
          end if
          if(ib == ib_r) then
            found_c(ind_c) = .true.
            found = .true.
          end if
        end if
      end do
#ifdef MPI
      do ind_c = 1, ind_c2
        i_node = int((num_c + ind_c - 1) / BlockSize_c)
        call MPI_Bcast(found_c(ind_c),1,MPI_Logical,i_node,MPI_Comm_World,MPIerror)
        if(found_c(ind_c)) then
          found = .true.
        end if
      end do
#endif
      if(found) then
        if((blk_ro .ne. blk_r) .or. (blk_co .ne. blk_c)) then
          deallocate(block_data)
          deallocate(block_data_s)
          allocate(block_data(1:blk_r,1:blk_c))
          allocate(block_data_s(1:blk_r,1:blk_c))
          blk_co = blk_c
          blk_ro = blk_r
        end if
        block_data(:,:) = 0.0_dp
        block_data_s(:,:) = 0.0_dp
        do ind_c = 1, ind_c2
          i_node = int((num_c + ind_c - 1) / BlockSize_c)
          if(Node == i_node) then
            io = mod((num_c + ind_c - 1), BlockSize_c) + 1
            ind = ind_u(io)
            jo = listh(ind_o(ind))
            call get_index(jo, row_blk_sizes, nblocks_r, ib, ind_r) 
            do while(ib == ib_r) 
              he = h_sparse(ind_o(ind), 1)
              se = s_sparse(ind_o(ind))
         !   print'(a,i5,a,i5,a,i5,a,i5,a,i5,a,f15.10)','H i_node = ',i_node,' ib_c=',ib_c,&
         !      ' ib_r=',ib_r,' ind_c=',ind_c,' ind_r=',ind_r,' h=',he
              block_data(ind_r,ind_c) = he
              block_data_s(ind_r,ind_c) = se
              ind_u(io) = ind + 1
              ind = ind_u(io)
              if(ind .le. (listhptr(io) + numh(io))) then
                jo = listh(ind_o(ind))
                call get_index(jo, row_blk_sizes, nblocks_r, ib, ind_r)
              else
                ib = 0
              end if 
            end do
          end if
        end do
#ifdef MPI      
        do ind_c = 1, ind_c2
          i_node = int((num_c + ind_c - 1) / BlockSize_c) 
          call MPI_Bcast(block_data(:,ind_c),blk_r,MPI_Double_Precision,i_node,MPI_Comm_World,MPIerror)
          call MPI_Bcast(block_data_s(:,ind_c),blk_r,MPI_Double_Precision,i_node,MPI_Comm_World,MPIerror)      
        end do
#endif      
        call m_set_element(H,ib_r,ib_c,block_data,0.0_dp)     
        call m_set_element(S,ib_r,ib_c,block_data_s,0.0_dp)
      end if
    end do
    num_c = num_c + ind_c2
  end do     
  deallocate(block_data_s) 
  if(ionode) print'(a)',    'H and S set'
  
  if(first_call) then
!    seed = omm_rand_seed()
!    do i = 1, nblocks_c
!      do j = 1, nblocks_r1
!        blk_c = col_blk_sizes(i)
!        blk_r = row_blk_sizes1(j)
!        if((blk_ro .ne. blk_r) .or. (blk_co .ne. blk_c)) then
!          deallocate(block_data)
!          allocate(block_data(1:blk_r,1:blk_c))
!          blk_co = blk_c
!          blk_ro = blk_r
!        end if        
!        block_data(:,:) = 0.0_dp
!        do io = 1, dims(2) 
!          do jo = 1, dims(1)
!            el = 0.0_dp
!            ind = (i-1) * dims(2) + io
!            ind1 = (j-1) * dims(1) + jo
!            if((ind .le. h_dim) .and. (ind1 .le. N_occ)) then
!              do k = 1, 2
!                call omm_bsd_lcg(seed, rn(k))
!              end do
!              el = sign(0.5_dp * rn(1), rn(2) - 0.5_dp)
!            end if
!            block_data(jo, io) = el 
!          end do
!        end do
!        call m_set_element(C_min,j,i,block_data(:,:),0.0_dp)
!      end do
!    end do
!    call m_scale(C_min,1.0d-2/sqrt(real(h_dim,dp)),m_operation)

    rmax = 0.0_dp
    do i = -1,1
      do j = -1,1
        do k = -1,1
          rr(1) = i*ucell(1,1) + j*ucell(1,2) + k*ucell(1,3)
          rr(2) = i*ucell(2,1) + j*ucell(2,2) + k*ucell(2,3)
          rr(3) = i*ucell(3,1) + j*ucell(3,2) + k*ucell(3,3)
          rrmod = sqrt( rr(1)**2 + rr(2)**2 + rr(3)**2 )
          if (rrmod .gt. rmax) rmax = rrmod
        enddo
      enddo
    enddo
    if(ionode) print'(a,f13.7)',    'rmax = ', rmax
    
    call mneighb(ucell,rcoor,na_u,xa,0,0,nna)
    if(ionode) print'(a,f13.7)',    'rcoor = ', rcoor
    
    seed = omm_rand_seed()

    ia_off = 0
    do ia =  1, na_u !! loop_ia
      if(2.0 * rcoor .lt. rmax) then
        call mneighb(ucell,rcoor,na_u,xa,ia,0,nna) ! 
      else
        nna = na_u
        do j = 1, na_u
          jan(j) = j
        end do
      end if
      allocate(neib(1:nna))
      neib(:) = 0
      index = 0
      do i = 1, nna
        found = .false.
        do j = 1, index
          if(neib(j) == jan(i)) found = .true.
        end do
        if(.not. found) then
          index = index + 1
          neib(index) = jan(j)
        end if
      end do
      nna = index
      do i = 2, nna
        neib0 = neib(i)
        j = i - 1
        do while (neib0 .lt. neib(j)) 
            neib(j + 1) = neib(j)
            neib(j) = neib0
            j = j - 1        
        end do
      end do
           
!      if(ionode) print'(a,i5,a,i5)',  'ia = ', ia,  '    nna = ', nna 
!      if(ionode) then 
!        do i = 1, nna
!          print'(a,i5,a,i5,a,i5)',  'ia = ', ia,  '    i_neib = ', i, ' neib = ', neib(i)
!        end do
!      end if 

      call get_number_of_lwfs_on_atom(ia, indexi)  ! gets indexi, nelectr       
!      if(ionode) print'(a,i5,a,i5)',  'ia = ', ia,  '    indexi = ', indexi
       
      ind_c = 0
      ib_r1 = 0
      ib_r2 = -1
      if(indexi .gt. 0) then
        call get_index(ia_off+1, row_blk_sizes1, nblocks_r1, ib_r1, ind_a1)
        call get_index(ia_off+indexi, row_blk_sizes1, nblocks_r1, ib_r2, ind_a2)
        ia_off = ia_off + indexi
      end if

 !    if(ionode) print'(a,i5,a,i5,a,i5,a,i5,a,i5)','ia =', ia,' ib_r1= ',ib_r1,' ib_r2= ',ib_r2 ,&
 !          ' ind_a1= ',ind_a1,' ind_a2= ',ind_a2 
      
      do ib_r = ib_r1, ib_r2
        ind_r1 = 1
        blk_r = row_blk_sizes1(ib_r)
        ind_r2 = blk_r
        if(ib_r == ib_r1) ind_r1 = ind_a1
        if(ib_r == ib_r2) ind_r2 = ind_a2
        ib_co = 0
        ind_co = blk_co
        do  j = 1, nna   !!loop_neighbor_atoms
          ja = neib(j)
          norb = lasto(ja) - lasto(ja-1)
          call get_index(lasto(ja-1)+1, col_blk_sizes, nblocks_c, ib_c, ind_c)
          blk_c = col_blk_sizes(ib_c)

 !         if(ionode) print'(a,i5,a,i5,a,i5,a,i5)','ia =', ia,' ib_r= ',ib_r,' ib_c= ',ib_c ,&
 !          ' ind_c= ',ind_c
          if(ib_c .gt. ib_co) then
            if(ind_co .lt. blk_co) then
              call m_set_element(C_min,ib_r,ib_co,block_data(:,:),0.0_dp)
!             do io = 1, blk_r
!                do jo = 1, blk_co
!                  if(ionode) print'(a,i5,a,i5,a,i5,a,i5,a,f15.10)','Block_set0 ib_r= ',ib_r,' ib_c= ',ib_c ,&
!                     ' ind_r= ',io, ' ind_c=', jo, ' bl = ', block_data(io,jo)
!                 end do
!              end do
            end if
            if((blk_ro .ne. blk_r) .or. (blk_co .ne. blk_c)) then
              deallocate(block_data)
              allocate(block_data(1:blk_r,1:blk_c))
              blk_co = blk_c
              blk_ro = blk_r
            end if            
            block_data(:,:) = 0.0_dp
            if(ind_r1 .gt. 1) then
              call get_block(C_min,ib_r,ib_c,blk_r,blk_c,block_data)
!              do io = 1, blk_r
!                 do jo = 1, blk_c
!                   if(ionode) print'(a,i5,a,i5,a,i5,a,i5,a,f15.10)','Block_read ib_r= ',ib_r,' ib_c= ',ib_c ,&
!                  ' ind_r= ',io, ' ind_c=', jo, ' bl = ', block_data(io,jo)
!                 end do
!              end do 
            end if
          end if 
          do iorb = 1, norb  !!loop_orbs_on_neighbor            
            do ind_r = ind_r1, ind_r2
              do k = 1, 2
                call omm_bsd_lcg(seed, rn(k))
              end do
              cgval = sign(0.5_dp * rn(1), rn(2) - 0.5_dp)
              block_data(ind_r,ind_c) = cgval
            end do
            if(ind_c == blk_c) then
              call m_set_element(C_min,ib_r,ib_c,block_data(:,:),0.0_dp)
 !             do io = 1, blk_r
 !                do jo = 1, blk_c
 !                  if(ionode) print'(a,i5,a,i5,a,i5,a,i5,a,f15.10)','Block_set ib_r= ',ib_r,' ib_c= ',ib_c ,&
 !                 ' ind_r= ',io, ' ind_c=', jo, ' bl = ', block_data(io,jo)
 !                end do
 !             end do
              ind_c = ind_c - blk_c
              ib_c = ib_c + 1
              blk_c = col_blk_sizes(ib_c)
              if((ib_c .le. nblocks_c) .and. (iorb .lt. norb)) then
                if((blk_ro .ne. blk_r) .or. (blk_co .ne. blk_c)) then
                  deallocate(block_data)
                  allocate(block_data(1:blk_r,1:blk_c))
                  blk_co = blk_c
                  blk_ro = blk_r
                end if
                block_data(:,:) = 0.0_dp
                if(ind_r1 .gt. 1) then
                  call get_block(C_min,ib_r,ib_c,blk_r,blk_c,block_data)
 !             do io = 1, blk_r
 !                do jo = 1, blk_c
 !                  if(ionode) print'(a,i5,a,i5,a,i5,a,i5,a,f15.10)','Block_read2 ib_r= ',ib_r,' ib_c= ',ib_c ,&
 !                 ' ind_r= ',io, ' ind_c=', jo, ' bl = ', block_data(io,jo)
 !                end do
 !             end do
                end if
              end if
            end if
            ind_co  = ind_c
            ind_c = ind_c + 1
            ib_co = ib_c
          end do
        end do
        if((ib_co .gt. 0) .and. (ind_co .lt. blk_co)) then
          call m_set_element(C_min,ib_r,ib_co,block_data(:,:),0.0_dp)
!              do io = 1, blk_r
!                 do jo = 1, blk_co
!                   if(ionode) print'(a,i5,a,i5,a,i5,a,i5,a,f15.10)','Block_set2 ib_r= ',ib_r,' ib_c= ',ib_c ,&
!                  ' ind_r= ',io, ' ind_c=', jo, ' bl = ', block_data(io,jo)
!                 end do
!              end do
        end if
      end do
      deallocate(neib)
    end do
    call m_scale(C_min,1.0d-2/sqrt(real(h_dim,dp)),m_operation)
  end if
  
  init_C = .true.
  if(.not. calcE) then
    call omm(h_dim,N_occ,H,S,.true.,e_min,D_min,.false.,0.0_dp,&
      C_min,init_C,T,0.0_dp,0,nspin,1,-1.0_dp,.true.,.false.,&
      m_storage,m_operation,row_blk_sizes,col_blk_sizes,&
      row_blk_sizes1,col_blk_sizes1)
    if(ionode) print'(a, f13.7)','e_min = ', e_min
  else
    call omm(h_dim,N_occ,H,S,.false.,e_min,D_min,.true.,0.0_dp,&
      C_min,init_C,T,0.0_dp,0,nspin,1,-1.0_dp,.true.,.true.,&
      m_storage,m_operation,row_blk_sizes,col_blk_sizes,& 
      row_blk_sizes1,col_blk_sizes1)
  end if
          
  do io = 1, nbasis
    ind_u(io) = listhptr(io) + 1
  end do
  num_c = 0
  do ib_c = 1, nblocks_c
    ind_c2 = col_blk_sizes(ib_c)
    if(ib_c == nblocks_c) ind_c2 = h_dim - num_c
    do ib_r = 1, nblocks_r
      blk_c = col_blk_sizes(ib_c)
      blk_r = row_blk_sizes(ib_r)
      if(blk_f .ne. blk_c) then
        deallocate(found_c)
        allocate(found_c(1:blk_c))
        blk_f = blk_c
      end if
      if((blk_ro .ne. blk_r) .or. (blk_co .ne. blk_c)) then
        deallocate(block_data)
        allocate(block_data(1:blk_r,1:blk_c))
        blk_co = blk_c
        blk_ro = blk_r
      end if
      block_data(:,:) = 0.0_dp
      found_c(:) = .false.
      found = .false.
      do ind_c = 1, ind_c2
        i_node = int((num_c + ind_c - 1) / BlockSize_c)
        if(Node == i_node) then
          io = mod((num_c + ind_c - 1), BlockSize_c) + 1
          ind = ind_u(io)
          if(ind .le. (listhptr(io) + numh(io))) then
            jo = listh(ind_o(ind))
            call get_index(jo, row_blk_sizes, nblocks_r, ib, ind_r)
          else
            ib = 0
          end if
          if(ib == ib_r) then
            found_c(ind_c) = .true.
            found = .true.
          end if
        end if
      end do
#ifdef MPI
      do ind_c = 1, ind_c2
        i_node = int((num_c + ind_c - 1) / BlockSize_c)
        call MPI_Bcast(found_c(ind_c),1,MPI_Logical,i_node,MPI_Comm_World,MPIerror)
        if(found_c(ind_c)) then
          found = .true.
        end if
      end do
#endif
      
      if(found) then
        if((blk_ro .ne. blk_r) .or. (blk_co .ne. blk_c)) then
          deallocate(block_data)
          allocate(block_data(1:blk_r,1:blk_c))
          blk_co = blk_c
          blk_ro = blk_r
        end if
        block_data(:,:) = 0.0_dp
        call get_block(D_min,ib_r,ib_c,blk_r,blk_c,block_data)
      
        do ind_c = 1, ind_c2
          i_node = int((num_c + ind_c - 1) / BlockSize_c)
          if(Node == i_node) then
            io = mod((num_c + ind_c - 1), BlockSize_c) + 1
            ind = ind_u(io)
            jo = listh(ind_o(ind))
            call get_index(jo, row_blk_sizes, nblocks_r, ib, ind_r)
            do while(ib == ib_r)
              d_sparse(ind_o(ind), 1) = 2.0_dp * block_data(ind_r, ind_c)
!            print'(a,i5,a,i5,a,i5,a,i5,a,f15.10)','Output  ib_c=',ib_c,&
!               ' ib_r=',ib_r,' ind_c=',ind_c,' ind_r=',ind_r,' h=',block_data(ind_r,ind_c)
              if(nspin == 2) d_sparse(ind_o(ind), 2) = d_sparse(ind_o(ind), 1)
              ind_u(io) = ind + 1
              ind = ind_u(io)
              if(ind .le. (listhptr(io) + numh(io))) then
                jo = listh(ind_o(ind))
                call get_index(jo, row_blk_sizes, nblocks_r, ib, ind_r)
              else
                ib = 0
              end if
            end do
          end if
        end do
      end if
    end do
    num_c = num_c + ind_c2
  end do

  if(first_call) first_call = .false.
   
!  call m_deallocate(D_min)
!  call m_deallocate(S)
!  call m_deallocate(H)
!  call m_deallocate(C_min)
!  call ms_dbcsr_finalize()

  deallocate(ind_o)
  deallocate(ind_u)
  deallocate(found_c)
  deallocate(block_data)

end subroutine omm_min_block

subroutine get_index(num, blk_sizes, num_blks, ib, ind)
  implicit none

  integer, intent(in)     :: num
  integer, intent(inout)  :: blk_sizes(:) 
  integer, intent(in)     :: num_blks
  integer, intent(out)    :: ib, ind

  ib = 0
  ind = num

  do while((ib .lt. num_blks) .and. (ind .gt. 0))
    ib = ib + 1
    ind = ind - blk_sizes(ib)
  end do
  if(ind .le. 0) then
    ind = ind + blk_sizes(ib)
  else    
    if(ionode) then
      write(6,*) 'omm_min: Wrong matrix index'
    endif
    call die()
  end if
end subroutine get_index


subroutine get_block(mat, ib_r, ib_c, blk_r, blk_c, block_data)
  implicit none

  type(matrix), intent(inout)   :: mat
  integer, intent(in)        :: ib_r, ib_c
  integer, intent(in)        :: blk_r, blk_c
  real(dp), intent(inout)    :: block_data(:,:)

  logical :: found
  integer :: i_node, ind_r, ind_c
  real(dp),pointer :: myblock(:,:)

#ifdef MPI
  integer :: MPIerror, i_mpi
#endif

  call m_get_element(mat,ib_r,ib_c,myblock,found)
  i_node = 0
  if(found) then
    i_node = Node
    do ind_c = 1, blk_c
      do ind_r = 1, blk_r
        block_data(ind_r,ind_c) = myblock(ind_r,ind_c)
      end do
    end do
  end if

#ifdef MPI
  call MPI_Allreduce(i_node,i_mpi,1,MPI_Integer,MPI_Sum,MPI_Comm_World,MPIerror)
  i_node = i_mpi
  do ind_c = 1, blk_c
    call MPI_Bcast(block_data(:,ind_c),blk_r,MPI_Double_Precision,i_node,MPI_Comm_World,MPIerror)
  end do
#endif
end subroutine get_block

subroutine get_number_of_lwfs_on_atom(ia, indexi)
  implicit none
  
  integer, intent(in)  :: ia
  integer, intent(out) :: indexi
  
  real(dp) :: tiny
  integer  :: nelectr
  logical, save :: secondodd = .false.
 
  tiny = 1.d-10
  nelectr = qa(ia) + tiny 
  if (abs(nelectr - qa(ia) + tiny) .gt. 1e-3) then
    if(ionode) then
      write(6,*) 'omm_min: Wrong atomic charge for atom ',ia
      write(6,*) '      qa = ',qa(ia),' must be an integer'
    endif
    call die()
  endif
  if ( (nelectr / 2) * 2 .ne. nelectr) then
    if (secondodd) then
      indexi = ( ( nelectr - 1 ) / 2 )
      secondodd = .false.
    else
      indexi = ( ( nelectr + 1 ) / 2 )
      secondodd = .true.
    endif
  else
    indexi = ( ( nelectr ) / 2 )
  endif

end subroutine get_number_of_lwfs_on_atom

end module m_lib_omm
