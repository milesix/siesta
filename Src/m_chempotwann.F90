! ----
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!
!> \brief General purpose of the m_chempotwann module:
!! Here we compute the matrix elements of the Hamiltonian in real space
!! between atomic orbitals in the basis (with the same sparsity as usual)
!! when a "chemical potential" is applied to a given Wannier
!!
!! The matrix elements are stored in an array called
!! H_chempotwann
!! defined in sparse_matrices.F90 under the name of H_chempotwann_2D
!! Then, these matrix elements are initialized in state_init.F
!! The pointer H_chempotwann is forced to point over H_chempotwann_2D
!! in setup_hamiltonian.F and siesta_analysis.F
!!
!! The analytical expressions can be found in
!!
!! <https://personales.unican.es/junqueraj/JavierJunquera_files/Notes/Wannier/wannier_in_nao.pdf>
!!

  module m_chempotwann

  use precision,        only: dp               ! Double precision
  use atomlist,         only: no_u             ! Number of orbitals in the unit cell
  use atomlist,         only: no_s             ! Number of atoms in the supercell
  use atomlist,         only: no_l             ! Number of orbitals (local)
  use sparse_matrices,  only: maxnh            ! Number of non-zero elements in 
                                               !   the sparse matrix 
                                               !   (local to MPI Node)
  use sparse_matrices,  only: numh             ! Number of non-zero elements 
                                               !   per row in the CSR matrix 
                                               !   (local to MPI Node)
  use sparse_matrices,  only: listhptr         ! Pointer to the start of 
                                               !   ach row (-1) of the
                                               !   hamiltonian matrix
  use sparse_matrices,  only: listh            ! Nonzero hamiltonian-matrix 
                                               !   element column indices for 
                                               !   each matrix row
  use sparse_matrices,  only: H                ! Hamiltonian-matrix in real space
  use m_spin,           only: spin             ! Spin configuration for SIESTA


#ifdef MPI
    use mpi_siesta
#endif

  implicit none

  private

  public :: chempotwann
  public :: add_Hamiltonian_chempotwann

  CONTAINS

!> \brief General purpose of the subroutine chempotwann
!!
!!
  subroutine chempotwann ( ispin, H_chempotwann )

  use parallel,         only: Node             ! Local processor number
  use parallel,         only: Nodes            ! Total number of processors in a
                                               !  parallel run
  use parallel,         only: BlockSize        ! Blocking factor used to 
                                               !   divide up the arrays over 
                                               !   the processes for the
                                               !   Scalapack routines.
  use parallelsubs,     only: WhichNodeOrb     ! Which node handles the 
                                               !   the information of a given
                                               !   orbital in the unit cell
  use parallelsubs,     only: GlobalToLocalOrb ! Transformation from a global
                                               !   index of a given orbital 
                                               !   in the unit cell to the 
                                               !   local index within a node
  use parallelsubs,     only: LocalToGlobalOrb ! Subroutine that transforms 
                                               !   a local orbital index in 
                                               !   a given node to the global 
                                               !   atomic index
  use m_wannier_in_nao, only: coeffs_wan_nao   ! Coefficients of the Wannier 
                                               !   functions in a basis of NAO
                                               !   First index: Wannier function
                                               !   Second index: NAO in the 
                                               !      supercell
  use w90_in_siesta_types, only: chempotwann_val
                                               ! Chemical potential
                                               !   applied to shift the energy 
                                               !   of a matrix elements in real
                                               !   space associated with a given
                                               !   Wannier function
  use w90_in_siesta_types,   only: num_proj_local
                                               ! Number of projections that will
                                               !  be  handled in the local node
  use w90_parameters,   only: num_proj         ! Number of projections
  use sys,              only: die              ! Termination routine
  use alloc,            only: re_alloc         ! Allocatation routines
  use alloc,            only: de_alloc         ! Deallocatation routines

#ifdef MPI
  use parallelsubs,     only: set_blocksizedefault
                                               ! Finds a sensible default value
                                               !   for the BlockSize default.
#endif

  integer,  intent(in)    :: ispin             ! Counter of spin components
  real(dp), intent(inout) :: H_chempotwann(maxnh,spin%H)   
                                               ! Extra term in the Hamiltonian 
                                               !   that accounts for a 
                                               !   rigid shift of the 
                                               !   bands associated
                                               !   with a given Wannier

  external :: timer

! Internal variables ................................................

  integer :: iuo       ! Counter on orbitals in the unit cell
  integer :: j         ! Counter of neighbours of a given orbital
  integer :: jneig     ! Index of the neighbour orbital (in supercell notation)
  integer :: ind       ! Index of the neighbour orbital in listh
  integer :: iproj_local  ! Local index within one node of the Wannier function
  integer :: iproj_global ! Global index of the Wannier function

  real(dp), pointer :: H_chempotwann_full(:,:) => null() 
                       ! New matrix elements coming from the chemical potential
                       !   applied to the Wanniers
                       !   Every nodes computes the full line in the 
                       !   sparse matrix, but summing only over the 
                       !   Wanniers known by the Node
  integer :: jo        ! Counter for the loops on neighbour orbitals

! Variable required to globalize the neighbour list
  integer :: io_local  ! Counter for the local atomic orbitals in a given node
  integer :: io_global ! Global index of an atomic orbital
  integer :: BNode     ! Node that contains the information of a given
                       !    orbital in the unit cell
  integer :: maxnhg    ! Compute the maximum number of 
                       !    interacting atomic orbitals neighbours
  integer,  pointer :: numhg(:) => null()
                       ! Temporal array used to broadcast the array numh,
                       ! with number of non-zero elements 
                       ! each orbital connects to
  integer,  pointer :: listhptrg(:) => null()
                       ! Temporal array used to broadcast the array listhptr
                       ! index pointer to listh such that listh(listhptr(1) + 1)
                       ! is the first non-zero element of orbital 1.
                       ! listh(listhptr(io) + 1) is thus the 
                       ! first non-zero element
                       ! of orbital 'io' while listh(listhptr(io) + numh(io)) 
                       ! is the last non-zero element of orbital 'io'.
  integer,  pointer :: listhg(:) => null()
                       ! Temporal array used to broadcast the array listh
                       !      the column indices for the non-zero elements

#ifdef MPI
  integer :: blocksize_save
  integer :: blocksizeincbands_wannier    
                       !  BlockSize when the number of projectors are
                       !    distributed over the nodes
  integer :: MPIerror  ! MPI code error
  integer :: maxnumh   ! Maximum value in numh
  real(dp), pointer :: H_loc(:) => null() 
                       ! Sum of all the contributions to the Hamiltonian 
                       !   from all the nodes (i.e. from all the projections)
                       !   It takes the full rows of the Hamiltonian for
                       !   the different nodes, adds them, and split
                       !   according to the parallel sparsity
#endif 

!  Start time counter
   call timer( 'chempotwann', 1 )

#ifdef MPI
!  Calculate the default value for BlockSize, when we distribute the
!  number of projectors over the nodes
   call set_blocksizedefault( Nodes, num_proj,                                &
 &                            blocksizeincbands_wannier )

!!  For debugging
!   blocksize_save = BlockSize
!   BlockSize = blocksizeincbands_wannier
!   do iproj_local = 1, num_proj_local
!     call LocalToGlobalOrb(iproj_local, Node, Nodes, iproj_global)
!     write(6,'(a,i2,a,f12.5,2i5)')      &
! &     'chempotwann: Node, Nodes, Chemical potential for Wannier ',           &
! &     iproj_global,                                                          &
! &     ' = ', chempotwann_val(iproj_global), Node, Nodes
!   enddo 
!   BlockSize = blocksize_save
!!   call MPI_barrier(MPI_Comm_world,MPIerror)
!!   call die()
!!  End debugging
#endif

!  Initialize the matrix elements of the Hamiltonian coming from the
!  chemical potential applied to the Wannier functions
   H_chempotwann(:,ispin) = 0.0_dp

!  H_chempotwann follows the same sparsity pattern as the Hamiltonian in 
!  real space.
!  The first thing to do is to globalise the neighbour list arrays, so
!  all the nodes will know the list of neighbours of all the orbitals in the
!  unit cell
!  Copy of the procedure implemented in subroutine diagkp
   nullify(numhg, listhptrg)

!  Allocate local memory for global list arrays
   call re_alloc( numhg, 1, no_u, name='numhg',                               &
 &                routine= 'chempotwann' )
   call re_alloc( listhptrg, 1, no_u, name='listhptrg',                       &
 &                routine= 'chempotwann' )

!  Globalise numh
!  Loop over all the orbitals in the unit cell
   do io_global = 1, no_u
!    Localize which node handles the information related with the orbital io
     call WhichNodeOrb(io_global,Nodes,BNode)
     if ( Node .eq. BNode ) then
!      Identify the local index for the orbital in the unit cell in the 
!      node that handles its information
       call GlobalToLocalOrb(io_global,Node,Nodes,io_local)
!      Assign the value of the number of neighbours of that particular orbital
       numhg(io_global) = numh(io_local)
     endif
!    Transfer the information from the node that contains the information
!    on the orbital io_global to all the other nodes
#ifdef MPI
     call MPI_Bcast( numhg(io_global),1,MPI_integer,BNode,                    &
 &                   MPI_Comm_World,MPIerror )
#endif
   enddo

!  Build global listhptr
   listhptrg(1) = 0
   do io_global = 2, no_u
     listhptrg(io_global) = listhptrg(io_global-1) + numhg(io_global-1)
   enddo

!  Globalise listh
!  Compute the maximum number of interacting atomic orbitals neighbours 
!  considering all the orbitals in the unit cell
   maxnhg = listhptrg(no_u) + numhg(no_u)
   nullify(listhg)
   call re_alloc( listhg, 1, maxnhg, name='listhg',                           &
 &                routine= 'chempotwann' )
   do io_global = 1, no_u
     call WhichNodeOrb(io_global,Nodes,BNode)
     if (Node.eq.BNode) then
       call GlobalToLocalOrb(io_global,Node,Nodes,io_local)
       do jo = 1, numhg(io_global)
         listhg(listhptrg(io_global)+1:listhptrg(io_global)+numhg(io_global))=&
 &         listh(listhptr(io_local)+1:listhptr(io_local)+numh(io_local))
       enddo
     endif
#ifdef MPI
     call MPI_Bcast( listhg(listhptrg(io_global)+1),numhg(io_global),         &
 &                   MPI_integer,BNode,MPI_Comm_World,MPIerror )
#endif 
   enddo

#ifdef MPI
!  Find maximum value in numh and create local storage
   maxnumh = 0
   do io_local = 1, no_l
     maxnumh = max(maxnumh,numh(io_local))
   enddo
   nullify(H_loc)
   call re_alloc( H_loc, 1, maxnumh, name='H_loc',                            &
 &                routine= 'chempotwann' )
#endif


!  Create new distribution of the Hamiltonian matrix containing the 
!  penalty for a given Wannier function
   nullify( H_chempotwann_full )
   call re_alloc( H_chempotwann_full, 1, maxnhg, 1, spin%H,                   &
 &                name='H_chempotwann_full', routine= 'chempotwann' )
   H_chempotwann_full = 0.0_dp

!!  For debugging
!   do io_global = 1, no_l
!     write(6,'(a,4i5)') &
! &     'Node, Nodes, io_global, numh = ', &
! &      Node, Nodes, io_global, numh(io_global) 
!   enddo
!   do io_global = 1, no_u
!     write(6,'(a,5i5)')                                                      & 
! &     'Node, Nodes, io, numhg(io), listhptrg(io) = ',                       &
! &      Node, Nodes, io_global, numhg(io_global), listhptrg(io_global)
!     do jo = 1, numhg(io_global)
!       ind = listhptrg(io_global) + jo
!       write(6,'(a,5i5)')                                                    & 
! &     'Node, Nodes, io, jo, listhg(io) = ',                                 &
! &      Node, Nodes, io_global, jo, listhg(ind)
!     enddo 
!   enddo 
!   do iproj_local = 1, num_proj_local
!     do io_global = 1, no_s
!       write(6,'(a,4i5,2f12.5)') &
! &       'Node, Nodes, iproj, iorb, coeffs_wan_nao(iproj_local,iorb) = ',   &
! &        Node, Nodes, iproj_local, io_global, coeffs_wan_nao(iproj_local,io_global) 
!     enddo
!   enddo
!#ifdef MPI
!   call MPI_barrier(MPI_Comm_world,MPIerror)
!#endif
!   call die()
!!  End debugging

!  Compute the Hamiltonian elements with the penalty for the Wannier functions
!  First loop on all the local orbitals in the unit cell
   do iuo = 1, no_u
 
!    Then, loop on all the neighbour orbitals
     do j = 1, numhg(iuo)

!      Identify the neighbour orbital (in supercell notation)
       ind   = listhptrg(iuo) + j
       jneig = listhg(ind)

!!      For debugging
!       write(6,'(a,7i5)')                                       &
! &       'Node, Nodes, iuo, j, ind, jneig, num_proj_local = ',  &
! &        Node, Nodes, iuo, j, ind, jneig, num_proj_local
!!      End debugging

!      Compute the extra term in the Hamiltonian
!      One node knows the Hamiltonian with the penalty of all the Wannier
!      projections handled locally in that node,
!      but it knows all the terms in the sparse matrix
#ifdef MPI
       blocksize_save = BlockSize
       BlockSize = blocksizeincbands_wannier
#endif
!      Loop on the Wanniers stored in the local node
       do iproj_local = 1, num_proj_local
!        Identify the global index of the Wannier
#ifdef MPI
         call LocalToGlobalOrb(iproj_local, Node, Nodes, iproj_global)
#else 
         iproj_global = iproj_local
#endif
         H_chempotwann_full(ind,ispin) = H_chempotwann_full(ind,ispin)     +   &
 &                                   coeffs_wan_nao(iproj_local,iuo)       *   &
 &                                   chempotwann_val(iproj_global)         *   &
 &                                   coeffs_wan_nao(iproj_local,jneig)  
       enddo 
#ifdef MPI
       BlockSize = blocksize_save
#endif

!!      For debugging
!       if( H_chempotwann_full(ind,ispin) .gt. 1.d-6)                    &
! &     write(6,'(a,7i7,f12.5)')                                         &
! & 'Node, Nodes, iuo, j, ind, jneig, ispin, H_chempotwann_full = ',     &
! &  Node, Nodes, iuo, j, ind, jneig, ispin,                             &
! &  H_chempotwann_new(ind,ispin) 
!!      End debugging

     enddo ! End loop on neighbours
   enddo   ! End loop on atomic orbitals

!  First loop on all the local orbitals in the unit cell
#ifdef MPI
   do io_global = 1, no_u
     call WhichNodeOrb( io_global,Nodes,BNode )
     call GlobalToLocalOrb(io_global,BNode,Nodes,io_local)
     call MPI_Reduce( H_chempotwann_full(listhptrg(io_global)+1,ispin),       &
 &                    H_loc(1),numhg(io_global),MPI_double_precision,         &
 &                    MPI_sum,BNode,MPI_Comm_World,MPIerror )

!    Then, loop on all the neighbour orbitals
     if ( Node .eq. BNode ) then
       do j = 1,numh(io_local)
         H_chempotwann(listhptr(io_local)+j,ispin) = H_loc(j)
       enddo ! End loop on neighbours
     endif

   enddo   ! End loop on atomic orbitals
#else
   H_chempotwann = H_chempotwann_full
#endif

!!  For debugging
!   do io_local = 1, no_l
!     do j = 1, numh(io_local)
!       ind = listhptr(io_local) + j
!       write(6,'(a,5i5,f12.5)')                                               &
! &       'Node, Nodes, io_local, j, ind, H_chempotwann = ',                   &
! &        Node, Nodes, io_local, j, ind, H_chempotwann(ind,ispin)
!     enddo 
!   enddo 
!
!#ifdef MPI
!   call MPI_barrier(MPI_Comm_world,MPIerror)
!#endif
!   call die()
!!  End debugging

#ifdef MPI
!  Free local memory
   call de_alloc( H_loc,     name='H_loc',     routine= 'chempotwann' )
   call de_alloc( listhg,    name='listhg',    routine= 'chempotwann' )
   call de_alloc( listhptrg, name='listhptrg', routine= 'chempotwann' )
   call de_alloc( numhg,     name='numhg',     routine= 'chempotwann' )
   call de_alloc( H_chempotwann_full, name='H_chempotwann_full', & 
 &                routine= 'chempotwann' )
#endif

        
!  Stop time counter
   call timer( 'chempotwann', 2 )

   end subroutine chempotwann


   subroutine add_Hamiltonian_chempotwann( H_chempotwann )
   implicit none 

   real(dp), intent(in) :: H_chempotwann(maxnh,spin%H)   
                                            ! Extra term in the Hamiltonian 
                                            !   that accounts for a 
                                            !   rigid shift of the 
                                            !   bands associated
                                            !   with a given Wannier

   integer :: io_local  ! Counter for the local atomic orbitals in a given node
   integer :: j         ! Counter of neighbours of a given orbital
   integer :: ispin     ! Counter of spin components
   integer :: ind       ! Index of the neighbour orbital in listh

!  Add the Hamiltonian elements
   do io_local = 1, no_l
     do j = 1, numh(io_local)
       ind = listhptr(io_local) + j
       do ispin = 1, spin%H
         H(ind,ispin) = H(ind,ispin) + H_chempotwann(ind,ispin)
       enddo 
     enddo 
   enddo 

   end subroutine add_Hamiltonian_chempotwann
!
  end module m_chempotwann
