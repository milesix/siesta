! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
  subroutine broadcast_w90_in_siesta()
!
! Globalizes the data related with the Wannier transformation
!
! Javier Junquera, December 2018, based on broadcast_basis written by
! Alberto Garcia, June 2000--

  use parallel,  only : Node, Nodes
  use files,     only: label_length       ! Number of characters in slabel
  use siesta_options, only: n_wannier_manifolds
                                          ! Number of bands manifolds
                                          !   that will be considered
                                          !   for Wannier transformation
  use w90_in_siesta_types,   only: manifold_bands_w90_in
                                          ! Variable where the initial
                                          !   and final band of each
                                          !   manifold are stored
  use w90_in_siesta_types,   only: compute_chempotwann
                                          ! Compute the Hamiltonian matrix
                                          !   elements between NAO
                                          !   if a chemical potential in 
                                          !   a Wannier is applied
  use w90_in_siesta_types,   only: first_chempotwann
                                          ! Is this the first time that 
                                          !   that the subroutines to call
                                          !   the chemical potential 
                                          !   are called?
  use w90_in_siesta_types,   only: chempotwann_val
                                          ! Chemical potential
                                          !   applied to shift the energy of 
                                          !   the matrix elements in real space
                                          !   associated with a given
                                          !   Wannier function
  use atomlist,       only: no_u          ! Number of orbitals in unit cell
                                          ! NOTE: When running in parallel,
                                          !   this is core independent
!
! Allocation/Deallocation routines
!
  use alloc,          only: re_alloc      ! Reallocation routines
  use alloc,          only: de_alloc      ! Deallocation routines

      
#ifdef MPI
  use mpi_siesta
  use m_mpi_utils, only: broadcast
#endif

  implicit none

#ifndef MPI
!
! Do nothing...
!
  end subroutine broadcast_w90_in_siesta
#else
  integer MPIerror

  integer i_man     ! Counter for loops on manifolds
  integer iorb      ! Counter for loops on orbitals
  integer nwannier90! Number of bands considered for Wannier transformation

#ifdef DEBUG
  call write_debug( '  PRE broadcast_w90_in_siesta' )
#endif

  if (Nodes.eq.1) return

  if (Node.ne.0) then
    allocate(manifold_bands_w90_in(n_wannier_manifolds))
  end if

  do i_man = 1, n_wannier_manifolds
    call MPI_Bcast(manifold_bands_w90_in(i_man)%seedname_w90_in,              &
 &                 label_length+3,MPI_character,0,                            &
 &                 MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%initial_band,                 &
 &                 1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%final_band,                   &
 &                 1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%number_of_bands,              &
 &                 1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%numbands_w90_in,              &
 &                 1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    nwannier90 = manifold_bands_w90_in(i_man)%numbands_w90_in
    call MPI_Bcast(manifold_bands_w90_in(i_man)%num_iter,                     &
 &                 1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%blocksizeincbands_w90_in,     &
 &                 1,MPI_Integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%nincbands_loc_w90_in,         &
 &                 1,MPI_Integer,0,MPI_Comm_World,MPIerror)

    call MPI_Bcast(manifold_bands_w90_in(i_man)%dis_win_min,                  &
 &                 1,MPI_double_precision,                                    &
 &                 0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%dis_win_max,                  &
 &                 1,MPI_double_precision,                                    &
 &                 0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%dis_froz_min,                 &
 &                 1,MPI_double_precision,                                    &
 &                 0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%dis_froz_max,                 &
 &                 1,MPI_double_precision,                                    &
 &                 0,MPI_Comm_World,MPIerror) 
    call MPI_Bcast(manifold_bands_w90_in(i_man)%disentanglement,              &
 &                 1,MPI_logical,                                             &
 &                 0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%wannier_plot,                 &
 &                 1,MPI_logical,                                             &
 &                 0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%wannier_plot_supercell,       & 
 &                 3,MPI_integer,                                             &
 &                 0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%fermi_surface_plot,           &
 &                 1,MPI_logical,                                             &
 &                 0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%write_hr,                     &
 &                 1,MPI_logical,                                             &
 &                 0,MPI_Comm_World,MPIerror)
    if (Node .ne. 0) then
      nullify( manifold_bands_w90_in(i_man)%orbital_indices )
      nullify( manifold_bands_w90_in(i_man)%isexcluded      )
      nullify( manifold_bands_w90_in(i_man)%orbexcluded     )
      nullify( manifold_bands_w90_in(i_man)%orb_in_manifold )

      call re_alloc( manifold_bands_w90_in(i_man)%orbital_indices,            &
 &                   1, manifold_bands_w90_in(i_man)%numbands_w90_in )
      call re_alloc( manifold_bands_w90_in(i_man)%isexcluded,                 &
 &                   1, no_u )
      call re_alloc( manifold_bands_w90_in(i_man)%orbexcluded,                &
 &                   1, no_u )
      call re_alloc( manifold_bands_w90_in(i_man)%orb_in_manifold,            &
 &                   1, no_u )
    end if

    call MPI_Bcast(manifold_bands_w90_in(i_man)%orbital_indices,              &
 &                 nwannier90,MPI_integer,0,                                  &
 &                 MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%isexcluded,                   &
 &                 no_u,MPI_logical,0,                                        &
 &                 MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%orbexcluded,                  &
 &                 no_u,MPI_logical,0,                                        &
 &                 MPI_Comm_World,MPIerror)
    call MPI_Bcast(manifold_bands_w90_in(i_man)%orb_in_manifold,              &
 &                 no_u,MPI_integer,0,                                        &
 &                 MPI_Comm_World,MPIerror)

  enddo

  call MPI_Bcast(compute_chempotwann,                                         &
 &                 1,MPI_logical,                                             &
 &                 0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(first_chempotwann,                                           &
 &                 1,MPI_logical,                                             &
 &                 0,MPI_Comm_World,MPIerror)

  if (Node.ne.0) then
    allocate(chempotwann_val(manifold_bands_w90_in(1)%numbands_w90_in))
  end if
  call MPI_Bcast(chempotwann_val,                                            &
 &               manifold_bands_w90_in(1)%numbands_w90_in,                   &
 &               mpi_double_precision,0,                                     &
 &               MPI_Comm_World,MPIerror)

#ifdef DEBUG
  call write_debug( '  POS broadcast_w90_in_siesta' )
#endif


!! For debugging
!  write(6,'(a,2i5)')                                                         &
! & 'broadcast_w90_in_siesta: Node, n_wannier_manifolds = ',                  &
! &  Node, n_wannier_manifolds 
!  do i_man = 1, n_wannier_manifolds
!    nwannier90 = manifold_bands_w90_in(i_man)%numbands_w90_in
!    write(6,'(a,2i5,2x,a)')                                                  &
! &   'broadcast_w90_in_siesta: Node, i_manifold, seedname        = ',        &
! &   Node, i_man, manifold_bands_w90_in(i_man)%seedname_w90_in
!    write(6,'(a,3i5)')                                                       &
! &   'broadcast_w90_in_siesta: Node, i_manifold, initial_band    = ',        &
! &   Node, i_man, manifold_bands_w90_in(i_man)%initial_band
!    write(6,'(a,3i5)')                                                       &
! &   'broadcast_w90_in_siesta: Node, i_manifold, final_band      = ',        &
! &   Node, i_man, manifold_bands_w90_in(i_man)%final_band
!    write(6,'(a,3i5)')                                                       & 
! &   'broadcast_w90_in_siesta: Node, i_manifold, number_of_bands = ',        &
! &   Node, i_man, manifold_bands_w90_in(i_man)%number_of_bands
!    write(6,'(a,3i5)')                                                       & 
! &    'broadcast_w90_in_siesta: Node, i_manifold, numbands_w90_in = ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%numbands_w90_in
!    write(6,'(a,3i5)')                                                       &
! &    'broadcast_w90_in_siesta: Node, i_manifold, num_iter        = ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%num_iter      
!    write(6,'(a,3i5)')                                                       &
! &    'broadcast_w90_in_siesta: Node, i_manifold, blocksize       = ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%blocksizeincbands_w90_in
!    write(6,'(a,3i5)')                                                       &
! &    'broadcast_w90_in_siesta: Node, i_manifold, nincbands_loc   = ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%nincbands_loc_w90_in
!    do iorb = 1, nwannier90
!      write(6,'(a,4i5)')                                                     &
! &     'broadcast_w90_in_siesta: Node, i_manifold, orbital_indices = ',      &
! &     Node, i_man, iorb,                                                    &
! &     manifold_bands_w90_in(i_man)%orbital_indices(iorb)
!    enddo
!    do iorb = 1, no_u
!      write(6,'(a,3i5,l5)')                                                  &
! &     'broadcast_w90_in_siesta: Node, i_manifold, isexcluded      = ',      &
! &     Node, i_man, iorb,                                                    &
! &     manifold_bands_w90_in(i_man)%isexcluded(iorb)
!    enddo
!    do iorb = 1, no_u
!      write(6,'(a,3i5,l5)')                                                  &
! &     'broadcast_w90_in_siesta: Node, i_manifold, orbexcluded     = ',      &
! &     Node, i_man, iorb,                                                    &
! &     manifold_bands_w90_in(i_man)%orbexcluded(iorb)
!    enddo
!    do iorb = 1, no_u
!      write(6,'(a,4i5)')                                                     &
! &     'broadcast_w90_in_siesta: Node, i_manifold, orb_in_manifold = ',      &
! &     Node, i_man, iorb,                                                    &
! &     manifold_bands_w90_in(i_man)%orb_in_manifold(iorb)
!    enddo
!    write(6,'(a,2i5,f12.5)')                                                 &
! &    'broadcast_w90_in_siesta: Node, i_manifold, dis_win_min     = ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%dis_win_min       
!    write(6,'(a,2i5,f12.5)')                                                 &
! &    'broadcast_w90_in_siesta: Node, i_manifold, dis_win_max     = ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%dis_win_max
!    write(6,'(a,2i5,f12.5)')                                                 &
! &    'broadcast_w90_in_siesta: Node, i_manifold, dis_froz_min    = ',       &
! &   Node, i_man, manifold_bands_w90_in(i_man)%dis_froz_min
!    write(6,'(a,2i5,f12.5)')                                                 &
! &    'broadcast_w90_in_siesta: Node, i_manifold, dis_froz_max    = ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%dis_froz_max
!    write(6,'(a,2i5,l5)')                                                    &
! &    'broadcast_w90_in_siesta: Node, i_manifold, disentanglement = ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%disentanglement
!    write(6,'(a,2i5,l5)')                                                    &
! &    'broadcast_w90_in_siesta: Node, i_manifold, wannier_plot    = ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%wannier_plot
!    write(6,'(a,5i5)')                                                       &
! &    'broadcast_w90_in_siesta: Node, i_manifold, wannier_plot_su = ',       &
! &    Node, i_man,                                                           &
! &    manifold_bands_w90_in(i_man)%wannier_plot_supercell(:)
!    write(6,'(a,2i5,l5)')                                                    &
! &    'broadcast_w90_in_siesta: Node, i_manifold, fermi_surface_pl= ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%fermi_surface_plot
!    write(6,'(a,2i5,l5)')                                                    &
! &    'broadcast_w90_in_siesta: Node, i_manifold, write_hr        = ',       &
! &    Node, i_man, manifold_bands_w90_in(i_man)%write_hr
!  enddo 
!! End debugging

  end subroutine broadcast_w90_in_siesta

#endif

