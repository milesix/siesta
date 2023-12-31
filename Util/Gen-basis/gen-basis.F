! 
! Copyright (C) 1996-2014	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      program gen_basis

C Stand-alone program to generate the PAOs, KB projectors, Vlocal, etc
C for a set of atoms. 
C
C The input is exactly the same as that for SIESTA.
C This program generates .ion (optionally .ion.nc) files readable
C by SIESTA
C
C RESTRICTION: Just a single species allowed
C
!!      use f2kcli
      use fdf

      use precision
      use parallel,     only : Node, Nodes, IOnode
      use basis_types, only: nsp, basis_specs_transfer
      use basis_types, only: iz, lmxkb, nkbl, 
     .           erefkb, lmxo, nzeta, rco, 
     .           lambda,
     .           atm_label, polorb, semic, nsemic,
     .           cnfigmx, charge, smass, basistype,
     $           rinn, vcte, split_norm
      use basis_types, only: write_basis_specs
      use basis_types, only: basis_def_t, basis_parameters
      use basis_specs, only: read_basis_specs
      use basis_io
      use old_atmfuncs, only : nsmax
      use old_atmfuncs, only : clear_tables, allocate_old_arrays

      use atom,         only : atom_main, prinput
      use xcmod,        only : setXC
      use sys,          only : die
      use atom_options, only : get_atom_options, write_ion_plot_files

      implicit none

      integer is
!     integer narg
      character(len=50) filein

      type(basis_def_t), pointer   :: basp

      external atm_transfer

c Reading input for the pseudopotentials and atomic orbitals 

      Node =  0
      Nodes = 1

      IOnode = (Node .eq. 0)

      if (Ionode) then
!
         filein = "stdin"
         call fdf_init(filein,'gen-basis_fdf.log')

         call get_atom_options()
         write_ion_plot_files = .true.  ! get them always

         call setXC

         call read_basis_specs()
         call basis_specs_transfer()

         nsmax = nsp            ! For old_atmfuncs
         call allocate_old_arrays
         call clear_tables

         if (nsp .gt. 1)
     $        call die("gen-basis: JUST A SINGLE SPECIES ALLOWED...")

         do is = 1,nsp
            call write_basis_specs(6,is)
            basp=>basis_parameters(is)
            call atom_main( iz(is), lmxkb(is),
     $           nkbl(0:,is), erefkb(1:,0:,is),lmxo(is),
     $           nzeta(0:,1:,is), rco(1:,0:,1:,is), 
     $           lambda(1:,0:,1:,is), atm_label(is),
     $           polorb(0:,1:,is), semic(is), nsemic(0:,is),
     $           cnfigmx(0:,is),charge(is),
     $           smass(is), basistype(is), is,
     $           rinn(0:,1:,is), vcte(0:,1:,is),
     $           split_norm(0:,1:,is), basp)
         enddo 

         call prinput(nsp)

!        Create the new data structures

         call atm_transfer

         call dump_basis_ascii
         call dump_basis_netcdf

      endif

      end program gen_basis
!






