! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
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
      use atom,         only : generate_all_atomic_info
C
C     It has to be MPI-aware since it calls routines down the line...
C
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

!     integer narg
      character(len=50) filein

#ifdef MPI
      integer  MPIerror
#endif
      external atm_transfer

c Reading input for the pseudopotentials and atomic orbitals 
C Initialise MPI and set processor number
#ifdef MPI
      call MPI_Init( MPIerror )
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      Node =  0
      Nodes = 1
#endif

      IOnode = (Node .eq. 0)

      if (Ionode) then
!
!        Problems with mpich command-line argument handling...
!$$$         narg = command_argument_count()
!$$$         if (narg == 0) then
!$$$            filein = "stdin"
!$$$         else if (narg == 1) then
!$$$            call get_command_argument(1,filein)
!$$$         else
!$$$            stop "too many command-line arguments"
!$$$         endif
!
         filein = "stdin"
         call fdf_init(filein,'gen-basis_out.fdf')
       
!        Create the new data structures
         
         call generate_all_atomic_info()
      endif

      end program gen_basis
!






