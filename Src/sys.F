! ---
! Copyright (C) 1996-2014	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! 
!
      module sys
!
!     Termination and messaging routines, MPI aware
!
      public :: die, bye, message

      CONTAINS

      subroutine message(str)

      use parallel, only : Node
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str


      if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
      endif

      end subroutine message
!
!--------------------------------------------------
      subroutine die(str)

      use parallel, only : Node
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str

#ifdef MPI
      integer MPIerror
#endif

! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

!!                                       if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
            write(0,'(a)') trim(str)
         endif
         write(6,'(a,i4)') 'ERROR STOP from Node: ', Node
         write(0,'(a,i4)') 'ERROR STOP from Node: ', Node
!!                                       endif
         if (Node .eq. 0) then
            call pxfflush(6)
            call pxfflush(0)
         endif

#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
      stop
#else
      call pxfabort()
#endif
      end subroutine die

!---------------------------------------------------------
      subroutine bye(str)

      use parallel, only : Node
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str

#ifdef MPI
      integer rc
#endif


! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

!!                                       if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a)') 'User-Requested End of Run. Bye!!'
!!                                       endif

         if (Node .eq. 0) call pxfflush(6)

#ifdef MPI
      call MPI_Finalize(rc)
      stop
#else
      stop
#endif
      end subroutine bye

      end module sys

