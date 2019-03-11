!
!   Broadcast complete fdf structure
!
    SUBROUTINE broadcast_fdf_struct(reading_node,comm)

#ifndef MPI
      implicit none

      integer, intent(in)       :: reading_node         ! Node which contains the struct
      integer, intent(in)       :: comm
#else
      use mpi

      use fdf, only: serialize_fdf_struct, recreate_fdf_struct
      use fdf, only: set_fdf_started
      
      implicit none

      integer, intent(in)       :: reading_node         ! Node which contains the struct
      integer, intent(in)       :: comm

      character, allocatable    :: bufferFDF(:)
      integer                   :: ierr, nchars, rank

      call MPI_Comm_Rank( Comm, rank, ierr )
      
      if (rank == reading_node) then
         call serialize_fdf_struct(bufferFDF)
         nchars = size(bufferFDF)
      endif

      call MPI_Bcast(nchars, 1,                                 &
                     MPI_INTEGER, reading_node, comm, ierr)
      if (ierr .ne. MPI_SUCCESS) then
        call die("Error broadcasting size of fdf struct")
      endif

      if (rank /= reading_node) then
         ALLOCATE(bufferFDF(nchars), stat=ierr)
         if (ierr .ne. 0) then
            call die("Error allocating buffer in non-root node")
         endif
      endif

      call MPI_Bcast(bufferFDF, nchars,              &
                     MPI_CHARACTER, reading_node, MPI_COMM_WORLD, ierr)
      if (ierr .ne. MPI_SUCCESS) then
        call die("Error Broadcasting bufferFDF")
      endif

      if (rank /= reading_node) then
         call recreate_fdf_struct(bufferFDF)
         call set_fdf_started(.true.)
      endif

      DEALLOCATE(bufferFDF)
#endif
    END SUBROUTINE broadcast_fdf_struct
