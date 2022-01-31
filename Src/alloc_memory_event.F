!--------------------------------------------------
! Stand-alone routine to capture memory events from
! the alloc module
!
! Each program using the alloc module needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta
!
!   - The use of the memory_log module for reports
!
!------------------------------------------------------

!
!     Version for libgridxc handler ... to be
!     deprecated when library moves to using
!     procedure pointers for handlers
!
      subroutine alloc_memory_event(bytes,name)
      use memory_log, only: memory_event

      integer, intent(in)           :: bytes
      character(len=*), intent(in)  :: name

      call memory_event(bytes,name)

      end subroutine alloc_memory_event


!---------------------------------------------------
!     Version for Siesta
!
      subroutine alloc_memory_event_long(bytes,name)
      use memory_log, only: memory_event

      integer, parameter :: i8b = selected_int_kind(18)

      integer(i8b), intent(in)      :: bytes
      character(len=*), intent(in)  :: name

      call memory_event(bytes,name)

      end subroutine alloc_memory_event_long
