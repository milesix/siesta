
      subroutine alloc_error_report(str,code)

      character(len=*), intent(in)  :: str
      integer, intent(in)  :: code

      external :: die
      
      if (code == 0) then
        call die(str)
      else
         write(0,*) trim(str)
         write(6,*) trim(str)
      endif

      end subroutine alloc_error_report

! The PSML library calls a "die" routine when it encounters an
! error. This routine should take care of carrying out any needed
! cleaning and terminating the program.  As the details would vary with
! the client program, each program has to provide its own.
! 

      subroutine psml_die(str)
      character(len=*), intent(in) :: str

      external :: die
      
      call die(str)

      end subroutine psml_die
