      subroutine die(str)

      character(len=*), intent(in) :: str

      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)

      stop 
      end subroutine die
!
      subroutine bye(str)

      character(len=*), intent(in)  :: str

      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)

      stop
      
      end subroutine bye

      subroutine message(level,str)

      ! One of INFO, WARNING, FATAL
      character(len=*), intent(in)  :: level
      character(len=*), intent(in)  :: str

      write(6,"(a)") trim(level) // ": " // trim(str)
      write(0,"(a)") trim(level) // ": " // trim(str)

      end subroutine message
     
