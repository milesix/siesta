!
! These are 'external' versions of the routines in module sys
!
subroutine die(str)
  use sys, only: die_sys => die
  character(len=*), intent(in) :: str
  call die_sys(str)
end subroutine die

subroutine bye(str)
  use sys, only: bye_sys => bye
  character(len=*), intent(in) :: str
  call bye_sys(str)
end subroutine bye

subroutine message(level,str)
  use sys, only: message_sys => message
  character(len=*), intent(in) :: level
  character(len=*), intent(in) :: str
  call message_sys(level, str)
end subroutine message

subroutine reset_messages_file()
  use sys, only: reset_messages_file_sys => reset_messages_file
  call reset_messages_file_sys()
end subroutine reset_messages_file


