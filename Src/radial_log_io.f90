module radial_log_io
  use precision, only : dp
  use radial_log_low_level, only : log_rad_func_t
  use flib_spline, only: generate_spline
  implicit none

  public log_rad_dump_file,log_rad_dump_funcs_ascii, &
       log_rad_dump_ascii, log_rad_read_ascii_unformatted, &
       log_rad_read_ascii_formatted,log_rad_write_ascii_formatted

  private
 
contains

  subroutine log_rad_dump_file(func,label)
    type(log_rad_func_t), intent(in) :: func
    character(len=*), intent(in) :: label
    
    integer :: ir,io
    
    io = 1024
    open(unit=io,file=label,status="replace")
    do ir=1,size(func%f)
       write(io,*) func%grid%r(ir),func%f(ir),func%d2(ir)
    end do
    close(io)

  end subroutine log_rad_dump_file

  !-----------------------------------------

   subroutine log_rad_dump_funcs_ascii(func1,func2,io)
    type(log_rad_func_t), intent(in) :: func1, func2
    integer, intent(in) :: io
    
    integer :: j
    
 9040    format(i4,1x,6f12.6)

    do j = 1, size(func1%Grid%r)
       write(io,9040) j, func1%Grid%r(j),func1%f(j),func2%f(j)
    enddo
  end subroutine log_rad_dump_funcs_ascii

  !-----------------------------------------------------

  subroutine log_rad_dump_ascii(func1,io)
    type(log_rad_func_t), intent(in) :: func1
    integer, intent(in) :: io
    
    integer :: j
    do j = 1, size(func1%f)
       write(io,'(2f20.7)') func1%Grid%r(j),func1%f(j)
    enddo
  end subroutine log_rad_dump_ascii

  !------------------------------------------------------

  subroutine log_rad_read_ascii_unformatted(io,func)
    integer, intent(in) :: io
    type(log_rad_func_t), intent(inout) :: func

    integer :: j
    read(io) (func%f(j),j=2,size(func%f))
    func%f(1)=func%f(2)
    call generate_spline(func%grid%r,func%f,size(func%f),func%d2,0.0_dp,0.0_dp)
    
  end subroutine log_rad_read_ascii_unformatted

  !---------------------------------------------
  subroutine log_rad_read_ascii_formatted(io,func)
    integer, intent(in) :: io
    type(log_rad_func_t), intent(inout) :: func

    integer :: j
   
8030 format(4(g20.12))
    read(io,8030) (func%f(j),j=2,size(func%f))
    func%f(1)=func%f(2) 
    call generate_spline(func%grid%r,func%f,size(func%f),func%d2,0.0_dp,0.0_dp)
    
  end subroutine log_rad_read_ascii_formatted

  !---------------------------------------------

  subroutine log_rad_write_ascii_formatted(io,func)
    integer, intent(in) :: io
    type(log_rad_func_t), intent(in) :: func

    integer :: j
   
8030 format(4(g20.12))
    write(io,8030) (func%f(j), j=2,size(func%f))
    
  end subroutine log_rad_write_ascii_formatted

  !---------------------------------------------
  
  
end module radial_log_io
