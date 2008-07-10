module radial_logGrid
  use precision, only : dp
  use sys, only: die
  implicit none

  type logGrid_t
     real(dp)                        :: a,b
     real(dp), dimension(:), pointer :: r   ! Log mesh
     real(dp), dimension(:), pointer :: drdi! 2nd derivative of logmesh
  end type logGrid_t

  public logGrid_t, log_grid_alloc, log_grid_get_r, log_grid_dealloc, &
       log_grid_copy, log_grid_get_number_of_points, log_grid_get_a, &
       log_grid_get_b, log_grid_dump_ascii_formatted, log_grid_get_ir

  private

  contains
    
  !--------------------------------------------

  function log_grid_alloc(n,a,b,r) result(grid)
    integer,intent(in)          :: n
    real(dp), intent(in)        :: a, b,r(1:n)
    type(logGrid_t)             :: grid

    !Internal vars
    real(dp) :: rpb, ea
    integer :: ir

    allocate(grid%r(1:n),grid%drdi(1:n))
    grid%a = a
    grid%b = b
    grid%r = r

    !    Calculate drdi
    !    drdi is the derivative of the radial distance respect to the mesh index
    !    i.e. rofi(ir)= b*[ exp( a*(i-1) ) - 1 ] and therefore 
    !    drdi=dr/di =a*b*exp(a*(i-1))= a*[rofi(ir)+b] 

    rpb=grid%b
    ea=exp(grid%a)
    
    do ir=1,size(grid%r)
       Grid%drdi(ir)=grid%a*rpb
       rpb=rpb*ea
    enddo

  end function log_grid_alloc

  !-------------------------------------------
  
  subroutine log_grid_dealloc(grid)
    type(logGrid_t),intent(inout) :: grid

    nullify(grid%r,grid%drdi)
    grid%a = 0.0_dp
    grid%b = 0.0_dp
  end subroutine log_grid_dealloc

  !-------------------------------------------
    
  subroutine log_grid_zero(grid)
    type(logGrid_t),intent(out) :: grid
    !grid%n = 0.0_dp
    grid%a = 0.0_dp
    grid%b = 0.0_dp
    grid%r = 0.0_dp
    grid%drdi = 0.0_dp
  end subroutine log_grid_zero

  !---------------------------------------------

  subroutine log_grid_copy(src,dest)
    type(logGrid_t), intent(in) :: src
    type(logGrid_t), intent(out):: dest

    dest = log_grid_alloc(size(src%r),src%a,src%b,src%r)
    
  end subroutine log_grid_copy

  !---------------------------------------------

  function log_grid_get_number_of_points(grid) result(length)
    type(logGrid_t), intent(in) :: grid
    integer                     :: length
    length = size(grid%r)
  end function log_grid_get_number_of_points

  !---------------------------------------------

  function log_grid_get_r(grid,ir) result(r)
    type(logGrid_t), intent(in) :: grid
    integer,intent(in)          :: ir

    real(dp) :: r
    r=0.0_dp
    if (ir .le. size(grid%r)) then
       r=grid%r(ir)
    else
       call die("log_grid_get_r: ir too big!")
    endif
    
  endfunction log_grid_get_r

  !---------------------------------------------

  function log_grid_get_ir(grid,r) result (ir)
    type(logGrid_t), intent(in) :: grid
    real(dp), intent(in) :: r
    
    integer :: ir
   
    ir = nint(log(r/grid%b+1.0_dp)/grid%a)+1

  endfunction log_grid_get_ir

  !-------------------------------------------------

  function log_grid_get_a(grid) result(a)
    type(logGrid_t), intent(in) :: grid
    real(dp) :: a
    a = grid%a
  endfunction log_grid_get_a

  !--------------------------------------------------

  function log_grid_get_b(grid) result(b)
    type(logGrid_t), intent(in) :: grid
    real(dp) :: b
    b = grid%b
  endfunction log_grid_get_b

  !--------------------------------------------------

  subroutine log_grid_dump_ascii_formatted(grid,io_ps)
    type(logGrid_t), intent(in) :: grid
    integer,         intent(in) :: io_ps
    integer :: j
    write(io_ps,8030) (Grid%r(j),j=2,size(grid%r))
8030 format(4(g20.12))

  end subroutine log_grid_dump_ascii_formatted

end module radial_logGrid
