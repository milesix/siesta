module radial_log_low_level
use precision, only : dp
use radial_logGrid, only : logGrid_t, log_grid_get_r, log_grid_copy,&
     log_grid_dealloc
use flib_spline, only : evaluate_spline, generate_spline
use sys, only : die
implicit none


  real(dp),public, parameter                :: eps=1.0d-6
  real(dp),public, parameter                :: min_func_val=1E-6
  integer, public :: maximum_length = 0
  logical, public :: restricted_grid

  type log_Rad_Func_t
     real(dp), dimension(:), pointer :: f   ! Actual data
     real(dp), dimension(:), pointer :: d2  ! 2nd Derivative of data
     type(logGrid_t), pointer        :: grid         
  end type log_Rad_Func_t

  public log_rad_func_t, log_rad_alloc, &
       log_rad_dealloc, log_rad_cutoff, log_rad_get_default_length, &
       log_rad_setup_d2, log_rad_set_default_length, log_rad_copy_grid, &
       log_rad_set_maximum_length, log_rad_min, log_rad_get, &
       log_rad_set_origin, log_rad_update, log_rad_copy, log_rad_get_length

  private

  contains

    subroutine log_rad_alloc(rad_func,values,grid,yp1)
    type (log_rad_func_t), intent(out) :: rad_func
    real(dp), intent(in) :: values(:)
    type(logGrid_t), intent(in) :: grid
    real(dp), optional, intent(in) :: yp1

    integer :: length
    length = size(values)
    
    allocate(rad_func%f(1:length),rad_func%d2(1:length))
    allocate(rad_func%grid)
    rad_func%grid = grid

    rad_func%f = values

    if (present(yp1))then
       call log_rad_setup_d2(rad_func,yp1)
    else
       call log_rad_setup_d2(rad_func)
    endif
    
  end subroutine log_rad_alloc

  !--------------------------------------------

  subroutine log_rad_dealloc(func)
    type(log_rad_func_t), intent(inout) :: func
    !func%n = 0
    !func%cutoff = 0.0_dp
    deallocate(func%f,func%d2)
    call log_grid_dealloc(func%grid)
    deallocate(func%grid) 
  end subroutine log_rad_dealloc

  !--------------------------------------------
  function log_rad_cutoff(func)
    real(dp)                         :: log_rad_cutoff
    type(log_rad_func_t), intent(in) :: func
    log_rad_cutoff = log_grid_get_r(func%grid,size(func%f))    
  End function log_rad_cutoff

 !--------------------------------------------

  function log_rad_get_default_length() result(length)
    !type(log_rad_func_t), intent(in) :: func
    integer :: length
    length = maximum_length
    !print *, "WARNING: radial_log: using log_rad_default_length"
  end function log_rad_get_default_length

  !------------------------------------------

  subroutine log_rad_set_default_length(length) 
    integer :: length
    if (length > maximum_length) maximum_length = length
  end subroutine  log_rad_set_default_length

  !------------------------------------------

  
  subroutine log_rad_copy_grid(src,dest)
    !     Copy src into dest.
    type (log_Rad_Func_t), intent(in) :: src
    type (log_Rad_Func_t), intent(out):: dest

    !call log_Rad_Alloc(dest,size(src%f))
    !call log_grid_set_ab(dest%grid,src%grid%a,src%grid%b)
    !dest%cutoff = src%cutoff
    call log_grid_copy(src%grid,dest%grid)
  end subroutine log_rad_copy_grid

!----------------------------------

  subroutine log_rad_set_maximum_length(length)
    integer, intent(in) :: length
    maximum_length = length
  end subroutine log_rad_set_maximum_length

!----------------------------------

  function log_rad_min(func) result(min)
    type(log_rad_func_t), intent(in) ::func
    real(dp) :: min
    min = minval(func%f)
  end function  log_rad_min

!----------------------------------
  subroutine log_rad_get(func,r,fr,dfdr)
    type(log_rad_func_t), intent(in) :: func
    real(dp), intent(in)         :: r
    real(dp), intent(out)        :: fr
    real(dp), optional,intent(out)        :: dfdr

    real(dp) :: df
    if (size(func%f) .eq. 0) then
       fr = 0._dp
       dfdr = 0._dp
    else
       if (r.gt.log_rad_cutoff(func)) then
          call die('radial: Attempt to evaluate beyond cutoff')
          fr = 0._dp
          if (present(dfdr)) dfdr = 0._dp
       else
          call evaluate_spline(func%grid%r,func%f,func%d2,size(func%f),r,fr,df)
          if (present(dfdr)) dfdr = df
       endif
    endif
  end subroutine log_rad_get

  !---------------------------------

 subroutine log_rad_set_origin(func,value)
    type(log_rad_func_t), intent(inout) :: func
    real(dp),             intent(in) :: value

    func%f(1) = value
    call log_rad_setup_d2(func)
  end subroutine log_rad_set_origin

  !--------------------------------------------

  !------------------------------------------------

  subroutine log_rad_setup_d2(func, yp1)
    type(log_rad_func_t) :: func
    real(dp),          intent(in),optional  :: yp1
    if (present(yp1)) then
       call generate_spline(func%grid%r,func%f,size(func%f),func%d2,YP1,0.0_dp)
    else
       call generate_spline(func%grid%r,func%f,size(func%f),func%d2)
    endif
  end subroutine log_rad_setup_d2

  !-------------------------------------------------

  function log_rad_update(func) result(new_func)
    type(log_rad_func_t), intent(in) :: func
    type(log_rad_func_t) :: new_func
    integer :: nr_new,ir, old_nr
    real(dp) :: val
    real(dp), pointer :: new_values(:)
    type(logGrid_t) :: grid
    nr_new = 0    
    old_nr = size(func%f)

    do ir=old_nr,2,-1
       val=func%f(ir)
       if((abs(val).gt.eps).and.(nr_new.eq.0)) then            
          if (ir==old_nr) then
             nr_new=old_nr
          else
             nr_new=ir+1
          endif
          exit
       endif
    enddo
    allocate(new_values(1:nr_new))
    grid = func%grid
    new_values(1:nr_new) = func%f(1:nr_new)
    call log_rad_alloc(new_func,new_values,grid)

    deallocate(new_values)
  end function log_rad_update

  !------------------------------------------------

  subroutine log_rad_copy(src,dest)
    type(log_rad_func_t), intent(in) :: src
    type(log_rad_func_t), intent(out) :: dest
    
    call log_rad_alloc(dest,src%f,src%grid)

  end subroutine log_rad_copy

  !---------------------------------------

  !function log_rad_get_grid(func) result(grid)
  !  type(log_rad_func_t), intent(in) :: func
  !  type(logGrid_t)                  :: grid
  !  grid=func%grid
  !end function log_rad_get_grid

  !-----------------------------------------

  function log_rad_get_length(func) result(length)
    type(log_rad_func_t), intent(in) :: func
    integer :: length
    length = size(func%f)
  end function log_rad_get_length

  !----------------------------------------

end module radial_log_low_level
