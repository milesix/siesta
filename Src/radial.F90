module radial

  use precision
  use m_recipes, only: spline,splint
  use radial_lin
  use radial_log
  use radial_log_dirty
  use radial_log_low_level
  use radial_log_io
  use radial_logGrid, only : logGrid_t, log_grid_get_a, log_grid_get_b, &
       log_grid_dump_ascii_formatted, log_grid_get_ir
  use sys

  implicit none

  type :: rad_kind_t
     private
     integer :: kind = 0
  end type rad_kind_t

  type(rad_kind_t), parameter, public :: &
       log_t = rad_kind_t(1), &
       lin_t = rad_kind_t(2)

  type, public :: rad_func_t
     private
     type (rad_kind_t)              :: kind
     type (lin_rad_func_t), pointer :: lin
     type (log_rad_func_t), pointer :: log
  end type rad_func_t

  type :: rad_grid_t
     private
     type (rad_kind_t)         :: kind
     type (logGrid_t), pointer :: log
     type (linGrid_t), pointer :: lin
  end type rad_grid_t

  interface operator (==)
     module procedure rad_kind_compare
  end interface

  real(dp), public ::rmax_radial_grid


contains

  !-------------------------------------------------------------------------

  subroutine rad_alloc(rad_func,values,grid,yp)
    !Allocate a given function
    type (rad_func_t), intent(out)   :: rad_func
    real(dp), dimension(:), intent(in) :: values
    type (rad_grid_t), intent(in)      :: grid
    real(dp), optional, intent(in)     :: yp

    if (grid%kind == lin_t) then
       allocate(rad_func%lin)
       call lin_rad_alloc(rad_func%lin,grid%lin%delta,values)
       rad_func%kind = lin_t
    elseif(grid%kind == log_t) then
       allocate (rad_func%log)
       if (present(yp)) then
          call log_rad_alloc(rad_func%log, values,grid%log,yp)
       else
          call log_rad_alloc(rad_func%log, values,grid%log)
       endif
       rad_func%kind = log_t
    else
       call die("radial: rad_alloc unknown type!")
    endif
  end subroutine rad_alloc

  !-------------------------------------------------------------------------

  function rad_allocated(rad_func) result(allocated)
    !Check if the function is allocated.
    type(rad_func_t), intent(in) :: rad_func
    logical :: allocated

    allocated = .false.
    if(associated(rad_func%log) .or. associated(rad_func%lin)) allocated = .true.

  end function rad_allocated

  !-------------------------------------------------------------------------

  subroutine rad_broadcast(rad_func)
    ! Broadcast (MPI)
    use parallel, only : Node,  nodes

#ifdef MPI
    use mpi_siesta
#endif

    type (rad_func_t), intent(inout)   :: rad_func
  
#ifndef MPI
  end subroutine rad_broadcast
#else

  if (Nodes .eq. 1) return

  call rad_kind_broadcast(rad_func)

  if (rad_func%kind == lin_t) then
     if(Node .ne. 0) allocate(rad_func%lin)
     call lin_rad_broadcast(rad_func%lin)
     rad_func%kind = lin_t
  elseif(rad_func%kind == log_t) then
     if(Node .ne. 0) allocate(rad_func%log)
     call log_rad_broadcast(rad_func%log)
     rad_func%kind = log_t
  else
     call die("radial: broadcast unknown type!")
  endif
  end subroutine rad_broadcast

#endif

  !------------------------------------------------------------------------- 

  subroutine rad_copy(src,dest)
    !Copy from src to dest.
    type (rad_func_t), intent(in)   :: src
    type (rad_func_t), intent(out)  :: dest

    if (src%kind == lin_t) then
       allocate(dest%lin)
       call lin_rad_copy(src%lin, dest%lin)
       dest%kind = lin_t
    elseif(src%kind == log_t) then
       allocate(dest%log)
       call log_rad_copy(src%log, dest%log)
       dest%kind = log_t
    else
       call die("radial: rad_copy unknown type!")
    endif
  end subroutine rad_copy

  !--------------------------------------------------------------------------

  function rad_cutoff(rad_func) result (cutoff)
    !Obtain the cutoff radius of a radial function
    type (rad_func_t), intent(in) :: rad_func
    real(dp)                      :: cutoff

    cutoff = 0.0_dp
    if (rad_func%kind == lin_t) then
       cutoff = lin_rad_cutoff(rad_func%lin)
    elseif(rad_func%kind == log_t) then
       cutoff = log_rad_cutoff(rad_func%log)
    else
       call die("radial: rad_cutoff unknown type!")
    endif
  end function rad_cutoff

  !-----------------------------------------------------------------------

  subroutine rad_dealloc(rad_func)
    !Deallocate a radial function
    type (rad_func_t), intent(inout) :: rad_func

    if (rad_func%kind == lin_t) then
       call lin_rad_dealloc(rad_func%lin)
       deallocate(rad_func%lin)
    elseif(rad_func%kind == log_t) then
       call log_rad_dealloc(rad_func%log)
       deallocate(rad_func%log)
    else
       call die("radial: rad_dealloc unknow function kind")
    endif

  end subroutine rad_dealloc

  !--------------------------------------------------------------------------

  subroutine rad_dealloc_kind(rad_func,kind)
    !Deallocate a kind data structure.
    type (rad_func_t), intent(inout) :: rad_func
    type (rad_kind_t), intent(in)    :: kind

    if (kind == lin_t) then
       call lin_rad_dealloc(rad_func%lin)
       rad_func%kind = log_t
       deallocate(rad_func%lin)
    elseif(kind == log_t) then
       call log_rad_dealloc(rad_func%log)
       deallocate(rad_func%log)
       rad_func%kind = lin_t
    else
       call die("radial: rad_dalloc_kind unknown type!")
    endif

  end subroutine rad_dealloc_kind

  !----------------------------------------------------------------------

  function rad_default_length(func) result(length)
    !Obtain the default number of points of a function (grid/log).
    type (rad_func_t), intent(in) :: func

    integer :: length

    length=-1

    if(func%kind == lin_t) then
       length = lin_rad_default_length()
    elseif(func%kind == log_t) then
       length = log_rad_get_default_length()
    else
       call die("radial: defualt length: unknown func type!")
    endif
  end function rad_default_length

  !--------------------------------------------------------------------------

  function rad_divide_by_4pir2(rad_func,update) result(divided)
    !Divide a function by 4pir2
    type (rad_func_t), intent(in)   :: rad_func
    logical, intent(in)             :: update

    type (rad_func_t) :: divided

    if (rad_func%kind == lin_t) then
       !call rad_divide_by_4pir2(rad_func%lin)
       call die("radial: rad_divide_by_4pir2 not implemented for lin funcs!")
    elseif(rad_func%kind == log_t) then
       allocate(divided%log)
       divided%kind = log_t
       divided%log = log_rad_divide_by_4pir2(rad_func%log,update)
    else
       call die("radial: rad_divide_by_4pir2 unknown type!")
    endif
  end function rad_divide_by_4pir2

  !------------------------------------------------------------------------

  subroutine rad_dump_ascii(rad_func,io, header)
    !Dump the function into a file
    type (rad_func_t), intent(in) :: rad_func
    integer, intent(in) :: io
    logical,intent(in),optional :: header

    if (rad_func%kind == lin_t) then
       if(present (header)) then
          call lin_rad_dump_ascii(rad_func%lin,io,header)
       else
          call lin_rad_dump_ascii(rad_func%lin,io)
       endif
    elseif(rad_func%kind == log_t) then
       call log_rad_dump_ascii(rad_func%log, io)         
    else
       call die("radial: rad_dump_ascii unknown type!")
    endif
  end subroutine rad_dump_ascii

  !--------------------------------------------------------------------------

  subroutine rad_dump_xml(rad_func,io)
    !Dump a function into a xml file
    type (rad_func_t), intent(in) :: rad_func
    integer, intent(in) :: io

    if (rad_func%kind == lin_t) then
       call lin_rad_dump_xml(rad_func%lin,io) 
    elseif(rad_func%kind == log_t) then
       !call log_rad_dump_xml(rad_func%log) 
       call die("radial: rad_dump_xml for log_t not implemented")
    else
       call die("radial: rad_dump_xml unknown type!")
    endif
  end subroutine rad_dump_xml

  !-------------------------------------------------------------------------

  subroutine rad_dump_fft_xml(rad_func,io)
    !Dump the fft of a function into a xml file
    type (rad_func_t), intent(out) :: rad_func
    integer, intent(in) :: io

    if (rad_func%kind == lin_t) then
       call lin_rad_dump_fft_xml(rad_func%lin,io) 
    elseif(rad_func%kind == log_t) then
       !call log_rad_dump_fft_xml(rad_func%log) 
       call die("radial: rad_dump_fft_xml for log_t not implemented")
    else
       call die("radial: rad_dump_fft_xml unknown type!")
    endif
  end subroutine rad_dump_fft_xml

  !--------------------------------------------------------------------------

  subroutine rad_dump_fft_ascii(rad_func,io)
    !Dump the fft of a function into an ascii files
    type (rad_func_t), intent(in) :: rad_func
    integer, intent(in) :: io

    if (rad_func%kind == lin_t) then
       call lin_rad_dump_fft_ascii(rad_func%lin,io) 
    elseif(rad_func%kind == log_t) then
       !call log_rad_dump_fft_ascii(rad_func%log) 
       call die("radial: rad_dump_fft_ascii for log_t not implemented")
    else
       call die("radial: rad_dump_fft_ascii unknown type!")
    endif
  end subroutine rad_dump_fft_ascii

  !-------------------------------------------------------------------------

   subroutine rad_dump_fft_file(rad_func,label)
     !Dump the fft of a function into a file
    type (rad_func_t), intent(in) :: rad_func
    character(len=*), intent(in)   :: label

    if (rad_func%kind == lin_t) then
       call lin_rad_dump_fft_file(rad_func%lin,label)
    elseif(rad_func%kind == log_t) then
       !call log_rad_dump_fft_file(rad_func%log,label) 
       call die("radial: rad_dump_fft_file for log not implemented")
    else
       call die("radial: rad_dump_fft_file unknown type!")
    endif
  end subroutine rad_dump_fft_file

  !-------------------------------------------------------------------------

  subroutine rad_dump_file(rad_func,label)
    !Dump a function to a file
    type (rad_func_t), intent(in) :: rad_func
    character(len=*), intent(in)   :: label

    if (rad_func%kind == lin_t) then
       call lin_rad_dump_file(rad_func%lin,label)
    elseif(rad_func%kind == log_t) then
       call log_rad_dump_file(rad_func%log,label) 
    else
       call die("radial: rad_dump_file unknown type!")
    endif
  end subroutine rad_dump_file

  !-------------------------------------------------------------------------
  subroutine rad_dump_funcs_ascii(rad_func1,rad_func2,io)
    !Dump two functions in ascii
    type(rad_func_t), intent(in) :: rad_func1, rad_func2
    integer, intent(in)          :: io

    if (rad_func1%kind == lin_t .and. rad_func2%kind == lin_t) then
       call die("radial: rad_dump_funcs_ascii: not implemented for lin funcs.")
    elseif(rad_func1%kind == log_t .and. rad_func2%kind == log_t) then
       call log_rad_dump_funcs_ascii(rad_func1%log,rad_func2%log,io)
    else
       call die("radial: rad_dump_funcs_ascii: unknown kind.")
    endif
  end subroutine rad_dump_funcs_ascii

  !--------------------------------------------------------------------------

  function rad_energy_deriv(rad_func,vps,ve,l,e) result(ederiv)
    !To be moved
    type(rad_func_t) , intent(in) :: rad_func
    type(rad_func_t) , intent(in) :: vps, ve
    integer, intent(in)           :: l
    real(dp), intent(in)       :: e

    type(rad_func_t) :: ederiv

    if (rad_func%kind == lin_t) then
       call die("radial: rad_energy_deriv lin_t not implemented")        
    elseif(rad_func%kind == log_t ) then
       allocate(ederiv%log)
       ederiv%kind = log_t
       ederiv%log = log_rad_energy_deriv(rad_func%log,vps%log,ve%log,l,e)
    else
       call die("radial: rad_energy_deriv unknown type!")
    endif
  endfunction rad_energy_deriv

  !---------------------------------------------------------------------------

  subroutine rad_fft(func,l)
    !Find fft of a radial function
    type(rad_func_t), intent(inout) :: func
    integer, intent(in) :: l

    if(func%kind == lin_t) then
       call lin_rad_fft(func%lin,l)
    elseif(func%kind == log_t) then
       call die("radial: log_fft not implemented")
    else
       call die("radial: fft unknown func kind")
    endif
    
  endsubroutine rad_fft
  
  !---------------------------------------------------------------------------

  function rad_get_filter_cutoff(rad_func,l,etol) result(kc)
    !Given a tolerance intent the kinetic energy this function
    !returns the corresponding reciprocal space cutoff.
    !See module filter.f90

    type(rad_func_t), intent(in) :: rad_func
    integer, intent(in) :: l
    real(dp), intent(in) :: etol 
    real(dp) :: kc 

    type(rad_func_t) :: rad_tmp
    
    kc = 0.0_dp 

    if (rad_func%kind == lin_t) then
       kc = lin_rad_get_filter_cutoff(rad_func%lin,l,etol)
    elseif(rad_func%kind == log_t) then
       call rad_copy(rad_func,rad_tmp)
       call rad_log_to_linear(rad_tmp)
       kc = lin_rad_get_filter_cutoff(rad_tmp%lin,l,etol)
       call rad_dealloc(rad_tmp)
    else
       call die("radial: rad_get_filter_cutoff unknown type!")
    endif

  end function rad_get_filter_cutoff

  !---------------------------------------------------------------------------

  function rad_filter(rad_func,l,factor,norm_opt,kc) result(filtered)
    !Remove the components of a function in reciprocal space.
    ! The threshold is defined in kc. 
    ! The factor scales the kc 
    ! The norm_opt specifies the normalization option of the filtered function
    ! See notes in module filter.f90
    type(rad_func_t), intent(in) :: rad_func
    integer, intent(in) :: l,norm_opt
    real(dp), intent(in) :: factor
    real(dp), intent(in) :: kc !Filter cutoff in reciprocal space

    type(rad_func_t) :: filtered,rad_tmp

    allocate(filtered%lin)
    filtered%kind=lin_t

    if (rad_func%kind == lin_t) then
       filtered%lin = lin_rad_filter(rad_func%lin,l,factor,norm_opt,kc) 
       call rad_fft(filtered,l)
    elseif(rad_func%kind == log_t) then
       !Convert to lin, filter and then convert to log the filtered func.
       call rad_copy(rad_func,rad_tmp)
       call rad_log_to_linear(rad_tmp)
       !call rad_fft(rad_tmp,l)
       filtered%lin = lin_rad_filter(rad_tmp%lin,l,factor,norm_opt,kc)
       !call rad_fft(filtered,l)
       !call rad_dump_file(filtered,"test-filt.dat")
       call rad_dealloc(rad_tmp)
       call rad_lin_to_log(filtered,rad_func%log%grid)
       !call rad_dump_file(filtered,"test-filt-log.dat")
    else
       call die("radial: rad_filter unknown type!")
    endif

  end function rad_filter

  !------------------------------------------------

  subroutine rad_find_parabola_parameters(func,l,spln,rmatch,const1,const2,dpbug)
    !Obtain the parameters defining a parabola
    !To be moved?
    type(rad_func_t), intent(in) :: func
    real(dp), intent(in) :: spln
    integer, intent(in) :: l
    real(dp), intent(out) :: const1,const2,rmatch
    logical,intent(in) :: dpbug

    if(func%kind == lin_t) then
       call die("rad_parabola not implemented for lin")
    elseif(func%kind == log_t) then
       call log_rad_parab_params(func%log,l,spln,rmatch,const1,const2,dpbug)
    else
       call die("rad_parabola unknown function kind")
    end if
  end subroutine rad_find_parabola_parameters

  !------------------------------------------------------------------------


  subroutine rad_fit_parabola(func,rmatch,l,const1,const2,parab_norm)
    !Fit a parabola
    type(rad_func_t), intent(in) :: func
    real(dp), intent(in) :: rmatch
    integer, intent(in) :: l
    real(dp), intent(out) :: const1,const2
    real(dp), intent(out) :: parab_norm

    if(func%kind == lin_t) then
       call die("rad_parabola not implemented for lin")
    elseif(func%kind == log_t) then
       call log_rad_fit_parabola(func%log,rmatch,l,const1,const2,parab_norm)
    else
       call die("rad_parabola unknown function kind")
    end if
  end subroutine rad_fit_parabola

  !------------------------------------------------------------------------

  subroutine rad_get(rad_func,r,fr,dfdr)
    !Obtain the the value of a function at a given r
    type (rad_func_t), intent(in) :: rad_func
    real(dp), intent(in)         :: r
    real(dp), intent(out)        :: fr
    real(dp), intent(out)        :: dfdr

    if (rad_func%kind == lin_t) then
       call lin_rad_get(rad_func%lin,r,fr,dfdr)
    elseif(rad_func%kind == log_t) then
       call log_rad_get(rad_func%log,r,fr,dfdr)
    else
       call die("radial: rad_get unknown type!")
    endif
  end subroutine rad_get

  !---------------------------------------------------------

  function rad_get_grid(rad_func) result(grid)
    !Obtain the grid
    type(rad_func_t), intent(in) :: rad_func
    type(rad_Grid_t) :: grid
    if(rad_func%kind == lin_t) then
       allocate(grid%lin)
       grid%lin = lin_rad_get_grid(rad_func%lin)
       grid%kind = lin_t
    elseif(rad_func%kind == log_t)then
       allocate(grid%log)
       grid%log = rad_func%log%grid
       grid%kind = log_t
    else
       call die("radial: rad_get_grid: unknown grid kind")
    endif
  end function rad_get_grid

  !-----------------------------------------------------------------------

  subroutine rad_grid_dump_ascii_formatted(grid,io_ps)
    !save a function in ascii formatted.
    type(rad_grid_t), intent(in) :: grid
    integer, intent(in) :: io_ps

    if(grid%kind==lin_t) then
       call die("radial: rad_grid_dump_ascii is not implemented for lin funcs")
    elseif(grid%kind==log_t)then
       call log_grid_dump_ascii_formatted(grid%log,io_ps)
    else
       call die("radial: rad_grid_get_dump_ascii unknown function kind!")
    endif
  end subroutine rad_grid_dump_ascii_formatted

  !-----------------------------------------------------------------------

  function rad_grid_get_a(grid) result(a)
    !Obtain the a from the grid definition: r(i)= b*[ exp( a*(i-1) ) - 1 ]
    type(rad_grid_t), intent(in) :: grid
    real(dp) :: a

    a = 0.0_dp
    if(grid%kind==lin_t) then
       call die("radial: rad_grid_get_a doesn't make sense for lin funcs")
    elseif(grid%kind==log_t)then
       a = log_grid_get_a(grid%log)
    else
       call die("radial: rad_grid_get_a unknown function kind!")
    endif
       
  end function rad_grid_get_a

  !-----------------------------------------------------------------------

  function rad_grid_get_b(grid) result(b)
    !Obtain the b from the grid definition: r(i)= b*[ exp( a*(i-1) ) - 1 ]
    type(rad_grid_t), intent(in) :: grid
    real(dp) :: b

    b = 0.0_dp
    if(grid%kind==lin_t) then
       call die("radial: rad_grid_get_b doesn't make sense for lin funcs")
    elseif(grid%kind==log_t)then
       b = log_grid_get_b(grid%log)
    else
       call die("radial: rad_grid_get_b unknown function kind!")
    endif
       
  end function rad_grid_get_b

  !-----------------------------------------------------------------------

  function rad_grid_get_length(grid) result(length)
    !Obtain the number of points in the grid
    !To be deprecated
    type(rad_grid_t), intent(in) :: grid
    integer :: length

    length = 0
    if(grid%kind==lin_t)then
       call die("radial: rad_grid_get_length not implemented for lin funcs")
    elseif(grid%kind==log_t) then
       length=log_grid_get_number_of_points(grid%log)
    else
       call die("radial: rad_grid_get_length unknown kind of grid!")
    endif
  end function rad_grid_get_length

  !-----------------------------------------------------------------------

  function rad_get_r_from_ir(rad_func,ir) result (r)
    !Obtain the r from the grid index ir.
    !To be deprecated.
    type (rad_func_t), intent(in)  :: rad_func
    integer, intent(in)            :: ir
    real(dp)                       :: r

    r = 0.0_dp

    if (rad_func%kind == lin_t) then
       !r = lin_rad_get_r_from_ir(rad_func%lin,ir) 
       call die("radial: rad_get_r_from_ir lin_t not implemented")
    elseif(rad_func%kind == log_t) then
       r = log_rad_get_r(rad_func%log,ir)
    else
       call die("radial: rad_get_r_from_ir unknown type!")
    endif
  end function rad_get_r_from_ir

  !---------------------------------------------------------------------------------

  function rad_get_r_from_value(rad_func,value) result (r)
    !Obtain the r corresponding to a given value of a radial function
    type (rad_func_t), intent(in)  :: rad_func
    real(dp), intent(in)            :: value
    real(dp)                       :: r

    r = 0.0_dp
    if (rad_func%kind == lin_t) then
       !r = lin_rad_get_r_from_ir(rad_func%lin,ir) 
       call die("radial: rad_get_r_from_ir lin_t not implemented")
    elseif(rad_func%kind == log_t) then
       r = log_rad_get_r_from_value(rad_func%log,value)
    else
       call die("radial: rad_get_r_from_ir unknown type!")
    endif
  end function rad_get_r_from_value

  !---------------------------------------------------------------------------------

  function rad_get_value_from_ir(rad_func,ir) result (value)
    !Obtain the value of a function from the array index
    !To be deprecated.
    type (rad_func_t), intent(in) :: rad_func
    integer,  intent(in)         :: ir
    real(dp)                      :: value

    value = 0.0_dp

    if (rad_func%kind == lin_t) then
       value = lin_rad_get_value_from_ir(rad_func%lin,ir)
    elseif(rad_func%kind == log_t) then
       value = log_rad_get_value_from_ir(rad_func%log,ir)
    else
       call die("radial: rad_get_value_from_ir unknown type!")
    endif
  end function rad_get_value_from_ir

  !----------------------------------------------------------------------------

  subroutine rad_ghost(func,vps,vlocal,ve,l,e,zval,ighost) 
    !Check if there are "ghost" states
    !To be moved to  kb.f90 
    type(rad_func_t) , intent(in) :: func, vps, vlocal, ve
    integer, intent(in)           :: l
    real(dp), intent(in)          :: e, zval
    integer, intent(out)          :: ighost

    if (func%kind == lin_t) then
       call die("radial: rad_ghost lin_t not implemented")        
    elseif(func%kind == log_t ) then
       call log_rad_ghost(func%log,vps%log,vlocal%log,ve%log,l,e,zval,ighost)
    else
       call die("radial: rad_ghost unknown type!")
    endif
  end subroutine rad_ghost

  !-----------------------------------------------------------------

  function rad_grid_alloc(length,delta,a,b,r) result(grid)
    !Allocated a grid (linear or logarithmic)
    integer, intent(in), optional :: length
    real(dp), intent(in), optional :: delta
    real(dp), intent(in), optional :: a,b 
    real(dp), dimension(:), intent(in), optional :: r

    type(rad_grid_t) :: grid

    if (present(delta) .and. .not. present(a) .and. .not. present(b) .and. &
         .not. present(r) )then
       allocate(grid%lin)
       call lin_grid_alloc(grid%lin,length,delta)
       grid%kind = lin_t
    elseif(.not. present(delta) .and. present(a) .and. present(b) .and. present(r))then
       allocate(grid%log)
       grid%log = log_grid_alloc(length,a,b,r)
       grid%kind = log_t
    else
       call die("")
    endif
  end function rad_grid_alloc

  !--------------------------------------------------------------------------

  subroutine rad_grid_dealloc(grid)
    !Deallocate a grid
    type(rad_grid_t), intent(inout) :: grid

    if(grid%kind== lin_t) then
       !call lin_rad_grid_dealloc(grid%lin)
       deallocate(grid%lin)
    elseif(grid%kind == log_t) then
       call log_grid_dealloc(grid%log)
       deallocate(grid%log)
    else
       call die("radial: rad_grid_dealloc unknown kind of grid!")
    end if
  end subroutine rad_grid_dealloc

  !--------------------------------------------------------------------------
  
  function rad_is_log(func) result(is)
    !Check if a given function is defined in a logarithmic grid.
    type(rad_func_t), intent(in) :: func
    logical :: is

    if (func%kind == log_t) then
       is = .true.
    else
       is = .false.
    endif
  end function rad_is_log

  !--------------------------------------------------------------------------

  function rad_is_lin(func) result(is)
    !Check if a given function is defined in a linear grid.
    type(rad_func_t), intent(in) :: func
    logical :: is

    if (func%kind == lin_t) then
       is = .true.
    else
       is = .false.
    endif
  end function rad_is_lin

  !---------------------------------------------------------------------------
  function rad_kbproj(func,pseudo,vlocal,l,dkbcos,ekb) result(kb_proj)
    !Calculate the kb projectors
    !To be moved to kb.f90
    type(rad_func_t), intent(in) :: func, pseudo,vlocal
    integer, intent(in) :: l
    real(dp), intent(out) :: dkbcos,ekb

    type(rad_func_t) :: kb_proj

    if(func%kind == lin_t) then
       call die("rad_kbproj not implemented for lin functions")
    elseif(func%kind == log_t) then
       allocate(kb_proj%log)
       kb_proj%kind = log_t
       kb_proj%log = log_rad_kb_proj(func%log,pseudo%log,vlocal%log,l,dkbcos,ekb)
    else
       call die("rad_kbproj unknown rad func kind")
    endif
  end function rad_kbproj

  !------------------------------------------------------------------

  function rad_kind(rad_func) result(kind)
    !Query the kind of radial function, it can be either logarithmic or linear.
    type(rad_func_t), intent(in) :: rad_func

    type(rad_kind_t) :: kind

    if (rad_func%kind == log_t) then
       kind = log_t
    elseif(rad_func%kind == lin_t) then
       kind = lin_t
    else
       call die("rad_kind unknown rad func kind")
    endif
  end function rad_kind

  !------------------------------------------------------------------

  subroutine rad_kind_broadcast(func)
    !Broad cast (MPI) the kind of function

    use parallel, only : Node,  nodes

#ifdef MPI
    use mpi_siesta
#endif

    type (rad_func_t), intent(inout)   :: func
    

#ifndef MPI
  end subroutine rad_kind_broadcast
#else

  integer :: mpiError
  if (Nodes .eq. 1) return


  call MPI_Bcast(func%kind%kind,1,MPI_Integer,0,MPI_Comm_World,MPIerror)

  end subroutine rad_kind_broadcast

#endif

  !------------------------------------------------------------------

  function rad_kind_compare(k1,k2) result (compare)
    !Check if two radial kinds are the same.
    type(rad_kind_t), intent(in) :: k1,k2
    logical                  :: compare 
    compare = .false.
    if(k1%kind == k2%kind) compare = .true.
  end function rad_kind_compare

  !------------------------------------------------------------------

  function rad_min(func) result(min)
    !Obtain the minimum value of a function 
    type(rad_func_t), intent (in) :: func
    real(dp) :: min

    min  = 0.0_dp

    if (func%kind == lin_t) then
       call die("rad_min not implemented for lin_t")
    elseif(func%kind == log_t) then
       min = log_rad_min(func%log)
    else
       call die("rad_min: unknown function type")
    endif
  end function rad_min

  !-----------------------------------------------------------------------

  function rad_parabolic_split(func,rmatch,const1,const2,l,nsm,shells,lambdas) result(parabola)
    !Split function.
    !To be moved to multiple-zeta generation function.
    type(rad_func_t), intent(in) :: func
    real(dp), intent(in)  :: rmatch
    real(dp), intent(in) :: const1,const2
    integer, intent(in) :: l,nsm
    type(rad_func_t), intent(in) :: shells(:)
    real(dp), intent(in) :: lambdas(:)

    type(rad_func_t)             :: parabola

    type(log_rad_func_t), pointer :: log_shells(:)

    integer :: i

    if(func%kind == lin_t) then
       call die("rad_parabola not implemented for lin")
    elseif(func%kind == log_t) then
       allocate(log_shells(1:nsm))
       do i=1,nsm
          call log_rad_copy(shells(i)%log,log_shells(i))
       end do
       allocate(parabola%log)
       parabola%kind = log_t
       parabola%log = log_rad_parabolic_split(func%log,rmatch,const1,const2,l,nsm,log_shells,lambdas)

       do i=1,nsm
          call log_rad_dealloc(log_shells(i))
       enddo
       deallocate(log_shells)

    else
       call die("rad_parabolic unknown function kind")
    end if
  end function rad_parabolic_split


  !------------------------------------------------------------------------

  subroutine rad_read_ascii(rad_func,io)
    !Read a function in ascii
    type (rad_func_t), intent(out) :: rad_func
    integer, intent(in) :: io

    if (rad_func%kind == lin_t) then
       call lin_rad_read_ascii(rad_func%lin,io) 
    elseif(rad_func%kind == log_t) then
       !call log_rad_ascii(rad_func%log) 
       call die("radial: rad_read_ascii for log_t not implemented")
    else
       call die("radial: rad_read_ascii unknown type!")
    endif
  end subroutine rad_read_ascii

  !---------------------------------------------------------------------------------
  function rad_reparametrize(rad_func,a,b) result(new_func)
    !Re-parametrize a function 
    type (rad_func_t), intent(in) :: rad_func
    real(dp), intent(in) :: a,b

    type(rad_func_t) :: new_func

    if (rad_func%kind == lin_t) then
       call die("radial: reparametrizing the log grid of a lin rad_func")
    elseif(rad_func%kind == log_t) then
       allocate(new_func%log)
       new_func%kind = log_t
       new_func%log = log_rad_reparametrize(rad_func%log,rmax_radial_grid,a,b)     
    else
       call die("radial: rad_reparam unknown type!")
    endif
  end function rad_reparametrize

  !--------------------------------------------------------------------------
  subroutine rad_read_ascii_unformatted(rad_func,io)
    !Read a function in ascii and unformatted
    type (rad_func_t), intent(out) :: rad_func
    integer, intent(in) :: io

    if (rad_func%kind == lin_t) then
       !call lin_rad_read_ascii_unformatted(rad_func%lin) 
       call die("radial: rad_read_ascii_unformatted for lin_t not implemented")
    elseif(rad_func%kind == log_t) then
       call log_rad_read_ascii_unformatted(io,rad_func%log)         
    else
       call die("radial: rad_read_ascii unknown type!")
    endif
  end subroutine rad_read_ascii_unformatted

  !---------------------------------------------------------------------------
  subroutine rad_read_ascii_formatted(rad_func,io)
     !Read a function in ascii and formatted
    type (rad_func_t), intent(inout) :: rad_func
    integer, intent(in) :: io

    if (rad_func%kind == lin_t) then
       !call lin_rad_read_ascii_formatted(rad_func%lin) 
       call die("radial: rad_read_ascii_formatted for lin_t not implemented")
    elseif(rad_func%kind == log_t) then
       call log_rad_read_ascii_formatted(io,rad_func%log)         
    else
       call die("radial: rad_read_ascii unknown type!")
    endif
  end subroutine rad_read_ascii_formatted

  !---------------------------------------------------------------------------

  subroutine rad_write_ascii_formatted(rad_func,io)
    !Write a function in ascii and unformatted
    type (rad_func_t), intent(in) :: rad_func
    integer, intent(in) :: io

    if (rad_func%kind == lin_t) then
       !call lin_rad_write_ascii_formatted(rad_func%lin) 
       call die("radial: rad_read_ascii_formatted for lin_t not implemented")
    elseif(rad_func%kind == log_t) then
       call log_rad_write_ascii_formatted(io,rad_func%log)         
    else
       call die("radial: rad_read_ascii unknown type!")
    endif
  end subroutine rad_write_ascii_formatted

  !----------------------------------------------------------------------

  subroutine rad_set_origin(rad_func,value)
    !Set the value of a function at the origin.
    type (rad_func_t), intent(inout) :: rad_func
    real(dp), intent(in) :: value

    if (rad_func%kind == lin_t) then
       call lin_rad_set_origin(rad_func%lin,value)         
    elseif(rad_func%kind == log_t) then
       call log_rad_set_origin(rad_func%log,value)         
    else
       call die("radial: rad_set_origin(rad_func,r) unknown type!")
    endif
  end subroutine rad_set_origin

  !--------------------------------------------------------------------------

  subroutine rad_set_default_length(grid,length)
    !Establish the default length of a grid. 
    type(rad_grid_t), intent(in) :: grid
    integer, intent(in) :: length
    if(grid%kind == lin_t) then
       call die("rad_set_default_length not implemented for lin funcs")
    elseif(grid%kind == log_t) then
       call log_rad_set_default_length(length)
    else
       call die("radial: rad_set_default_length unknown type!")
    endif
  end subroutine rad_set_default_length

  !---------------------------------------------------------------------------
 
  subroutine rad_lin_to_log(func,grid)
    !Convert a function from a linear to logarithmic grid.
    !The a, b parameters of the grid should have been specified beforehand.
    type(rad_func_t), intent(inout)         :: func
    type(logGrid_t), intent(in)                 :: grid
    
    real(dp), pointer :: logValues(:)
    integer :: ir, size
    real(dp) :: r, val, dv,rc

    rc = rad_cutoff(func)
    size = log_grid_get_ir(grid,rc)
    allocate(logValues(1:size))
    do ir=1, size
       r = grid%r(ir)
       call rad_get(func,r,val,dv)
       logvalues(ir)=val
    enddo
    logvalues(size)=0.0_dp

    allocate(func%log)
    call log_rad_alloc(func%log,logvalues,grid)

    call rad_dealloc_kind(func, lin_t)
    func%kind = log_t
    deallocate(logvalues)

  end subroutine rad_lin_to_log

 !---------------------------------------------------------------------------

  subroutine rad_log_to_linear(func,yp1)
    !     Convert a logarithmic function into a linear one.

    type(rad_func_t), intent(inout)         :: func
    real(dp),          intent(in),optional  :: yp1

    !     Internal variables
    !     
    integer ::  ntbmax, ir
    real(dp) :: r,rc, delta, v,dv
    real(dp), pointer   :: linvalues(:)
    real(dp), parameter :: deltmax=0.05d0
   
    !Check that default sizes are big enough
    rc = rad_cutoff(func)
    ntbmax = lin_rad_default_length()
    delta = rc /dble(ntbmax-1)  

    if(delta .gt.deltmax) then 
       write(6,'(a)') 'radial: WARNING It might be a good idea to increase' 
       write(6,'(a)') 'radial: WARNING parameter ntbmax (in file atmparams.f) '
       write(6,'(a,i6)') 'radial: WARNING to at least ntable = ', nint(rc/deltmax)+2 
    endif

    allocate(linvalues(1:ntbmax))
    allocate(func%lin)
    linvalues(:)=0.0_dp
    r=0.0_dp
    do ir=1,ntbmax-1
       call rad_get(func,r,v,dv)
       linvalues(ir)=v
       r=r+delta
    enddo
    linvalues(ntbmax)=0.0_dp

    !if (present(f1)) linvalues(1) = f1

    if (present(yp1))then
       call lin_rad_alloc(func%lin,delta,linvalues,yp1)
    else
       call lin_rad_alloc(func%lin,delta,linvalues)
    endif

    call rad_dealloc_kind(func, log_t)
    func%kind = lin_t
    deallocate(linvalues)
  end subroutine rad_log_to_linear

  !--------------------------------------------------------------------------

  subroutine rad_normalize_r_l_1(rad_func,l) 
    !Normalize the function, taking into account the r*(l+1) factor
    !To be deprecated.

    type (rad_func_t), intent(inout) :: rad_func
    integer, intent(in)            :: l

    if (rad_func%kind == lin_t) then
       call die("radial: rad_normalize_r_l_1 lin_t not implemented")
    elseif(rad_func%kind == log_t) then
       call log_rad_normalize_r_l_1(rad_func%log,l)
    else
       call die("radial: rad_normalize_r_l_1 unknown type!")
    endif
  end subroutine rad_normalize_r_l_1

  !-------------------------------------------------------------------

  function rad_divide_by_r_l_1(rad_func,l,lambda) result(div)
    !Divide the function by r*(l+1) (matel spects all the functions
    ! divided by r**l, the +1 comes from the solution of radial part 
    ! of the Schrodinger Equation).
    type (rad_func_t), intent(in) :: rad_func
    integer, intent(in)            :: l
    real(dp), intent(in)           :: lambda

    type(rad_func_t) :: div

    if (rad_func%kind == lin_t) then
       call die("radial: rad_divide_r_l_1 lin_t not implemented")
    elseif(rad_func%kind == log_t) then
       allocate(div%log)
       div%kind = log_t
       div%log = log_rad_divide_by_r_l_1(rad_func%log,l,lambda)
    else
       call die("radial: rad_normalize_r_l_1 unknown type!")
    endif
  end function rad_divide_by_r_l_1

  !-------------------------------------------------------------------

  function rad_kinetic_energy(rad_func,l) result(ekin)
    ! Obtain the kinetic energy of a radial function.
    type (rad_func_t), intent(in) :: rad_func
    integer, intent(in)            :: l

    real(dp) :: ekin

    ekin = 0.0_dp

    if (rad_func%kind == lin_t) then
       call die("radial: rad_kinetic_energy lin_t not implemented")
    elseif(rad_func%kind == log_t) then
       ekin = log_rad_kinetic_energy(rad_func%log,l)
    else
       call die("radial: rad_kinetic_energy unknown type!")
    endif

  end function rad_kinetic_energy

  !-------------------------------------------------------------------

  function rad_potential_energy(rad_func,pot,l) result(epot)
    !Obtain the potential energy 
    !To be moved.
    type (rad_func_t), intent(in) :: rad_func
    type (rad_func_t), intent(in) :: pot
    integer, intent(in)           :: l

    real(dp) :: epot

    if (rad_func%kind == lin_t .and. pot%kind == lin_t ) then
       call die("radial: rad_potential_energy lin_t not implemented")
    elseif(rad_func%kind == log_t .and. pot%kind == log_t ) then
       epot = log_rad_potential_energy(rad_func%log,pot%log,l)
    else
       call die("radial: rad_potential_energy unknown type!")
    endif

  end function rad_potential_energy

  !-------------------------------------------------------------------

  function rad_self_energy(func) result(energy)
    !Self energy, to be moved.
    type (rad_func_t), intent(in) :: func
    real(dp)                      :: energy

    energy = 0.0_dp
    if (func%kind == lin_t) then
       call die("radial: rad_self_energy lin_t not implemented")
    elseif(func%kind == log_t) then
       energy = log_rad_self_energy(func%log)
    else
       call die("radial: rad_self_energy unknown type!")
    endif

  end function rad_self_energy

  !------------------------------------------------------------------------

  function rad_split_gauss(first_z,rc,lambda,l) result(gauss)
    !To be moved
    type(rad_func_t), intent(in)   :: first_z
    real(dp), intent(in)           :: rc, lambda
    integer, intent(in)            :: l

    type(rad_func_t) :: gauss

    if(first_z%kind == lin_t) then
       call die("radial: split gauss not implemented for lin functions")
    elseif(first_z%kind == log_t) then
       allocate(gauss%log)
       gauss%kind = log_t
       gauss%log  = log_rad_split_gauss(first_z%log,rc,lambda,l)
    else
       call die("radial: rad_split_gauss: unknown kind!")
    endif
  end function rad_split_gauss

  !-------------------------------------------------------------------

  function rad_sum(func1,func2,rmax) result(sum)
    !Sum two functions
    type (rad_func_t), intent(in) :: func1, func2
    real(dp), intent(in), optional :: rmax
    type (rad_func_t)             :: sum
     
    if (func1%kind == lin_t .and. func2%kind == lin_t) then
       call die("radial: rad_sum lin_t not implemented")
    elseif(func1%kind == log_t .and. func2%kind == log_t) then
       allocate(sum%log)
       sum%kind = log_t
       if(present(rmax)) then
          sum%log = log_rad_sum(func1%log,func2%log,rmax)
       else
          sum%log = log_rad_sum(func1%log,func2%log)
       endif
    else
       call die("radial: rad_sum unknown type!")
    endif

  end function rad_sum

  !--------------------------------------------------------------------

  function rad_multiply(func1,func2) result(mult)
    !Multiply two functions, at every point of the grid.
    type (rad_func_t), intent(in) :: func1, func2
    type (rad_func_t)             :: mult

    if (func1%kind == lin_t .and. func2%kind == lin_t) then
       call die("radial: rad_sum lin_t not implemented")
    elseif(func1%kind == log_t .and. func2%kind == log_t) then
       allocate(mult%log)
       mult%kind = log_t
       mult%log = log_rad_multiply(func1%log,func2%log)
    else
       call die("radial: rad_sum unknown type!")
    endif

  end function rad_multiply

  !--------------------------------------------------------------------

  function rad_integral(func) result(integral)
    !Integrate a function
    type (rad_func_t), intent(in) :: func
    real(dp)                      :: integral

    integral = 0.0_dp
    if (func%kind == lin_t) then
       call die("radial: rad_sum lin_t not implemented")
    elseif(func%kind == log_t) then
       integral = log_rad_integral(func%log)
    else
       call die("radial: rad_sum unknown type!")
    endif

  end function rad_integral

  !--------------------------------------------------------------------

  function rad_sum_function(func1,func2,r0,rmax) result(sum)
    !Sum two functions, from r0 to rmax.
    type (rad_func_t), intent(in) :: func1
    real(dp), intent(in) :: r0, rmax

    interface
       function func2(x) result(y)
         use precision
         real(dp), intent(in) :: x
         real(dp)             :: y
       end function func2
    end interface

    type (rad_func_t)             :: sum

    if (func1%kind == lin_t) then
       call die("radial: rad_sum lin_t not implemented")
    elseif(func1%kind == log_t) then
       allocate(sum%log)
       sum%kind = log_t
       sum%log = log_rad_sum_function(func1%log,func2,r0,rmax)
    else
       call die("radial: rad_sum unknown type!")
    endif

  end function rad_sum_function

  !--------------------------------------------------------------------

  function rad_schro(pseudo,ve,rc,l,nnodes,nprin,chg,energy) result(integral)
    !Integrate the radial part of the Schrodinger equation
    type(rad_func_t), intent(in) :: pseudo
    type(rad_func_t), intent(in) :: ve     !pseudo of valence electrons.
    real(dp), intent(in) :: rc
    integer, intent(in) :: l
    integer, intent(in) :: nnodes !number of nodes
    integer, intent(in) :: nprin  ! l + number of shells with this n
    real(dp), intent(in) :: chg   !charge
    real(dp), intent(out) :: energy !energy

    type(rad_func_t) :: integral

    if (pseudo%kind == lin_t .and. ve%kind == lin_t) then
       call die("radial: rad_schro lin_t not implemented")        
    elseif(pseudo%kind == log_t .and. ve%kind == log_t) then
       allocate(integral%log)
       integral%kind = log_t
       integral%log = log_rad_schro(pseudo%log,ve%log,rc,l,nnodes,nprin,chg,energy)
    else
       call die("radial: rad_schro unknown function kind!")
    endif
  end function rad_schro


  !----------------------------------------------------------------------------- 

  function rad_multiply_by_rl(rad_func,l) result (mult)
    !Multiply a function by rl
    type(rad_func_t), intent(in) :: rad_func
    integer, intent(in) :: l

    type(rad_func_t) :: mult

    if (rad_func%kind == lin_t) then
       call die("radial: rad_multiply_by_rl lin_t not implemented")        
    elseif(rad_func%kind == log_t ) then
       allocate(mult%log)
       mult%kind = log_t
       mult%log = log_rad_multiply_by_rl(rad_func%log,l)
    else
       call die("radial: rad_multiply_rl unknown type!")
    endif

  end function rad_multiply_by_rl

  !--------------------------------------------------------

  function rad_sum_value(rad_func,value) result (sum)
    !At every point of the grid sum a given value to the function.
    type(rad_func_t), intent(in) :: rad_func
    real(dp), intent(in) :: value

    type(rad_func_t) :: sum
    

    if (rad_func%kind == lin_t) then
       call die("radial: rad_sum_value lin_t not implemented")        
    elseif(rad_func%kind == log_t ) then
    
       allocate(sum%log)
       sum%log = log_rad_sum_value(rad_func%log,value)
       sum%kind = log_t

    else
       call die("radial: rad_sum_value unknown type!")
    endif

  end function rad_sum_value

  !--------------------------------------------------------

  function rad_get_length(rad_func) result(length)
    !Get the number of points where the function is evaluated.
    type(rad_func_t), intent(in) :: rad_func
    integer :: length

    length = 0
    if (rad_func%kind == lin_t) then
       length = lin_rad_get_length()     
    elseif(rad_func%kind == log_t) then
       length = log_rad_get_length(rad_func%log)       
    else
       call die("radial: get_length unknown type!")
    endif
  end function rad_get_length

  !----------------------------------------------------------

  function rad_get_ir_from_r(rad_func,r) result(ir)
    !Get the index of
    type(rad_func_t), intent(in) :: rad_func
    real(dp), intent(in)  :: r

    integer :: ir

    ir = -1
    if (rad_func%kind == lin_t) then
       call die("radial: rad_get_r_from_ir not implemented for lin")
    elseif(rad_func%kind == log_t) then
       ir = log_rad_get_ir(rad_func%log,r)       
    else
       call die("radial: rad_get_ir_from_r unknown type!")
    endif

  end function rad_get_ir_from_r

  !-----------------------------------------------------------------

  function rad_matching_radius(func,ref) result(rmatch)
    ! Obtain the radius where two functions become identical
    ! the threshold is defined by the min_func_val parameter.
  

    !    This routine returns the maximum radius for the
    !    Kleinman-Bylander projectors with a standard choice
    !    of the local potential.
    !    Check also at which radius the asymptotic 2*Zval/r
    !    behaviour is achieved.  
    !    D. Sanchez-Portal, Aug. 1998

    ! This routine studies several functions and returns the matching radious: 
    ! beyond that radious all the functions have the same value.
    ! E. Anglada 2007 (After the work of D. Sanchez-Portal. 1998)


    type(rad_func_t), intent(in) ::func
    type(rad_func_t), intent(in) ::ref
    real(dp) :: rmatch

    !**   Iterate over the possible local potentials**

    integer :: ir,nr
    real(dp) :: rmax,r,v1,vref,diff
    real(dp), parameter :: eps = 1E-4

    rmax = 0.0_dp
    rmatch = 0.0_dp

    rmax = rad_cutoff(func)

    nr = rad_get_length(func)
    r=rmax
    do ir=nr,2,-1
       r = rad_get_r_from_ir(func,ir)
       v1 = rad_get_value_from_ir(func,ir)
       vref = rad_get_value_from_ir(ref,ir)
       diff = abs(abs(v1)-abs(vref))
       if(diff .gt. eps) exit
    enddo

    rmatch=r 

  end function rad_matching_radius

  !---------------------------------------------------------------------------

  function rad_get_radius_from_value(func,value) result(r_match)
    !Obtain the r corresponding to a given function value.
    type(rad_func_t), intent(in) :: func
    real(dp), intent(in)         :: value

    real(dp) :: r_match

    integer :: ir,nr
    real(dp) :: r,v1,diff
    real(dp), parameter :: dr = 0.1
    real(dp), parameter :: eps = 1E-4

    nr = rad_get_length(func)

    do ir=nr,2,-1
       r = rad_get_r_from_ir(func,ir)
       v1 = rad_get_value_from_ir(func,ir)
       diff = abs(v1*r+value)       
       if(diff .gt. eps) exit
    enddo
    r_match=r 

  endfunction rad_get_radius_from_value

  !----------------------------------------------------------------------------

  function rad_integral_vs_value(vps,ve,l,e,rmax) result(integ_vs_val)
    !
    type(rad_func_t) , intent(in) :: vps, ve
    integer, intent(in)           :: l
    real(dp), intent(in)          :: e
    real(dp), intent(in)          :: rmax

    type(rad_func_t) :: integ_vs_val

    if (vps%kind == lin_t) then
       call die("radial: rad_integral_vs_value lin_t not implemented")        
    elseif(vps%kind == log_t ) then
       allocate(integ_vs_val%log)
       integ_vs_val%log = log_rad_integral_vs_value(vps%log,ve%log,l,e,rmax)
       integ_vs_val%kind = log_t
    else
       call die("radial: rad_integral_vs_value unknown type!")
    endif
  end function rad_integral_vs_value

  !-------------------------------------------------------------------

  subroutine rad_smooth(func,zval,rc,smooth,integral) 
    type(rad_func_t), intent(in) :: func

    real(dp), intent(in)          :: zval, rc
    type(rad_func_t), intent(out) :: smooth, integral

    if (func%kind == lin_t) then
       call die("radial: rad_smooth lin_t not implemented")
    elseif(func%kind == log_t) then
       allocate(smooth%log,integral%log)
       call log_rad_smooth(func%log,zval,rc,smooth%log,integral%log)
       smooth%kind = log_t
       integral%kind = log_t
    else
       call die("radial: rad_smooth unknown type!")
    endif

  end subroutine rad_smooth

  !-----------------------------------------------------------------

  subroutine rad_smooth_large(func,zval,rc,smooth,integral) 
    type(rad_func_t) , intent(in) :: func

    real(dp), intent(in)          :: zval, rc
    type(rad_func_t), intent(out) :: smooth, integral

    if(func%kind == lin_t) then
       call die("radial: rad_smooth_large")
    elseif(func%kind == log_t) then
       allocate(smooth%log,integral%log)
       call log_rad_smooth_large(func%log,zval,rc,smooth%log,integral%log)
       smooth%kind = log_t
       integral%kind = log_t
    else
       call die("radial: rad_smooth_large: unknown kind")
    endif

  end subroutine rad_smooth_large

  !--------------------------------------------------------

  function rad_multiply_each_value(func,value) result(multiplied)
    type(rad_func_t), intent(in) :: func
    real(dp)                        :: value

    type(rad_func_t)     :: multiplied

    if (func%kind == lin_t) then
       call die("radial: rad_multiply_each_value lin_t not implemented")      
    elseif(func%kind == log_t ) then
       allocate(multiplied%log)
       multiplied%log = log_rad_multiply_each_value(func%log,value)
       multiplied%kind = log_t
    else
       call die("radial: rad_multiply_each_value unknown type!")
    endif
  end function rad_multiply_each_value

  !--------------------------------------------------------

  function rad_vhartree(func) result(vhartree)
    type(rad_func_t), intent(in) :: func

    type(rad_func_t)     :: vhartree

    if (func%kind == lin_t) then
       call die("radial: rad_vhartree lin_t not implemented")      
    elseif(func%kind == log_t ) then
       allocate(vhartree%log)
       vhartree%log = log_rad_hartree(func%log)
       vhartree%kind = log_t
    else
       call die("radial: rad_vhartree unknown type!")
    endif
  end function rad_vhartree

  !--------------------------------------------------------

  function rad_vxc(func,irel,exc) result(vxc)
    type(rad_func_t), intent(in) :: func
    integer, intent(in) :: irel
    real(dp), intent(out) :: exc !Exchange-correlation energy

    type(rad_func_t)     :: vxc

    if (func%kind == lin_t) then
       call die("radial: rad_vxc lin_t not implemented")      
    elseif(func%kind == log_t ) then
       allocate(vxc%log)
       vxc%log = log_rad_vxc(func%log,irel,exc)
       vxc%kind = log_t
    else
       call die("radial: rad_vxc unknown type!")
    endif
  end function rad_vxc

  !--------------------------------------------------------

  function rad_rc_vs_e(func,pseudo,ve,l,el,nnodes) result(rc)
    type(rad_func_t), intent(in) :: func, pseudo,ve
    real(dp), intent(in) :: el
    integer, intent(in) :: l, nnodes

    real(dp) :: rc

    rc = 0.0_dp

    if(func%kind == lin_t) then
       call die("rad_rc_vs_e not implemented for lin functions")
    elseif(func%kind == log_t) then
       rc = log_rad_rc_vs_e(func%log,pseudo%log,ve%log,l, el, nnodes)
    else
       call die("rad_rc_vs_e: unknown rad func kind")
    endif
  end function rad_rc_vs_e


  !---------------------------------------------------------------------------

  function rad_split_scan_tail(func) result(split)
    type(rad_func_t), intent(in) :: func

    type(rad_func_t) :: split

    if(func%kind == lin_t) then
       call die("rad_split_scan_tail not implemented for lin functions")
    elseif(func%kind == log_t) then        
       allocate(split%log)
       split%log = log_rad_split_scan_tail(func%log)        
       split%kind = log_t
    else
       call die("rad_split_scan_tail unknown rad func kind")
    endif
  end function rad_split_scan_tail

  !--------------------------------------------------------------

  function rad_polarization(func,pseudo,ve,rc,l,energy) result(pol)
    type(rad_func_t), intent(in) :: func
    type(rad_func_t), intent(in) :: pseudo,ve
    real(dp), intent(in) :: rc
    integer, intent(in) :: l
    real(dp), intent(in) :: energy

    type(rad_func_t) :: pol

    if(func%kind == lin_t) then
       call die("rad_polarization not implemented for lin functions")
    elseif(func%kind == log_t) then        
       allocate(pol%log)
       pol%log = log_rad_polarization(func%log,pseudo%log,ve%log,rc,l,energy) 
       pol%kind = log_t
    else
       call die("rad_polarization unknown rad func kind")
    endif
  end function rad_polarization

  !---------------------------------------------------------------------------

  function rad_split_scan_tail_parabola(func,l,fix_split_table,label) result(split)
    type(rad_func_t), intent(in) :: func
    integer, intent(in) :: l
    logical, intent(in) ::  fix_split_table

    character(len=*), intent(in), optional :: label

    type(rad_func_t) :: split

    if(func%kind == lin_t) then
       call die("rad_split_scan_tail_parabola not implemented for lin functions")
    elseif(func%kind == log_t) then
       allocate(split%log)
       if(present(label)) then
          split%log = log_rad_scan_tail_parabola(func%log,l,fix_split_table,label)
       else
          split%log = log_rad_scan_tail_parabola(func%log,l,fix_split_table)
       end if
       split%kind = log_t
    else
       call die("rad_split_scan_tail_parabola unknown rad func kind")
    endif
  end function rad_split_scan_tail_parabola

  !--------------------------------------------------------------------------

  subroutine rad_zero(rad_func)
    type (rad_func_t), intent(inout) :: rad_func

    rad_func%kind = lin_t
    call lin_rad_zero(rad_func%lin) 
    !elseif(rad_func%kind == log_t) then
    !   call log_rad_zero(rad_func%log) 
    !else
    !   call die("radial: rad_zero unknown type!")
    !endif
  end subroutine rad_zero

  !---------------------------------------------------------------------------

end module radial










