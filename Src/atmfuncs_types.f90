!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

!Definition of the data structure used to identify which kind of radial
!function we are dealing with.

module atmfuncs_types
  use precision
  
  type, public :: func_t
     private
     integer     :: func_type = 0
  end type func_t

  type(func_t), parameter, public ::  &
       orb_f     = func_t(1) , &
       kbpj_f    = func_t(2) , &
       vlocal_f  = func_t(3), &
       vna_f     = func_t(4), &
       chlocal_f = func_t(5), &
       core_f    = func_t(6), &
       ldaupj_f  = func_t(7)

  INTERFACE OPERATOR (==)
     MODULE PROCEDURE func_t_compare
  END INTERFACE

  contains
    
    function get_func_type(code) result(func_type)
      type(func_t), intent(in) :: code

      integer :: func_type
      func_type = code%func_type
    end function get_func_type

  !-------------------------------------------------------------------------
    
    function func_t_compare(f1,f2) result(compare)
      type(func_t), intent(in) :: f1,f2
      logical                  :: compare 
      compare = .false.
      if(f1%func_type == f2%func_type) compare = .true.
    end function func_t_compare
  
  !--------------------------------------------------------------------------

end module atmfuncs_types
