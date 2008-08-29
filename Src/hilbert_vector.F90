module hilbert_vector_m
use radial
use precision
implicit none

type hilbert_vector_t
   private
   type(rad_func_t), pointer     :: rad_func !the radial function
   integer                       :: l,n !quantum numbers
   real(dp)                      :: energy
   real(dp)                      :: pop
   integer                       :: zeta
   logical                       :: pol
end type hilbert_vector_t


contains

subroutine broadcast_hilbert_vector(vector)
  use parallel, only : Node, Nodes

#ifdef MPI
  use mpi_siesta
#endif

  implicit none

  type(hilbert_vector_t), intent(inout) :: vector

  integer :: mpiError

  if (Nodes .eq. 1) return

#ifndef MPI
!
!     Do nothing...
!

#else

  call MPI_Bcast(vector%l,1,MPI_integer,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(vector%n,1,MPI_integer,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(vector%energy,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(vector%pop,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(vector%zeta,1,MPI_integer,0,MPI_Comm_World,MPIerror)

  call MPI_Bcast(vector%pol,1,MPI_logical,0,MPI_Comm_World,MPIerror)
  if (Node .ne. 0) allocate(vector%rad_func)

  call rad_broadcast(vector%rad_func)

#endif
end subroutine broadcast_hilbert_vector


!------------------------------------------------

subroutine init_vector(vector,func,n,l,z,pop,energy,pol)
  type(hilbert_vector_t) , intent (inout) :: vector
  type(rad_func_t), intent(in) :: func 
  integer, intent(in) :: l,n,z
  real(dp), intent(in) :: pop, energy
  logical, intent(in) :: pol

  type(rad_func_t) :: func_lin
  
  allocate(vector%rad_func)
  call rad_copy(func,func_lin)
  if (rad_is_log(func_lin))  call rad_log_to_linear(func_lin)
  !call rad_fft(func_lin,l)
  call rad_copy(func_lin,vector%rad_func)
  vector%l = l
  vector%n = n
  vector%energy = energy
  vector%zeta = z
  vector%pop = pop
  vector%pol = pol
  call rad_dealloc(func_lin)
end subroutine init_vector

!------------------------------------------------

subroutine copy_vector(src,dest)
  type(hilbert_vector_t), intent(in) :: src
  type(hilbert_vector_t), intent(out) :: dest

  call init_vector(dest,src%rad_func,src%n,src%l,src%zeta,src%pop,src%energy,src%pol)
end subroutine copy_vector

!------------------------------------------------

subroutine destroy_vector(vector)
  type(hilbert_vector_t), intent(inout) :: vector
  call rad_dealloc(vector%rad_func)
  deallocate(vector%rad_func)
end subroutine destroy_vector

!------------------------------------------------

function get_l_v(vector) 
  type(hilbert_vector_t) , intent (in) :: vector
  integer :: get_l_v

  get_l_v = vector%l
end function get_l_v

!------------------------------------------------

subroutine set_l_v(vector,l)
    type(hilbert_vector_t) , intent (inout) :: vector
    integer, intent(in) ::  l

    vector%l = l
end subroutine set_l_v

!------------------------------------------------

!function get_m_v(vector) 
!  type(hilbert_vector_t) , intent (in) :: vector
!  integer :: get_m_v
  
!  get_m_v = vector%m
!end function get_m_v

!------------------------------------------------

!subroutine set_m_v(vector,m)
!  type(hilbert_vector_t) , intent (inout) :: vector
!  integer, intent(in) ::  m
  
!  vector%m = m
!end subroutine set_m_v

!------------------------------------------------

function get_n_v(vector) 
  type(hilbert_vector_t) , intent (in) :: vector
  integer :: get_n_v

  get_n_v = vector%n
end function get_n_v

!------------------------------------------------

subroutine set_n_v(vector,n)
  type(hilbert_vector_t) , intent (inout) :: vector
  integer, intent(in) ::  n
  
  vector%n = n
end subroutine set_n_v

!------------------------------------------------

function get_zeta_v(vector) 
  type(hilbert_vector_t) , intent (in) :: vector
  integer :: get_zeta_v

  get_zeta_v = vector%zeta
end function get_zeta_v

!------------------------------------------------

subroutine set_zeta_v(vector,zeta)
  type(hilbert_vector_t) , intent (inout) :: vector
  integer, intent(in) ::  zeta
  
  vector%zeta = zeta
end subroutine set_zeta_v

!------------------------------------------------

function get_rad_func_v(vector) 
  type(hilbert_vector_t), intent (in) :: vector
  type(rad_func_t), pointer:: get_rad_func_v
  get_rad_func_v => vector%rad_func
end function get_rad_func_v

!------------------------------------------------

subroutine set_rad_func_v(vector,func)
  type(hilbert_vector_t) , intent (inout) :: vector
  type(rad_func_t), intent(in) :: func
  allocate(vector%rad_func)
  call rad_copy(func,vector%rad_func) 
end subroutine set_rad_func_v

!-------------------------------------------------

function get_pop_v(vector) 
  type(hilbert_vector_t) , intent (in) :: vector
  real(dp) :: get_pop_v

  get_pop_v = vector%pop
end function get_pop_v

!------------------------------------------------

subroutine set_pop_v(vector,pop)
  type(hilbert_vector_t) , intent (inout) :: vector
  real(dp), intent(in) :: pop
  
  vector%pop = pop
end subroutine set_pop_v

!-------------------------------------------------

subroutine set_energy_v(vector,energy)
  type(hilbert_vector_t) , intent (inout) :: vector
  real(dp), intent(in) :: energy
  
  vector%energy = energy
end subroutine set_energy_v

!-------------------------------------------------


function get_energy_v(vector) 
  type(hilbert_vector_t) , intent (in) :: vector
  real(dp) :: get_energy_v

  get_energy_v = vector%energy
end function get_energy_v

!------------------------------------------------

function get_pol_v(vector) 
  type(hilbert_vector_t) , intent (in) :: vector
  logical :: get_pol_v

  get_pol_v = vector%pol
end function get_pol_v

!------------------------------------------------

subroutine set_pol_v(vector,pol)
  type(hilbert_vector_t) , intent(inout) :: vector
  logical :: pol
  vector%pol = pol
end subroutine set_pol_v

!------------------------------------------------

function get_cutoff_v(vector)
  type(hilbert_vector_t) , intent (in) :: vector
  real(dp) :: get_cutoff_v

  get_cutoff_v = rad_cutoff(vector%rad_func)
end function get_cutoff_v

!-------------------------------------------------

end module  hilbert_vector_m
