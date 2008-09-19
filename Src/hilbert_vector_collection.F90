module hilbert_vector_collection
use precision
use hilbert_vector_m
implicit none

  type hilbert_vector_collection_t
     private
     !Data structure which contains a collection of rad_funcs
     !of the same type (orbitals, kb_proj, lda+u_proj, ...)
     !and the list of degenerate copies (those who share the
     !same l but a different m).
     integer                         ::  lmax ! l max
      
     !Aggregate numbers of orbitals and projectors (including 2l+1
     !        copies for each "nl"), and index arrays keeping track of
     !        which "nl" family they belong to, and their n, l, and m (to avoid
     !        a further dereference)

     integer                         ::  n_funcs = -1! num of funcs including 2l+1 copies.

     integer, pointer, dimension(:)  ::  m => Null() ! m of each nl funcs
     integer, pointer, dimension(:)  ::  index => Null() !Used to identify the nl the function belong.
    
     !The orbitals/projectors themselves (ie the "nl" in the old format):
     type(hilbert_vector_t), dimension(:), pointer :: vec => Null()
     
  end type hilbert_vector_collection_t
  
  
contains

  subroutine allocate_collection(data,length)
    type(hilbert_vector_collection_t), intent(inout) :: data
    integer, intent(in) :: length
    allocate(data%vec(1:length))    
  end subroutine allocate_collection

  !-------------------------------------------------------------------------

  subroutine broadcast_hilbert_vector_collection(data)
    !Distribute the hilbert vector collection among the MPI nodes.
   
    use parallel, only : Node, Nodes

#ifdef MPI
    use mpi_siesta
#endif

    implicit none

    type(hilbert_vector_collection_t), intent(inout) :: data

    integer :: mpiError,n,i

    if (Nodes .eq. 1) return

#ifndef MPI
!
!     Do nothing...
!
    end subroutine broadcast_hilbert_vector_collection
#else

    call MPI_Bcast(data%lmax,1,MPI_integer,0,MPI_Comm_World,MPIerror)

    call MPI_Bcast(data%n_funcs,1,MPI_integer,0,MPI_Comm_World,MPIerror)
   
    if (data%n_funcs .ge. 1) then
       if (Node .ne. 0) allocate(data%m(1:data%n_funcs),data%index(1:data%n_funcs))

       call MPI_Bcast(data%m,data%n_funcs,MPI_integer,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(data%index,data%n_funcs,MPI_integer,0,MPI_Comm_World,MPIerror)

       if (Node .eq. 0) n = size(data%vec)

       call MPI_Bcast(n,1,MPI_integer,0,MPI_Comm_World,MPIerror)

       if (Node .ne. 0) call allocate_collection(data,n)

       do i=1,n
          call broadcast_hilbert_vector(data%vec(i))
       enddo
    endif

  end subroutine broadcast_hilbert_vector_collection

#endif
  !-------------------------------------------------------------------------

  function get_size(data)
    type(hilbert_vector_collection_t), intent(in) :: data
    integer :: get_size
    get_size = data%n_funcs
  end function get_size

  !-----------------------------------------------------------------------
 
  subroutine set_deg(data)
    !Take into account the 2l+1 degeneracy of some vectors (orbs, kbs)
      type(hilbert_vector_collection_t), intent(inout) :: data

      integer :: n_funcs,i,l,m, length

      length = 0
      do i=1,size(data%vec)
         l = get_l_v(data%vec(i))
         length = length + 2*l +1
      end do

      allocate(data%index(1:length))
      allocate(data%m(1:length))

      n_funcs = 0
      do i=1,size(data%vec)
         l = get_l_v(data%vec(i))
         do m = -l,l        
            n_funcs = n_funcs+1
            data%m(n_funcs) = m 
            data%index(n_funcs) = i              
         enddo
      enddo
      data%n_funcs = n_funcs

    end subroutine set_deg

   !--------------------------------------------------------

    function get_lmax(data) 
      type(hilbert_vector_collection_t), intent(in) :: data
      integer get_lmax
      get_lmax = data%lmax
    end function get_lmax

    !--------------------------------------------------------

    subroutine set_lmax(data, l_max) 
      type(hilbert_vector_collection_t), intent(inout) :: data
      integer, intent(in) :: l_max
      data%lmax = l_max
    end subroutine set_lmax

    !--------------------------------------------------------

    function get_n_funcs(data)
      type(hilbert_vector_collection_t), intent(in) :: data
      integer get_n_funcs
      get_n_funcs = data%n_funcs
    end function get_n_funcs

    !--------------------------------------------------------

    subroutine set_n_funcs(data, n)
      type(hilbert_vector_collection_t), intent(inout) :: data
      integer, intent(in) :: n
      data%n_funcs = n
    end subroutine set_n_funcs

    !--------------------------------------------------------

    function get_number_of_vectors(data) result(n)
      type(hilbert_vector_collection_t), intent(in) :: data
      integer n
      n = size(data%vec)
    end function get_number_of_vectors

    !--------------------------------------------------------

    function get_l(data,i) 
      type(hilbert_vector_collection_t), intent(in) :: data
      integer get_l, i
      get_l = get_l_v(data%vec(data%index(i)))
    end function get_l

    !--------------------------------------------------------       

    subroutine set_l(data, i, l)
      type(hilbert_vector_collection_t), intent(inout) :: data
      integer, intent(in) :: i, l 
      call set_l_v(data%vec(i),l)
    end subroutine  set_l

    !--------------------------------------------------------

    function get_n(data,i) 
      type(hilbert_vector_collection_t), intent(in) :: data
      integer get_n, i
      get_n = get_n_v(data%vec(data%index(i)))
    end function get_n

    !--------------------------------------------------------

    subroutine set_n(data, i, n)
      type(hilbert_vector_collection_t), intent(inout) :: data
      integer, intent(in) :: i, n 
      call set_n_v(data%vec(i),n)
    end subroutine  set_n

    !--------------------------------------------------------

    function get_m(data,i) 
      type(hilbert_vector_collection_t), intent(in) :: data
      integer get_m, i
      get_m = data%m(i)
    end function get_m

    !------------------------------------------------------

    !subroutine set_m(data, i, m)
    !  type(hilbert_vector_collection_t), intent(inout) :: data
    !  integer, intent(in) :: i, m 
    !  data%m(i) = m
    !end subroutine set_m

    !--------------------------------------------------------

    function get_zeta(data,i) 
      type(hilbert_vector_collection_t), intent(in) :: data
      integer get_zeta, i
      get_zeta = get_zeta_v(data%vec(data%index(i)))
    end function get_zeta

    !------------------------------------------------------

    subroutine set_zeta(data, i, zeta)
      type(hilbert_vector_collection_t), intent(inout) :: data
      integer, intent(in) :: i, zeta 
      call set_zeta_v(data%vec(i),zeta)
    end subroutine set_zeta

    !--------------------------------------------------------

    function get_energy(data,i) 
      type(hilbert_vector_collection_t), intent(in) :: data
      integer, intent(in) :: i
      real(dp):: get_energy
      get_energy = get_energy_v(data%vec(data%index(i))) 
    end function get_energy

    !------------------------------------------------------

    subroutine set_energy(data, i, energy)
      type(hilbert_vector_collection_t), intent(inout) :: data
      integer, intent(in) :: i 
      real(dp), intent(in) :: energy
       call set_energy_v(data%vec(i), energy) 
    end subroutine set_energy

    !--------------------------------------------------------

    function get_pop(data,i) 
      type(hilbert_vector_collection_t), intent(in) :: data
      integer, intent(in) :: i
      real(dp) :: get_pop
      get_pop = get_pop_v(data%vec(data%index(i)))
    end function get_pop

    !------------------------------------------------------

    subroutine set_pop(data, i, pop)
      type(hilbert_vector_collection_t), intent(inout) :: data
      integer, intent(in)                              :: i 
      real(dp), intent(in)                             :: pop
      call set_pop_v(data%vec(i), pop) 
    end subroutine set_pop

    !-------------------------------------------------------

    function get_pol(data,i) 
      type(hilbert_vector_collection_t), intent(in) :: data
      integer, intent(in) :: i
      logical :: get_pol
      get_pol = get_pol_v(data%vec(data%index(i)))
    end function get_pol

    !------------------------------------------------------

    subroutine set_pol(data,i,pol)
      type(hilbert_vector_collection_t), intent(inout) :: data
      integer, intent(in) :: i
      logical, intent(in) :: pol
      call set_pol_v(data%vec(i),pol)
    end subroutine set_pol

    !-----------------------------------------------------

    function get_cutoff(data,i)
      type(hilbert_vector_collection_t), intent(in)    :: data
      integer, intent(in)                              :: i
      real(dp)                                         :: get_cutoff

      get_cutoff = get_cutoff_v(data%vec(data%index(i)))
    end function get_cutoff

    !--------------------------------------------------------

    function get_rad_func_p(data,i)
      type(hilbert_vector_collection_t), intent(in) :: data
      integer, intent(in)                           :: i !Includes all the copies!
      type(rad_func_t), pointer                     :: get_rad_func_p
      get_rad_func_p => get_rad_func_p_v(data%vec(data%index(i)))
    end function get_rad_func_p

    !--------------------------------------------------------

     function get_rad_func(data,i)
      type(hilbert_vector_collection_t), intent(in) :: data
      integer, intent(in)                           :: i !Includes all the copies!
      type(rad_func_t)                              :: get_rad_func
      type(rad_func_t) :: rad_tmp
      rad_tmp = get_rad_func_v(data%vec(data%index(i)))
      call rad_copy(rad_tmp , get_rad_func)
      call rad_dealloc(rad_tmp)
    end function get_rad_func

    !--------------------------------------------------------

    subroutine set_vector(data,vec,i)
      type(hilbert_vector_collection_t), intent(inout) :: data
      type(hilbert_vector_t), intent(in)               :: vec
      integer, intent(in)                              :: i

      if(vector_is_initialized(data%vec(i))) call destroy_vector(data%vec(i))
      call copy_vector(vec,data%vec(i))
    end subroutine set_vector

    !---------------------------------------------------------
    
    function get_vector(data,i) result(vec)
      type(hilbert_vector_collection_t), intent(in)    :: data
      integer, intent(in)                              :: i

      type(hilbert_vector_t)                           :: vec     
      call copy_vector(data%vec(i),vec)
    end function get_vector

    !---------------------------------------------------------
    

  end module hilbert_vector_collection
