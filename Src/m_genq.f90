!
module m_genq_lua
  use flook

  private
  
  integer, parameter :: dp = selected_real_kind(10,100)
  
  type(luaState), public :: lua
  public :: genq_lua_initialize

  integer, public :: natoms, nqs
  
  integer, allocatable, public  :: isa(:)
  real(dp), allocatable, public    :: q(:) ! could use setter
  real(dp), allocatable, public :: xa(:,:)
  
  logical, public :: genq_lua_initialized = .false.

CONTAINS

  subroutine genq_lua_initialize()
    ! This routine reads the Lua script, registers a few
    ! functions to exchange data (Fortran arrays/Lua tables)
    ! and gets the number of generalized coordinates and atoms
    ! that the Lua side knows about
    
  if (genq_lua_initialized) RETURN
    
  ! Initialize the @lua environment
  call lua_init(lua)

  ! Register a couple of functions to pass information back and
  ! forth between @lua.
  call lua_register(lua, 'fortran_set_dimensions', script_get_dimensions )
  call lua_register(lua, 'fortran_get_qs', script_set_qs )
  call lua_register(lua, 'fortran_set_xa', script_get_xa )
  call lua_register(lua, 'fortran_set_isa', script_get_isa )

  call lua_run(lua, 'genq_flook.lua' )

  call lua_run(lua, code = 'fortran_set_dimensions()')
  
  allocate(q(nqs))
  allocate(xa(3,natoms), isa(natoms))

  genq_lua_initialized = .true.
  
end subroutine genq_lua_initialize

  function script_set_qs(state) result(nret)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int
    
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(luaState) :: lua
    type(luaTbl) :: tbl

    call lua_init(lua,state)

    ! open global table in variable genq
    tbl = lua_table(lua,'genq')

    ! Set the variables to the genq table:
    call lua_set(tbl,'q',q)

    call lua_close_tree(tbl)

    ! this function returns nothing
    nret = 0
    
  end function script_set_qs

  function script_get_xa(state) result(nret)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int
    
    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(luaState) :: lua
    type(luaTbl) :: tbl

    call lua_init(lua,state)

    ! open global table in variable genq
    tbl = lua_table(lua,'genq')

    ! Get the variables from the genq table:
    call lua_get(tbl,'xa',xa)

    call lua_close_tree(tbl)

    ! this function returns nothing
    nret = 0
    
  end function script_get_xa

function script_get_isa(state) result(nret)
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int

  ! Define the state
  type(c_ptr), value :: state
  ! Define the in/out
  integer(c_int) :: nret

  type(luaState) :: lua
  type(luaTbl) :: tbl

  call lua_init(lua,state)

  ! open global table in variable genq
  tbl = lua_table(lua,'genq')

  ! Get the variables from the genq table:
  call lua_get(tbl,'species',isa)

  call lua_close_tree(tbl)

  ! this function returns nothing
  nret = 0

end function script_get_isa

  function script_get_dimensions(state) result(nret)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int
    
    ! Define the state
    type(c_ptr), value :: state
    ! Define the in/out
    integer(c_int) :: nret

    type(luaState) :: lua
    type(luaTbl) :: tbl

    call lua_init(lua,state)

    ! open global table in variable genq
    tbl = lua_table(lua,'genq')

    ! Get the variables from the genq table:
    call lua_get(tbl,'nqs',nqs)
    call lua_get(tbl,'natoms',natoms)

    call lua_close_tree(tbl)

    ! this function returns nothing
    nret = 0
    
  end function script_get_dimensions

end module m_genq_lua

module m_genq

  
  use precision, only: dp
  use sys,       only: die

  implicit none

  type, public :: genq_t
     integer                        :: nvars = 0
     real(dp), allocatable          :: q(:)
     character(len=32), allocatable :: qname(:)
  end type genq_t

  public :: geom_from_genq
  public :: cartF_to_genqF
  public :: read_genq

  type(genq_t), public, save  :: genq

  private

  CONTAINS

  subroutine geom_from_genq(M,q,xa,spec_no)
    ! Given the values of M generalized coordinates
    ! q(1), q(2), ... q(M)
    ! this routine specifies the values of the crystal coordinates xa

    use m_genq_lua, only: lua
    use m_genq_lua, only: genq_q => q, genq_xa => xa
    use m_genq_lua, only: genq_isa => isa
    use m_genq_lua, only: genq_lua_nqs => nqs
    use m_genq_lua, only: genq_lua_initialize
    use m_genq_lua, only: genq_lua_initialized

    use flook, only: lua_run
    
    integer, parameter    :: p = selected_real_kind(10,100)

    integer, intent(in)                         :: M
    real(p), intent(in)                         :: q(M)
    real(p), intent(inout)                      :: xa(:,:)
    integer, intent(inout), optional            :: spec_no(:)


    logical :: want_species_numbers

    if (.not. genq_lua_initialized) then
       call genq_lua_initialize()
       if (genq_lua_nqs /= M) then
          call die("nqs mismatch in genq")
       endif
    endif

    want_species_numbers = .false.
    if (present(spec_no)) then
       want_species_numbers = .true.
    endif

    ! Set the Fortran array and have Lua read it and put it
    ! in its genq.q table
    genq_q(:) = q(:)
    call lua_run(lua, code = 'fortran_get_qs()')

    ! Execute Lua function to generate xa and species tables
    ! in genq table
    call lua_run(lua, &
         code = 'genq.xa, genq.species = coords_from_qs( genq.q )' )

    ! Get xa array in the genq_lua module and copy it
    call lua_run(lua, code = 'fortran_set_xa()')
    xa(:,:) = genq_xa(:,:)

    if (want_species_numbers) then
       ! Get isa array in the genq_lua module and copy it
       call lua_run(lua, code = 'fortran_set_isa()')
       spec_no(:) = genq_isa(:)
    endif
    
  end subroutine geom_from_genq

  subroutine cartF_to_genqF(M,q,ucell,na,fa,gF,ftol,epsgF)

    ! Given the values of M generalized coordinates
    ! q(1), q(2), ... q(M)
    ! and a routine geom_from_genq
    ! this routine computes \partial xa \partial q
    ! and the projection of the cartesian forces on the gen coords

    ! We assume linear dependencies
    ! In most (all?) cases, the derivatives will be simple integers

    Integer, intent(in)                    :: M
    real(dp), intent(in)                   :: q(M)
    real(dp), intent(in)                   :: ucell(3,3)
    integer, intent(in)                    :: na
    real(dp), intent(in)                   :: fa(3,na)
    real(dp), intent(out)                  :: gF(M)
    real(dp), intent(in)                   :: ftol
    real(dp), intent(out)                  :: epsgF(M)


    real(dp), allocatable           :: xa(:,:)
    real(dp), allocatable           :: xaf(:,:)
    real(dp), allocatable           :: dxdq(:,:)
    real(dp), allocatable           :: qq(:)

    real(dp)       :: h  = 0.01
    real(dp)       :: dxdq_cart(3), epsf(3)

    integer       :: ia, i, j

    allocate(xa(3,na))
    allocate(xaf(3,na))
    allocate(dxdq(3,na))
    allocate(qq(M))

    ! Compute partial derivatives by simple rule

    call geom_from_genq(M,q,xa)

    qq = q

    do i = 1, M
       gF(i) = 0.0_dp
       epsgF(i) = 0.0_dp
       !
       qq(i) = q(i) + h
       call geom_from_genq(M,qq,xaf)
!       dxdq(1:3,1:na) = nint((xaf-xa) / h)
       dxdq(1:3,1:na) = (xaf-xa) / h
       !
       !  Project the cartesian forces to the generalized coords
       !  to get the "generalized forces"
       !
       do ia = 1, na
          dxdq_cart(:) = matmul(ucell,dxdq(:,ia))
          gF(i) = gF(i) + dot_product(fa(:,ia),dxdq_cart(:))
          ! Estimate the force when almost converged
          do j=1,3
             epsf(j) = sign(ftol,fa(j,ia))
          enddo
          ! Estimate a tolerance for gF(i)
          epsgF(i) = epsgF(i) + dot_product(epsf(:),dxdq_cart(:))
       enddo
       !
       qq(i) = q(i)
    enddo

    deallocate(xa,xaf,dxdq,qq)

  end subroutine cartF_to_genqF

  subroutine read_genq(genq)
    ! For now, this routine reads from file GENQ

    use m_mpi_utils, only:  broadcast
    use parallel,    only:  Node

    type(genq_t) :: genq
    
    integer :: var_u, iostat, i

    if (Node .eq. 0) then
       open(unit=var_u,file="GENQ",form="formatted",status="old", &
                  position="rewind",action="read")

       read(var_u, fmt=*, iostat=iostat) genq%nvars
       if (iostat /= 0) then
          call die("Cannot read nvars from GENQ")
       endif
    endif

    call broadcast(genq%nvars)

    allocate(genq%q(genq%nvars),genq%qname(genq%nvars))
  
    if (Node .eq. 0) then
       do i = 1, genq%nvars
          read(var_u, fmt=*, iostat=iostat) genq%q(i), genq%qname(i)
          if (iostat /= 0) then
             call die("ERROR while reading q")
          endif
       enddo
       close(var_u)
    endif

    call broadcast(genq%q)
    !! No support yet.. call broadcast(genq%qname)

end subroutine read_genq

end module m_genq

    


