!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

!Radial (1D) functions defined in a linnear grid
module radial_lin
  use precision, only : dp
  use m_recipes, only : spline, splint
  use m_radfft, only : radfft
  use m_filter, only : filter, kcPhi
  use xml
  use sys

  Implicit none

  public lin_rad_func_t, linGrid_t
  public lin_grid_alloc
  public lin_rad_alloc, lin_rad_broadcast, lin_rad_copy, lin_rad_cutoff, lin_rad_dealloc
  public lin_rad_default_length, lin_rad_dump_file, lin_rad_dump_fft_file
  public lin_rad_dump_ascii, lin_rad_dump_xml
  public lin_rad_filter,  lin_rad_get
  public lin_rad_FFT, lin_rad_dump_fft_ascii, lin_rad_dump_fft_xml
  public lin_rad_get_value_from_ir
  public lin_rad_get_filter_cutoff, lin_rad_get_grid, lin_rad_get_length
  public lin_rad_read_ascii, lin_rad_set_origin, lin_rad_zero, lin_rad_dump_netcdf
  private

  type linGrid_t
     !private
     integer  :: n
     real(dp) :: delta
  end type linGrid_t

  type lin_rad_func_t
     private
     real(dp)                        :: delta   ! real space distance between points
     real(dp)                        :: kdelta  ! 
     real(dp), dimension(:), pointer :: f  => Null()  ! real space data 
     real(dp), dimension(:), pointer :: fft => Null() ! Fft of the data
     real(dp), dimension(:), pointer :: d2 => Null()  ! Second derivative in real space
  end type lin_rad_func_t
  
  ! INTEGER NTBMAX    : Maximum number of points in the tables defining
  !                     orbitals,projectors and local neutral-atom 
  !                     pseudopotential.

  integer, parameter, public  :: ntbmax =  500 !It must be 2**n for the fft grid.

  
contains

  subroutine lin_rad_alloc(func,delta,values,yp1)
    type(lin_rad_func_t), intent(inout) :: func
    real(dp), intent(in) :: delta
    real(dp), intent(in) :: values(:)
    real(dp), optional, intent(in) :: yp1

    func%delta = delta
    allocate(func%f(ntbmax+1),func%d2(1:ntbmax+1),func%fft(1:ntbmax+1))
    func%f(1:ntbmax) = values(1:ntbmax)
    func%f(ntbmax+1) = 0.0_dp
    func%fft(ntbmax+1) = 0.0_dp
    func%d2(1:ntbmax+1) = 0.0_dp

    if (present(yp1)) then
       call lin_rad_update(func,yp1)
    else
       call lin_rad_update(func)
    endif
  end subroutine lin_rad_alloc

  !------------------------------------------

  subroutine lin_grid_alloc(grid,n,delta)
    type(linGrid_t), intent(inout) :: grid
    integer, intent(in) :: n
    real(dp), intent(in) :: delta

    grid%n = n
    grid%delta = delta
  end subroutine lin_grid_alloc

  !------------------------------------------

  subroutine lin_rad_copy(src,dest)
    type(lin_rad_func_t), intent(in) :: src
    type(lin_rad_func_t), intent(out) :: dest
    
    call lin_rad_alloc(dest,src%delta,src%f)
    dest%fft = src%fft
  end subroutine lin_rad_copy

  !--------------------------------------------------

  subroutine lin_rad_broadcast(func)
    use parallel, only : Node,  nodes
#ifdef MPI
  use mpi_siesta
#endif

  type(lin_rad_func_t), pointer :: func

#ifndef MPI
  end subroutine lin_rad_broadcast
#else

  integer :: funcLength, MPIerror

  if (Nodes .eq. 1) return

  if (node .eq. 0) funcLength = size(func%f)

  call MPI_Bcast(funcLength,1,MPI_integer,0,MPI_Comm_World,MPIError)

  if (node .ne. 0) allocate(func%f(1:funcLength),func%fft(1:funcLength),func%d2(1:funcLength))

  call MPI_Bcast(func%delta,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)

  call MPI_Bcast(func%kdelta,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)

  call MPI_Bcast(func%f,funclength,MPI_double_precision,0,MPI_Comm_World,MPIerror)

  call MPI_Bcast(func%fft,funclength,MPI_double_precision,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(func%d2,funclength,MPI_double_precision,0,MPI_Comm_World,MPIerror)

  end subroutine lin_rad_broadcast
#endif
  !--------------------------------------------------

   function lin_rad_cutoff(func)
    real(dp)                         :: lin_rad_cutoff
    type(lin_rad_func_t), intent(in) :: func
    lin_rad_cutoff = dble( (ntbmax-1)*func%delta )
  end function lin_rad_cutoff

  !------------------------------------------------------------

  subroutine lin_rad_dealloc(func)
    type(lin_rad_func_t), intent(inout)::func
    deallocate(func%f)
    deallocate(func%d2)
    deallocate(func%fft)
  endsubroutine lin_rad_dealloc

  !-----------------------------------------

  function lin_rad_default_length() result(leng)
    integer :: leng
    leng = ntbmax
  end function lin_rad_default_length

  !---------------------------------------------

   subroutine lin_rad_dump_ascii(op,lun,header)
    type(lin_rad_func_t)    :: op
    integer lun
    logical,intent(in),optional :: header

    integer j
    logical print_header

    print_header = .true.
    if(present(header))then
       print_header = header
    endif

    if(print_header)then
       write(lun,'(i4,2g22.12,a)') ntbmax,op%delta,dble(op%delta*(ntbmax-1)),&
            " #npts, delta, cutoff"
    endif

    do j=1,ntbmax
       write(lun,'(2g22.12)') (j-1)*op%delta, op%f(j) !,op%d2(j)
    enddo
  end subroutine lin_rad_dump_ascii

  !----------------------------------------------------------

  subroutine lin_rad_dump_file(func,label)
    type(lin_rad_func_t), intent(in) :: func
    character(len=*), intent(in) :: label

    integer :: ir,io
    real(dp) :: r
    io=70
    open(unit=io,file=label,status="replace")
    do ir=1, ntbmax
       r= dble(ir-1)*func%delta
       !print *, "y=",func%f(ir)
       write(70,'(2f15.9)') r,func%f(ir)
    end do
    close(io)
  end subroutine lin_rad_dump_file

  !----------------------------------------------------------

  subroutine lin_rad_dump_netcdf(op,lun,netcdf_id)

#ifndef CDF
    type(lin_rad_func_t), intent(in) :: op
    integer, intent(in)              :: lun, netcdf_id

  end subroutine lin_rad_dump_netcdf

#else
    use netcdf
    type(lin_rad_func_t), intent(in) :: op
    integer, intent(in)              :: lun, netcdf_id

    integer :: iret
    print *, "Size of op%f in dump_netcdf: ", size(op%f)
    iret = nf90_put_var(lun,netcdf_id,op%f(1:), start=(/1/),count=(/size(op%f)/))
    call check(iret)
  contains
    subroutine check(status)

      integer, intent(in):: status
      if (status .ne. nf90_noerr) then
         print  *, trim(nf90_strerror(status))
         call die()
      endif
    end subroutine check

  end subroutine lin_rad_dump_netcdf

#endif
  !----------------------------------------------------------

   subroutine lin_rad_dump_fft_ascii(op,lun,header)
    type(lin_rad_func_t)    :: op
    integer :: lun
    logical,intent(in),optional :: header

    integer :: j,ik,nk
    logical :: print_header
    real(dp) :: k(1:ntbmax),dk,rmax,pi,kcut
    pi = acos(-1._dp)

    print_header = .true.
    if(present(header))then
       print_header = header
    endif

    if(print_header) write(lun,'(i4,g22.12,a)') ntbmax,op%delta," #npts, delta"

    rmax = lin_rad_cutoff(op)
    dk = pi / rmax
    k = (/(ik,ik=0,ntbmax-1)/)*dk
    kcut=lin_rad_kcutoff(op)
    nk=2*kcut/dk+2

    write(lun,'(2g22.12)') (k(j), op%fft(j),j=1,ntbmax/2)
    
  end subroutine lin_rad_dump_fft_ascii

  !-----------------------------------------------------------
  
  subroutine lin_rad_dump_fft_file(op,file)
    type(lin_rad_func_t), intent(in) :: op
    character(len=*), intent(in)     :: file

    real(dp) :: pi,dk, rmax, k(1:ntbmax), kcut
    integer :: j, nk, ik

    pi = acos(-1._dp)
    
    rmax = lin_rad_cutoff(op)
    kcut = lin_rad_kcutoff(op)
    dk = pi / rmax
    k = (/(ik,ik=0,ntbmax-1)/)*dk
    nk=2*kcut/dk+2

    open(unit=99,file=file,status="replace")

    write(99,'(2g22.12)') (k(j), op%fft(j),j=1,ntbmax/2)
    close(99)
  end subroutine lin_rad_dump_fft_file

  !-----------------------------------------------------------

  subroutine lin_rad_dump_fft_xml(op,lun)
    type(lin_rad_func_t), intent(in) :: op
    integer, intent(in)              :: lun

    integer :: j

    write(lun,'(a)') '<fft radfunc>'
    call xml_dump_element(lun,'npts',str(ntbmax))
    call xml_dump_element(lun,'delta',str(op%delta))
    write(lun,'(a)') '<data>'
    do j=1,ntbmax/2
       write(lun,'(2g22.12)') (j-1)*op%delta, op%fft(j)
    enddo
    write(lun,'(a)') '</data>'
    write(lun,'(a)') '</fft radfunc>'
  end subroutine lin_rad_dump_fft_xml

  !-------------------------------------------------------------


  subroutine lin_rad_dump_xml(op,lun)
    type(lin_rad_func_t)    :: op
    integer lun
    integer j

    write(lun,'(a)') '<radfunc>'
    call xml_dump_element(lun,'npts',str(ntbmax))
    call xml_dump_element(lun,'delta',str(op%delta))
    call xml_dump_element(lun,'cutoff',str(ntbmax*(op%delta-1)))
    write(lun,'(a)') '<data>'
    do j=1,ntbmax
       write(lun,'(2g22.12)') (j-1)*op%delta, op%f(j)
    enddo
    write(lun,'(a)') '</data>'
    write(lun,'(a)') '</radfunc>'
  end subroutine lin_rad_dump_xml

  !----------------------------------------------------------

  subroutine lin_rad_fft(func,l)
      type(lin_rad_func_t), intent(inout) :: func
      integer, intent(in) :: l

    !Internal vars
    real(dp) :: pi,rc
    integer  :: i


    pi = 4.0_dp * ATAN( 1.0_dp )

    !Siesta stores the orbs divided by r**l
    func%fft = func%f

    if(l/=0) forall(i=1:ntbmax) func%fft(i)=func%fft(i)*(func%delta*(i-1))**l
    rc = lin_rad_cutoff(func)

    !Call fft subroutine
    call radfft(l,ntbmax,rc,func%fft,func%fft)

  end subroutine lin_rad_fft

  !-----------------------------------------------------------------

  function lin_rad_filter(rad_func,l,factor,norm_opt,kc) result(filtered)
    type(lin_rad_func_t), intent(in) :: rad_func
    integer, intent(in) :: l,norm_opt
    real(dp), intent(in) :: factor 
    real(dp), intent(in) :: kc !filter cutoff

    type(lin_rad_func_t) :: filtered

    real(dp), allocatable, dimension(:) :: r,y

    integer :: i,length

    length = ntbmax

    allocate(r(1:length),y(1:length))

    do i=2,length
       r(i)=rad_func%delta*(i-1)
       y(i)=(rad_func%f(i)*r(i)**l)*factor
    enddo
    y(1)=y(2)
    r(1)=0.0_dp

    call filter(l,length,r,y,kc,norm_opt)

    do i=2,length
       y(i)=y(i)/(r(i)**l*factor)
    enddo
    y(1)=y(2)
    
    call lin_rad_alloc(filtered,rad_func%delta,y)

    deallocate(r,y)

  end function lin_rad_filter

  !-----------------------------------------------------
 
  function lin_rad_get_filter_cutoff(rad_func,l,factor,etol) result(kc)
    !Given a tolerance in the kinetic energy this function
    !returns the corresponding reciprocal space cutoff.
    !See module filter.f90
    type(lin_rad_func_t), intent(in) :: rad_func
    integer, intent(in)              :: l
    real(dp), intent(in)             :: factor, etol
    real(dp) :: kc 

    real(dp), allocatable, dimension(:) :: r,y

    integer :: i,length

    length = ntbmax

    allocate(r(1:length),y(1:length))

    do i=1,length
       r(i)=rad_func%delta*(i-1)
       y(i)=(rad_func%f(i)*r(i)**(l))*factor
    enddo

    kc=kcPhi(l,length,r,y,etol)
    
    deallocate(r,y)

  end function lin_rad_get_filter_cutoff

  !-----------------------------------------------------

  subroutine lin_rad_get(func,r,fr,dfdr)
    type(lin_rad_func_t), intent(in) :: func
    real(dp), intent(in)         :: r
    real(dp), intent(out)        :: fr
    real(dp), intent(out)        :: dfdr

    if (size(func%f) .eq. 0) then
       fr = 0._dp
       dfdr = 0._dp
    else
       if (r.gt. lin_rad_cutoff(func)) then
          write(6,*) 'Attempt to evaluate beyond cutoff', r
          fr = 0._dp
          dfdr = 0._dp
       else
          call splint(func%delta,func%f,func%d2,size(func%f),r,fr,dfdr)
       endif
    endif

  end subroutine lin_rad_get

  !--------------------------------------------

  function lin_rad_get_grid(func) result(grid)
    type(lin_rad_func_t) :: func

    type(linGrid_t) :: grid
    grid%n     = ntbmax
    grid%delta = func%delta
  end function lin_rad_get_grid

  !--------------------------------------------

  function lin_rad_kcutoff(func) result(kcut)
    type(lin_rad_func_t), intent(in) :: func
    real(dp) :: kcut

    real(dp) :: rmax,pi, dk

    pi = acos(-1._dp)
    
    rmax = lin_rad_cutoff(func)
    dk = pi/rmax
    kcut = dble(dk*(ntbmax-1))

  end function lin_rad_kcutoff

  !--------------------------------------------

  function lin_rad_get_length() result(length)
    integer :: length

    length = ntbmax
  end function lin_rad_get_length

  !--------------------------------------------

    function lin_rad_get_value_from_ir(func,ir) result(y)
    integer,intent(in) :: ir
    type(lin_rad_func_t)   :: func

    real(dp)           :: y
    
    y=0.0_dp

    if (ir > ntbmax) then
       call die("radial_lin: log_get_value_from_ir: ir too big")
    else
       y=func%f(ir)
    endif

  end function lin_rad_get_value_from_ir

  !--------------------------------------------

  !subroutine lin_rad_interpol_values(func,rmin,rmax,values) 
  !  type(lin_rad_func_t)   :: func
  !  real(dp), intent(in) :: rmin, rmax
  !  real(dp), intent(inout) :: values(:)
    
  !  stop "lin_interpol_values not implemented"

  !end subroutine lin_rad_interpol_values

  !---------------------------------------------

  subroutine lin_rad_read_ascii(op,lun)
    type(lin_rad_func_t)    :: op
    integer, intent(in) :: lun

    integer :: j, npts
    real(dp) :: r1,r2,new_delta,delta,dummy
    real(dp), pointer :: values(:)

    read(lun,*) npts, delta, dummy !size(op%f)*(op%delta-1)
    allocate(values(1:npts))
    r1=0.0_dp
    r2=0.0_dp
    read(lun,*) r1, values(1)
    read(lun,*) r2, values(2)
    delta = r2-r1
    do j=3,npts
       read(lun,*) r2,values(j)
       new_delta = r2-r1
       if(abs(abs(new_delta)-abs(delta)) .gt. 1E-4)then
          call die("lin_rad_read_ascii: not a lin rad func!")
       end if
       r1=r2
       delta=new_delta
    enddo
    call lin_rad_alloc(op,delta,values)
  end subroutine lin_rad_read_ascii

  !--------------------------------------------------------------------    

  !     Do not use yet... interface in need of fuller specification
  !
  function lin_rad_rvals(func) result (r)
    real(dp), dimension(:), pointer :: r
    type(lin_rad_func_t), intent(in) :: func

    integer i

    nullify(r)
    if (size(func%f) .eq. 0) return
    allocate(r(size(func%f)))
    do i=1,size(func%f)
       r(i) = func%delta *(i-1)
    enddo
  end function lin_rad_rvals

  !-------------------------------------------------------------------

  subroutine lin_rad_setup_d2(func,yp1)
    !Set up second derivative (for spline interpol) of a radial function
    type(lin_rad_func_t), intent(inout) :: func
    real(dp),intent(in),optional  :: yp1

    real(dp) :: yp1def, ypn

    if (size(func%f) .eq. 0) return

    if (present(yp1))then
       yp1def = yp1
    else
       yp1def = huge(1._dp)
    endif

    ypn = huge(1._dp)
    call spline(func%delta,func%f,ntbmax,yp1def,ypn,func%d2)

  end subroutine lin_rad_setup_d2

  !-----------------------------------------------------------------

  subroutine lin_rad_set_origin(func,value)
    !Set value at origin
    type(lin_rad_func_t), intent(inout) :: func
    real(dp),             intent(in) :: value

    func%f(1) = value
  end subroutine lin_rad_set_origin

  !-------------------------------------------------------

   subroutine lin_rad_update(func,yp1)
    type(lin_rad_func_t), intent(inout) :: func
    real(dp), optional :: yp1

    if (present(yp1))then
       call lin_rad_setup_d2(func,yp1)
    else
       call lin_rad_setup_d2(func)
    endif

  end subroutine lin_rad_update

  !------------------------------------------

   subroutine lin_rad_update_rc(rad_func,rc)
    type(lin_rad_func_t), intent(inout) :: rad_func
    real(dp), intent(in) :: rc

    rad_func%delta = rc/ntbmax
  end subroutine lin_rad_update_rc

  !--------------------------------------------------

   subroutine lin_rad_zero(func)
    type(lin_rad_func_t), intent(out) :: func
    
    func%delta=0.0_dp
    func%f = 0.0_dp
    func%d2= 0.0_dp
    func%fft = 0.0_dp
  end subroutine lin_rad_zero

  !----------------------------------------------------------------
  
end module radial_lin
