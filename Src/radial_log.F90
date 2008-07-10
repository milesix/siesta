! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!

module radial_log
  use precision, only : dp
  use sys, only : die
  use flib_spline, only : evaluate_spline
  use radial_util, only : nrmpal
  use schro, only : schro_eq, polarization, rc_vs_e, rphi_vs_e,energ_deriv
  use radial_logGrid, only : logGrid_t, log_grid_alloc, log_grid_dealloc, &
       log_grid_copy

  
  use radial_log_low_level, only : log_rad_func_t, log_rad_alloc, &
       log_rad_dealloc, log_rad_cutoff, log_rad_get_default_length, &
       log_rad_setup_d2, log_rad_set_default_length, log_rad_copy_grid, &
       log_rad_set_maximum_length, log_rad_min, log_rad_get, &
       log_rad_set_origin, log_rad_update, log_rad_copy, eps, min_func_val,&
       maximum_length,restricted_grid, log_rad_get_length

  use radial_log_io, only : log_rad_dump_file,log_rad_dump_funcs_ascii, &
       log_rad_dump_ascii, log_rad_read_ascii_unformatted, &
       log_rad_read_ascii_formatted
  implicit none
 
contains

  subroutine log_rad_broadcast(func)

    use parallel, only : Node,  nodes
#ifdef MPI
  use mpi_siesta
#endif
  type(log_rad_func_t),pointer :: func
  
#ifndef MPI
  end subroutine log_rad_broadcast
#else

   integer :: funcLength, MPIerror

  if (Nodes .eq. 1) return

  if (node .eq. 0) funcLength = size(func%f)

  call MPI_Bcast(funcLength,1,MPI_integer,0,MPI_Comm_World,MPIError)

  if (node .ne. 0) allocate(func%f(1:funcLength),func%d2(1:funcLength),func%grid)

  call MPI_Bcast(func%grid%a,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(func%grid%b,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(func%grid%r,funclength,MPI_double_precision,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(func%grid%drdi,funclength,MPI_double_precision,0,MPI_Comm_World,MPIerror)

  call MPI_Bcast(func%f,funclength,MPI_double_precision,0,MPI_Comm_World,MPIerror)  
  call MPI_Bcast(func%d2,funclength,MPI_double_precision,0,MPI_Comm_World,MPIerror)


  end subroutine log_rad_broadcast
#endif
  !-----------------------------------------
    
  function log_rad_divide_by_4pir2(func,update) result(divided)
    type(log_rad_func_t),intent(in) :: func
    logical, intent(in)             :: update

    type(log_rad_func_t) :: divided,tmp
    real(dp) ::r2,pi
    integer :: ir

    call log_rad_copy(func,tmp)

    pi=acos(-1.0_dp) 
    
    do ir=2,size(tmp%f)
       r2=tmp%grid%r(ir)**2
       r2 = 4.0_dp*pi*r2
       tmp%f(ir)=tmp%f(ir) /(r2)
    enddo
    r2=tmp%Grid%r(2)/(tmp%Grid%r(3)-tmp%Grid%r(2))
    tmp%f(1)=tmp%f(2)-(tmp%f(3)-tmp%f(2))*r2

    if(update)then
       divided = log_rad_update(tmp)
    else
       call log_rad_setup_d2(tmp)
       call log_rad_copy(tmp,divided)
    endif

    call log_rad_dealloc(tmp)
  end function log_rad_divide_by_4pir2

  !-----------------------------------------
 function log_rad_divide_by_r_l_1(func,l,lambda) result(div)
    !Divide func by r**(l+1)
    type(log_rad_func_t), intent(in) :: func
    integer, intent(in) :: l
    real(dp), intent(in) :: lambda

    type(log_rad_func_t) :: div
    integer :: ir

    call log_rad_copy(func,div)

    do ir=2,size(div%f)
       div%f(ir)=div%f(ir)/(div%grid%r(ir)**(l+1)*sqrt(lambda**(2*l+3)))
    end do
    div%f(1)=div%f(2)
  end function log_rad_divide_by_r_l_1

  !-------------------------------------------------------

  function log_rad_get_value_from_ir(func,ir) result(y)
    integer,intent(in) :: ir
    type(log_rad_func_t)   :: func

    real(dp)           :: y
    
    y=0.0_dp
    if (ir > size(func%f)) then
       call die("radial_log: log_get_value_from_ir: ir too big")
    else
       y=func%f(ir)
    endif

  end function log_rad_get_value_from_ir

  !-------------------------------------------
  
  function log_rad_get_r_from_value(func,val) result(r)
    type(log_rad_func_t)   :: func
    real(dp), intent(in) :: val

    real(dp)           :: r
    integer            :: ir,nrc,nsp

    nrc = size(func%f)

    nsp =nrc
    do ir = nrc, 2, -1
       if (func%f(ir) > val) then
          if (ir == size(func%f)) then ! borderline case
             nsp = ir
          else
             ! Choose closest point
             if ( (func%f(ir)-val) > (val-func%f(ir+1)) ) then
                nsp = ir + 1
             else
                nsp = ir   
             endif
          endif
          exit
       endif
    enddo
    r = func%f(nsp)

  end function log_rad_get_r_from_value

  !-------------------------------------------
  
!!$  subroutine log_rad_set_value_at_ir(func,ir,y)
!!$    type(log_rad_func_t), intent(inout)   :: func
!!$    integer,intent(in)                    :: ir
!!$    real(dp), intent(in)                  :: y
!!$    
!!$    if (ir > size(func%f)) then
!!$       call die("radial_log: log_get_value_from_ir: ir too big")
!!$    else
!!$       func%f(ir) = y
!!$    endif
!!$
!!$  end subroutine log_rad_set_value_at_ir

  !-------------------------------------------

  function log_rad_get_r(func,ir) result(r)
    type(log_rad_func_t), intent(in) :: func
    integer,intent(in)               :: ir
   

    real(dp)           :: r

    r = 0.0_dp

    if (ir > size(func%Grid%r))then
       call die("radial_log: log_get_value_from_ir: ir too big")
    else
       r=func%grid%r(ir)
    endif

  end function log_rad_get_r

  !--------------------------------------------

 function log_rad_get_ir(func,r) result(ir)
   type(log_rad_func_t), intent(in)   :: func
   real(dp),intent(in)                :: r

   integer :: ir
   
   ir =  nint(log(r/func%grid%b+1.0_dp)/func%grid%a)+1

  end function log_rad_get_ir

  !--------------------------------------------
function log_rad_multiply_by_rl(func,l) result (mult)
    type(log_rad_func_t),intent(in) :: func
    integer, intent(in) :: l

    type(log_rad_func_t) :: mult,tmp

    real(dp) ::r
    integer :: ir
    
    call log_rad_copy(func,tmp)
    
    do ir=2,size(tmp%f)
       r=tmp%grid%r(ir)
       tmp%f(ir)=tmp%f(ir)*r**l
    enddo
    tmp%f(1) = tmp%f(2)
    mult = log_rad_update(tmp)
    call log_rad_dealloc(tmp)

  end function log_rad_multiply_by_rl

  !-----------------------------------------

  function log_rad_reparametrize(func,rmax,a,b) result(new_func)
    !
    !        Interpolate values into new grid, given by a and b
    !
    !        Typical new values:  a = 5x10-4, b=10

    type(log_rad_func_t), intent(in) :: func
    real(dp), intent(in) :: rmax
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    type(log_rad_func_t) :: new_func

    !Internal vars
    real(dp)  :: new_rmax, rpb, ea, ea2, rr
    integer   :: ir, new_nrval, j, old_n
    real(dp), dimension(:), pointer   :: new_r, old_r, old_f,new_f,old_d2
    type(logGrid_t) :: new_grid

    !First store old values.
    old_n = size(func%f)
    allocate(old_r(1:old_n),old_f(1:old_n),old_d2(1:old_n))
    old_r = func%Grid%r
    old_f = func%f
    old_d2 = func%d2

    !New grid
    new_rmax = func%grid%r(size(func%f)) !cutoff
    print *, "Reparametrization. new rmax: ",new_rmax
    if(new_rmax > rmax) new_rmax = rmax
    rpb=b
    ea=exp(a)
    ea2=1.0d0
    ir = 0
    do 
       rr = b*(ea2-1.0d0)
       if (rr > rmax) exit
       ir = ir + 1
       rpb=rpb*ea
       ea2=ea2*ea
    enddo
    new_nrval = ir
    
    !
    allocate(new_r(1:new_nrval),new_f(1:new_nrval))
    !print *, "Reparametrization. New nrval: ", new_nrval
    call log_rad_set_maximum_length(new_nrval)
    rpb=b
    ea=exp(a)
    ea2=1.0d0
    do ir = 1, new_nrval
       new_r(ir) = b*(ea2-1.0d0)
       rpb=rpb*ea
       ea2=ea2*ea
    enddo

    new_grid = log_grid_alloc(new_nrval,a,b,new_r)
         
   
    do j = 1, new_nrval
       call log_rad_get(func,new_r(j),new_f(j))
    enddo
    
    !print *, "reparam: rc, value(rc)=",new_r(new_nrval),new_f(new_nrval)
    call log_rad_alloc(new_func,new_f,new_grid)
    deallocate(old_r,old_f,old_d2,new_r,new_f)
  end function log_rad_reparametrize

  !-----------------------------------------  

  function log_rad_sum(f1,f2,rmax) result(sum)
    !Sums two functions. The grid used is the
    !longest one.
    type(log_rad_func_t), intent(in) :: f1
    type(log_rad_func_t), intent(in) :: f2
    real(dp), intent(in), optional :: rmax

    type(log_rad_func_t) :: sum,tmp

    integer  :: i, length, length1, length2,irmax
    real(dp) :: r,y1,y2 

    length1 = log_rad_get_length(f1)
    length2 = log_rad_get_length(f2)
    if(present(rmax)) irmax = log_rad_get_ir(f2,rmax)
   
    
    if (length1 .ge. length2) then
       length = length2
       call log_rad_copy(f1,tmp)
       do i=1,length
          r = log_rad_get_r(f1,i)
          y2 = log_rad_get_value_at_r(f2,r)
          tmp%f(i) = tmp%f(i)+y2
       enddo
    else
       length = length1
       call log_rad_copy(f2,tmp)
       do i=1,length
          r = log_rad_get_r(f2,i)
          y1 = log_rad_get_value_at_r(f1,r)
          tmp%f(i) = tmp%f(i)+y1
       enddo
    endif

    if(present(rmax)) then
       length=irmax
       call log_rad_alloc(sum,tmp%f(1:irmax),tmp%grid)
    else
       sum= log_rad_update(tmp)
    endif
    !call log_rad_setup_d2(sum)
    call log_rad_dealloc(tmp)
  end function log_rad_sum

  !-----------------------------------------

  function log_rad_sum_function(func,func2,r0,rmax) result(sum)
    type(log_rad_func_t), intent(in) :: func
    real(dp), intent(in) :: r0, rmax

    interface
       function func2(x) result(y)
         use precision
         real(dp), intent(in) :: x
         real(dp)             :: y
       end function func2
    end interface

    type(log_rad_func_t) :: sum,tmp

    real(dp) :: f2,r
    integer :: ir

    call log_rad_copy(func,tmp)

    do ir=1,size(tmp%f)
       r = log_rad_get_r(tmp,ir)
       if(r .ge. r0 .or. r .le. rmax) then
          f2 = func2(r)
          tmp%f(ir)=tmp%f(ir)+f2
       endif
    enddo

    sum = log_rad_update(tmp)
    call log_rad_dealloc(tmp)

  end function log_rad_sum_function

  !-----------------------------------------

 function log_rad_sum_value(func,value) result (sum)
    type(log_rad_func_t),intent(in) :: func
    real(dp), intent(in) :: value

    type(log_rad_func_t) :: sum,tmp

    real(dp) ::r
    integer :: ir

    call log_rad_copy(func,tmp)
    
    do ir=2,size(tmp%f)
       r=func%grid%r(ir)
       tmp%f(ir)=tmp%f(ir)+value
    enddo
   
    sum= log_rad_update(tmp)
    call log_rad_dealloc(tmp)
  end function log_rad_sum_value

  !-----------------------------------------

  function log_rad_integral(func) result (integral)
    type(log_rad_func_t),intent(in) :: func

    real(dp) :: integral
    integer :: ir
    
    integral = 0.0_dp

    do ir=1,size(func%f)
       integral=integral+func%grid%drdi(ir)*func%f(ir)
    enddo

  end function log_rad_integral

  !-----------------------------------------

  subroutine log_rad_interpol_values(func,rmin,rmax,values) 
    type(log_rad_func_t)   :: func
    real(dp), intent(in) :: rmin, rmax
    real(dp), intent(inout) :: values(:)
    
    integer :: nvalues,i
    real(dp) :: r,dfdr,dr

    nvalues = size(values)
    if (rmax > log_rad_cutoff(func)) then
       call die("radial_log: log_interpol_values: rmax too big")
    else
       dr = real((rmax-rmin)/size(values))
       r=rmin
       do i=1,nvalues
          call evaluate_spline(func%grid%r,func%f,func%d2,size(func%f),r,values(i),dfdr) 
          r = r + dr
       enddo
    endif

  end subroutine log_rad_interpol_values

  !----------------------------------------------------------------

  function log_rad_get_value_at_r(func,r) result(var)
    type(log_rad_func_t)    :: func
    real(dp), intent(in)    :: r

    real(dp) :: var
    
 
    real(dp) ::dfdr

    if (r > log_rad_cutoff(func)) then
       call die("radial_lin: lin_interpol_values: rmax too big")
    !elseif( r == 0.0_dp)then
    !   var = log_rad_get_value_from_ir(func,1)
    else
       call evaluate_spline(func%grid%r,func%f,func%d2,size(func%f),r,var,dfdr)
    endif

  end function log_rad_get_value_at_r

  !----------------------------------------------------------------

  subroutine log_rad_normalize_r_l_1(rad_func,l)
    type(log_rad_func_t), intent(inout) :: rad_func
    integer, intent(in) :: l

    !Internal vars
    real(dp) :: norm, r
    integer  :: ir

    type(logGrid_t), pointer :: grid

    grid => rad_func%grid

    norm=0.0_dp
    do ir=2,size(rad_func%f)
       r = grid%r(ir)
       norm=norm+grid%drdi(ir)*(rad_func%f(ir)*r**(l+1))**2
    enddo

    !Normalize (if they aren't)
    !eps=1.0E-4
    if(abs(norm-1.0_dp).gt.eps) then
       do ir=1,size(rad_func%f)
          rad_func%f(ir)=rad_func%f(ir)/sqrt(norm)            
       enddo
    endif

  end subroutine log_rad_normalize_r_l_1

  !-------------------------------------------------------

 

  function log_rad_kinetic_energy(rad_func,l) result(ekin)
    type(log_rad_func_t), intent(in) :: rad_func
    integer, intent(in) ::  l
    real(dp)            :: ekin

    real(dp)                          :: r,r1,r2,d2fdi2,dfdi,dr,d2fdr2
    integer                           :: ir

    type(logGrid_t), pointer :: logGrid

    logGrid => rad_func%grid
    ekin=0.0_dp
    !Kinetic
    do ir=2,size(rad_func%f)-1
       r=logGrid%r(ir)
       r1=logGrid%r(ir-1)
       r2=logGrid%r(ir+1)
       d2fdi2=(rad_func%f(ir-1)*r1**(l+1)+rad_func%f(ir+1)*r2**(l+1)-2.0d0*rad_func%f(ir)*r**(l+1))
       dfdi=0.5d0*(rad_func%f(ir+1)*r2**(l+1)-rad_func%f(ir-1)*r1**(l+1))
       dr=logGrid%drdi(ir)
       d2fdr2= ((-logGrid%a)*dfdi +  d2fdi2)/dr**2
       ekin=ekin+dr*rad_func%f(ir)*r**(l+1)*(-d2fdr2) &
            +dr*l*(l+1)*(rad_func%f(ir)*r**l)**2
    enddo
  end function log_rad_kinetic_energy

  !------------------------------------------------------------------------------

  function log_rad_potential_energy(rad_func,potential,l) result(energy)
    type(log_rad_func_t), intent(in) :: rad_func
    type(log_rad_func_t), intent(in) :: potential
    integer, intent(in) :: l

    real(dp) :: energy

    real(dp) :: r,f,pot,dfdr
    integer  :: ir
    type(logGrid_t), pointer :: logGrid

    logGrid => rad_func%grid

    !Potential:
    energy = 0.0_dp
    do ir=1,size(rad_func%f)
       r=logGrid%r(ir)
       f=rad_func%f(ir)
       call evaluate_spline(potential%grid%r,potential%f,potential%d2,size(potential%f),r,pot,dfdr)
       
       energy = energy + logGrid%drdi(ir) * (pot)*(f*r**(l+1))**2

    enddo
  end function log_rad_potential_energy

  !--------------------------------------------------------------------

  function log_rad_schro(pseudo,ve,rc,l,nnodes,nprin,chg,energy) result(integral)
    type(log_rad_func_t), intent(in) :: pseudo
    type(log_rad_func_t), intent(in) :: ve     !pseudo of valence electrons.
    real(dp), intent(in) :: rc
    integer, intent(in) :: l
    integer, intent(in) :: nnodes !number of nodes
    integer, intent(in) :: nprin  ! l + number of shells with this n
    real(dp), intent(in) :: chg   !charge
    real(dp), intent(out) :: energy !energy
    
    type(log_rad_func_t) :: integral

    !Internal vars
    real(dp), dimension(:), allocatable ::rphi, s
    real(dp) :: dnrm
    integer :: nrc,nrval,i,nrc_new, ir
    type(logGrid_t), pointer :: grid
   
    nrval = size(pseudo%f) 
        
    nrc = log_rad_get_ir(pseudo,rc)

    if(restricted_grid) nrc=nrc+1-mod(nrc,2)
    if(nrc > nrval) nrc=nrval
    if(nrc .gt. maximum_length) nrc=maximum_length
    if(nrval .gt. maximum_length) nrval=maximum_length

    allocate(rphi(1:nrc),s(1:nrc))
    s(2:nrc)=pseudo%grid%drdi(2:nrc)**2
    s(1) = s(2)
    grid => pseudo%grid
    
    do i=1,nrc
       write(22,*) grid%r(i),pseudo%f(i)
       write(23,*) grid%r(i),ve%f(i)
       write(24,*) grid%r(i),s(i)
       write(25,*) grid%r(i),grid%drdi(i)
    end do

    call schro_eq(chg,grid%r(1:nrval),pseudo%f(1:nrval),ve%f(1:nrval),s, &
         grid%drdi(1:nrval),nrc, l,grid%a,grid%b,nnodes,nprin,energy,rphi)

    dnrm=0.0
    do ir=1,nrc
       dnrm=dnrm+grid%drdi(ir)*rphi(ir)**2
    end do
    dnrm=sqrt(dnrm)
    
    do ir=1,nrc
       rphi(ir)=rphi(ir)/dnrm
    enddo
    nrc_new = nrc
    do i=nrc,2,-1
       if(abs(rphi(i)) .gt. min_func_val)then
          nrc_new=i+1
          exit
       endif
    enddo

    call log_rad_alloc(integral,rphi(1:nrc_new),grid)
    deallocate(rphi,s)
    nullify(grid)
    
  end function log_rad_schro

  !---------------------------------------------------------------------------

  function log_rad_energy_deriv(psi,vps,ve,l,energy) result(ederiv)
    integer, intent(in) :: l
    real(dp), intent(in)            :: energy
    type(log_rad_func_t), intent(in) ::psi, vps, ve

    type(log_rad_func_t) :: ederiv
    integer :: nrc
    
    real(dp), dimension(:), pointer :: psidev

    stop "lograd: energ_deriv: check grid is the same"

    nrc = size(psi%f)
    allocate(psidev(1:nrc))

    call energ_deriv(psi%grid%a,psi%grid%r,psi%f,vps%f,ve%f, &
         psi%grid%drdi,nrc,l,energy,psidev,nrc)

    call log_rad_alloc(ederiv,psidev,psi%grid)
    deallocate(psidev)

  endfunction log_rad_energy_deriv

  !----------------------------------------------------------------
 
  function log_rad_integral_vs_value(vps,ve,l,energy,rmax) result(func)
    type(log_rad_func_t), intent(in):: vps, ve
    integer, intent(in)  :: l
    real(dp), intent(in) ::energy, rmax

    type(log_rad_func_t) :: func

    !Internal vars
    type(logGrid_t) :: grid
    real(dp), pointer :: rphi(:)
    integer :: nr

    !nr = size(rad_func%f)
    stop "radial_log: rphi_vs_e: nr?"
    allocate(rphi(1:nr))
    
    call log_grid_copy(vps%grid,grid)
    call rphi_vs_e(grid%a,grid%b,grid%r,vps%f,ve%f,nr,l,energy,rphi,rmax)
    call log_rad_alloc(func,rphi,grid)
    deallocate(rphi)
  endfunction log_rad_integral_vs_value

!-----------------------------------------------------------------------

  function log_rad_multiply_each_value(func,val) result(multiplied)
    type(log_rad_func_t), intent(in) :: func
    real(dp), intent(in)             :: val

    type(log_rad_func_t) :: multiplied
    integer :: i, length

    length = size(func%f)

    call log_rad_copy(func,multiplied)

    do i=1,size(func%f)
       multiplied%f(i)=func%f(i)*val
    enddo
    !call log_rad_update(multiplied)
    call log_rad_setup_d2(func)
  end function log_rad_multiply_each_value

  !---------------------------------------

  function log_rad_multiply(fu1,fu2) result(mult)
    type(log_rad_func_t), intent(in) :: fu1,fu2

    type(log_rad_func_t) :: mult

    integer  :: i, length, length1, length2
    real(dp) :: r,y1,y2 

    length1 = log_rad_get_length(fu1)
    length2 = log_rad_get_length(fu2)
   
    if (length1 .ge. length2) then
       length = length2
       call log_rad_copy(fu1,mult)
       do i=1,length
          r = log_rad_get_r(fu1,i)
          y2 = log_rad_get_value_at_r(fu2,r)
          mult%f(i) = mult%f(i)*y2
       enddo
    else
       length = length1
       call log_rad_copy(fu2,mult)
       do i=1,length
          r = log_rad_get_r(fu2,i)
          y1 = log_rad_get_value_at_r(fu1,r)
          mult%f(i) = mult%f(i)*y1
       enddo
    endif
    !mult = log_rad_update(mult)
    call log_rad_setup_d2(mult)
  end function log_rad_multiply

  !---------------------------------------

  function log_rad_hartree(func) result(vhartree)
    type(log_rad_func_t), intent(in) :: func

    type(log_rad_func_t) :: vhartree

    real(dp), pointer :: potential(:),s(:),rho(:)
    real(dp):: rpb,ea
    integer :: length,i

    !length=size(func%f)
    length=size(func%grid%r)
    allocate(rho(1:length),potential(1:length),s(1:length))

    rpb=func%grid%b
    ea=exp(func%grid%a)
    do i=1,length
       s(i)=sqrt(func%grid%a*rpb)
       rpb=rpb*ea       
       if(i<=size(func%f)) then
          rho(i)=func%f(i)
       else
          rho(i)=0.0_dp
       endif
    enddo
   
    call vhrtre(rho,potential,func%grid%r,func%grid%drdi,s,length,func%grid%a)
    call log_rad_alloc(vhartree,potential,func%grid)
    deallocate(potential,s,rho)
  end function log_rad_hartree

  !-------------------------------------------------------

  function log_rad_vxc(func,irel,exc) result(vxc)
    type(log_rad_func_t), intent(in) :: func
    integer, intent(in) :: irel
    real(dp), intent(out) :: exc 

    type(log_rad_func_t) :: vxc

    integer :: length,ispin,i
    real(dp) :: dx, dc, ec, ex
    real(dp), pointer :: vxc_values(:),rho(:)

    length = size(func%grid%r)
    allocate(vxc_values(1:length),rho(1:length))
    ispin = 1

    do i=1,length
       if(i<=size(func%f)) then
          rho(i)=func%f(i)
       else
          rho(i)=0.0_dp
       endif
    enddo

    call atomxc(irel,length,length,func%grid%r,ispin,rho,ex,ec,dx,dc,vxc_values)

    exc=ex+ec

    vxc_values(1)=vxc_values(2)
    call log_rad_alloc(vxc,vxc_values,func%grid)

    !do i=1,length
    !   print *, i, vxc%grid%drdi(i), vxc%f(i),vxc%d2(i)
    !enddo
    deallocate(vxc_values,rho)
  end function log_rad_vxc

  !-------------------------------------------------------------

  function log_rad_rc_vs_e(func,pseudo,ve,l,el,nnodes) result(rc)
    type(log_rad_func_t), intent(in) :: func, pseudo,ve
    real(dp), intent(in) :: el
    integer, intent(in) :: l, nnodes

    real(dp) :: rc

    call rc_vs_e(func%grid%a,func%grid%b,func%grid%r,pseudo%f,&
         ve%f,log_rad_get_length(func), l, el, nnodes,rc)
  end function log_rad_rc_vs_e

  !--------------------------------------------------------------

  function log_rad_split_gauss(first_z,rc,lambda,l) result(gauss)
    type(log_rad_func_t), intent(in) :: first_z
    real(dp), intent(in) :: rc, lambda
    integer, intent(in) :: l
    type(log_rad_func_t) :: gauss

    real(dp) :: fac, cons, gexp, phi, pi,r
    integer  :: i, nrc,ir
    real(dp), allocatable :: g(:)

    pi = acos(-1.0_dp)

    !With spligauss option, compression factor must be taken
    !as the gaussian exponent

    if(lambda.le.0.0d0) then 
       write(6,'(/a,/a,a)') &
            'SPLITGAUSS: ERROR: within SPLITGAUSS option the compression ', &
            'SPLITGAUSS: ERROR: factors for all the augmentation functions', &
            ' must be explicitely specified' 
       call die
    endif


    gexp = lambda
    gexp = 1.0_dp/(gexp**2)

    fac=1.0d0
    do i=0,l
       fac=(2*i+1)*fac
    enddo


    cons=sqrt(pi)*fac/(2.0d0**(l+2))
    cons=cons/((2.0d0*gexp)**(l+1.5d0))
    cons=1.0d0/sqrt(cons)

    !Generate the gaussian
    nrc = log_rad_get_ir(first_z,rc)-1
    allocate(g(1:nrc))
    do ir=1,nrc
       r=first_z%grid%r(ir)
       phi=cons*exp((-gexp)*r**2)
       !dnrm=dnrm+drdi(ir)*(phi*r**(l+1))**2
       g(ir)=phi
    enddo
 
    call log_rad_alloc(gauss,g,first_z%grid)
  end function log_rad_split_gauss

  !--------------------------------------------------------------

  function log_rad_split_scan_tail(rphi) result(split)
    type(log_rad_func_t), intent(in) :: rphi

    integer :: ir, nrc

    real(dp), pointer :: drdi(:), rnrm(:),split_table(:)
    type(log_rad_func_t) :: split

    nrc = size(rphi%f)
    allocate(rnrm(1:nrc),split_table(1:nrc))
    drdi => rphi%grid%drdi

    do ir = 1, nrc
      rnrm(ir)=drdi(ir)*rphi%f(ir)**2
    enddo

    split_table(1:nrc) = sqrt(max(1.0_dp-rnrm(1:nrc),0.0_dp))

    call log_rad_alloc(split,split_table,rphi%grid)
    deallocate(rnrm,split_table)
  end function log_rad_split_scan_tail

  !--------------------------------------------------------------
  
  function log_rad_polarization(func,pseudo,ve,rc,l,energy) result(pol)
    type(log_rad_func_t), intent(in) :: func
    type(log_rad_func_t), intent(in) ::pseudo,ve
    real(dp), intent(in) :: rc
    integer, intent(in) :: l
    real(dp), intent(in) :: energy

    type(log_rad_func_t) :: pol
    real(dp), pointer :: psipol(:)
    integer :: ir, nrc, nrc_new
    real(dp) :: dnrm

    nrc = log_rad_get_ir(func,rc)-1
    allocate(psipol(1:nrc))
    call polarization(func%grid%a,func%grid%r,func%f,pseudo%f,ve%f, &
         func%grid%drdi,nrc,l,energy,psipol,nrc)

    dnrm=0.0
    do ir=1,nrc
       dnrm=dnrm+func%grid%drdi(ir)*psipol(ir)**2
    end do
    dnrm=sqrt(dnrm)
    
    do ir=1,nrc
       psipol(ir)=psipol(ir)/dnrm
    enddo

    do ir=nrc,2,-1
       if(abs(psipol(ir)) .gt. min_func_val)then
          nrc_new=ir+1
          exit
       endif
    enddo

    call log_rad_alloc(pol,psipol(1:nrc_new),func%grid)
    deallocate(psipol)
  end function log_rad_polarization

  !---------------------------------------------------------------------------

  subroutine log_rad_set_restricted_grid(restricted)
    logical, intent(in) :: restricted
    restricted_grid = restricted
  end subroutine log_rad_set_restricted_grid

  !---------------------------------------------------------------------------

  subroutine log_rad_fit_parabola(func,rmatch,l, cons1,cons2,parab_norm)
     ! Fits C1 and C2 to match a parabola to rphi at
     ! point of index nsp
     type(log_rad_func_t), intent(in) :: func
     real(dp), intent(in) :: rmatch
     integer, intent(in)   :: l
     real(dp), intent(out) :: cons1, cons2
     real(dp), intent(out) :: parab_norm
    

     real(dp), pointer  :: r(:), drdi(:)
     real(dp), pointer  :: rphi(:)
  

     integer   :: nsp
     real(dp) :: rsp, frsp, dfrsp
     
     nsp = log_rad_get_ir(func,rmatch)

     if(restricted_grid) nsp=nsp+1-mod(nsp,2)

     r=>func%grid%r
     drdi=>func%grid%drdi
     rphi=>func%f

     rsp = r(nsp)
     frsp=rphi(nsp)/rsp
     dfrsp=0.5d0*(rphi(nsp+1)/r(nsp+1) &
          -rphi(nsp-1)/r(nsp-1))
     dfrsp=dfrsp/drdi(nsp)
     
     cons1= 0.5d0*(dfrsp*rsp-l*frsp)/(rsp**(l+2))
     cons2= frsp/(rsp**l)-cons1*rsp**2
     call nrmpal(cons1,cons2,rsp,l,parab_norm)

  end subroutine log_rad_fit_parabola

  !--------------------------------------------------------

  subroutine log_rad_zero(func)
    type(log_rad_func_t), intent(out) :: func
    !func%cutoff = 0.0_dp
    !func%n      = 0.0_dp
    func%f      = 0.0_dp
    func%d2     = 0.0_dp
    !call log_grid_zero(func%grid)
  end subroutine log_rad_zero

  !-------------------------------------------

end module radial_log

