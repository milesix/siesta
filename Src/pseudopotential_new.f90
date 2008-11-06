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
module pseudopotential_new

  use sys, only: die
  use precision, only: dp
  use radial, only : rad_func_t, rad_grid_t, rad_grid_alloc, rad_alloc, &
       rad_dealloc, rad_copy, rad_reparametrize, rad_dump_funcs_ascii, &
       rad_get_grid, rad_write_ascii_formatted, rad_grid_get_a, &
       rad_grid_dealloc, rad_grid_get_b, rad_set_default_length, &
       rad_cutoff, rad_grid_get_length,rad_grid_dump_ascii_formatted

  use flib_spline, only: generate_spline, evaluate_spline
  use atom_options, only: write_ion_plot_files
  
  implicit none

  external :: io_assign, io_close

  private

  public :: pseudopotential_new_t
  public :: get_ve_val, get_ve_val_scaled, pseudo_convert
  public :: get_grid,get_pseudo_down,  pseudo_header_print, get_pseudo_valence_charge
  public :: pseudo_has_core_charge, get_core_charge, get_vlocal, pseudo_copy_valence_charge
  public :: get_pseudo_gen_valence_charge, set_rho_val, set_rho_val_scaled, set_ve_val
  public :: set_ve_val_scaled, pseudo_copy_scaled_ve, get_pseudo_method, get_pseudo_npotd
  public :: get_pseudo_nicore, pseudo_dealloc_vdown, pseudo_set_vdown, get_irel, get_icorr
  public :: get_pseudo_l, get_rho_val, get_rho_val_scaled, pseudo_copy_scaled_charge
  public :: set_pseudo_vlocal, pseudo_reparametrize

  !pubic :: pseudo_read,
  !public :: pseudo_write_formatted, pseudo_reparametrize
  !public :: read_ps_conf, pseudo_dump, 
 

  type pseudopotential_new_t
     private
     character(len=2)        :: name
     real(dp)                :: zval
     real(dp)                :: gen_zval ! Valence charge used during generation
     logical                 :: relativistic
     character(len=10)       :: correlation
     character(len=2)        :: icorr
     character(len=3)        :: irel
     character(len=4)        :: nicore = 'nc '
     character(len=10)       :: method(6)
     character(len=70)       :: text
     integer                 :: npotu
     integer                 :: npotd
     type(rad_func_t)        :: chcore !pseudo_rho_core for non linear core corrections
     type(rad_func_t)        :: chval  !pseudo_rho_val
     type(rad_func_t), pointer :: vdown(:)
     type(rad_func_t), pointer :: vup(:)

     integer, pointer        :: ldown(:)
     integer, pointer        :: lup(:)
     type(rad_func_t)        :: rho_val
     type(rad_func_t)        :: rho_val_scaled
     type(rad_func_t)        :: ve_val     !potential due to rho used KB genaration
     type(rad_func_t)        :: ve_val_scaled !potential due to rho_pao used in PAO genaration
     type(rad_func_t)        :: vlocal
  end type pseudopotential_new_t

CONTAINS

  function get_core_charge(pseudo) result(core)
    type(pseudopotential_new_t), intent(in) :: pseudo
    type(rad_func_t) :: core

    call rad_copy(pseudo%chcore,core)

  end function get_core_charge

  !-------------------------------------------------------------

  function get_grid(pseudo) result(grid)
    type(pseudopotential_new_t), intent(in) :: pseudo
    type(rad_grid_t) :: grid

    grid = rad_get_grid(pseudo%vdown(0))
    
  end function get_grid

  !-------------------------------------------------------------

  function get_irel(pseudo) result(irel)
    type(pseudopotential_new_t), intent(in) :: pseudo
    character(len=3) :: irel
    irel = pseudo%irel
  end function get_irel

  !-------------------------------------------------------------

   function get_icorr(pseudo) result(icorr)
    type(pseudopotential_new_t), intent(in) :: pseudo
    character(len=2) :: icorr
    icorr = pseudo%icorr
  end function get_icorr

  !-------------------------------------------------------------

  function get_pseudo_down(pseudo,l) result(vps)
    type(pseudopotential_new_t), intent(in) :: pseudo
    integer, intent(in) :: l
    
    type(rad_func_t) :: vps
    call rad_copy(pseudo%vdown(l),vps)
    
  end function get_pseudo_down

  !-------------------------------------------------------------

  function get_pseudo_method(pseudo) result(method)
    type(pseudopotential_new_t), intent(in) :: pseudo
    character(len=10) :: method(6)
    method = pseudo%method
  end function get_pseudo_method

  !-------------------------------------------------------------

  function get_pseudo_npotd(pseudo) result(npotd)
    type(pseudopotential_new_t), intent(in) :: pseudo
    integer :: npotd
    npotd = pseudo%npotd
  end function get_pseudo_npotd

  !-------------------------------------------------------------

  function get_pseudo_nicore(pseudo) result(nicore)
    type(pseudopotential_new_t), intent(in) :: pseudo
    character(len=4):: nicore
    nicore = pseudo%nicore
  end function get_pseudo_nicore

  !------------------------------------------------------------

  function get_ve_val(pseudo) result(ve)
    type(pseudopotential_new_t), intent(in) :: pseudo
    type(rad_func_t) :: ve
    call rad_copy(pseudo%ve_val,ve)
  end function get_ve_val

  !-------------------------------------------------------------

  function get_ve_val_scaled(pseudo) result(ve)
    type(pseudopotential_new_t), intent(in) :: pseudo
    type(rad_func_t) :: ve
    call rad_copy(pseudo%ve_val_scaled,ve)
  end function get_ve_val_scaled

  !-------------------------------------------------------------

  function get_vlocal(pseudo) result(vlocal)
    type(pseudopotential_new_t), intent(in) :: pseudo
    type(rad_func_t) :: vlocal
    call rad_copy(pseudo%vlocal,vlocal)
  end function get_vlocal

  !-------------------------------------------------------------

  function get_pseudo_l(pseudo,n) result(l)
    type(pseudopotential_new_t), intent(in) :: pseudo
    integer, intent(in)                     :: n
    integer                                 :: l

    l = pseudo%ldown(n)
  end function get_pseudo_l

  !-------------------------------------------------------------

  function get_pseudo_valence_charge(pseudo) result(zval)
    type(pseudopotential_new_t), intent(in) :: pseudo
    real(dp) :: zval
    zval = pseudo%zval
  end function get_pseudo_valence_charge

  !-------------------------------------------------------------

  function get_pseudo_gen_valence_charge(pseudo) result(gen_zval)
    type(pseudopotential_new_t), intent(in) :: pseudo
    real(dp) :: gen_zval
    gen_zval = pseudo%gen_zval
  end function get_pseudo_gen_valence_charge

  !-------------------------------------------------------------

  function get_rho_val(pseudo) result(rho)
    type(pseudopotential_new_t), intent(in) :: pseudo
    
    type(rad_func_t) :: rho
    call rad_copy(pseudo%rho_val,rho)
    
  end function get_rho_val

  !-------------------------------------------------------------

  function get_rho_val_scaled(pseudo) result(rho)
    type(pseudopotential_new_t), intent(in) :: pseudo
    
    type(rad_func_t) :: rho
    call rad_copy(pseudo%rho_val_scaled,rho)
    
  end function get_rho_val_scaled

  !-------------------------------------------------------------

  subroutine pseudo_dealloc_vdown(pseudo,l)
    type(pseudopotential_new_t), intent(inout) :: pseudo
    integer, intent(in) :: l
    call rad_dealloc(pseudo%vdown(l))
  end subroutine pseudo_dealloc_vdown

  !-------------------------------------------------------------

  subroutine pseudo_set_vdown(pseudo,l,func)
    type(pseudopotential_new_t), intent(inout) :: pseudo
    integer, intent(in) :: l
    type(rad_func_t) :: func
    call rad_copy(func,pseudo%vdown(l))
  end subroutine pseudo_set_vdown

  !-------------------------------------------------------------

  subroutine set_pseudo_vlocal(pseudo,vloc)
    type(pseudopotential_new_t), intent(inout) :: pseudo
    type(rad_func_t) :: vloc
    call rad_copy(vloc,pseudo%vlocal)
  end subroutine set_pseudo_vlocal

  !-------------------------------------------------------------


  function pseudo_has_core_charge(pseudo) result(core)
    type(pseudopotential_new_t), intent(in) :: pseudo
    logical :: core

    if(pseudo%nicore == 'nc ')then
       core = .false.
    else
       core = .true.
    endif

  end function pseudo_has_core_charge

  !-------------------------------------------------------------

  subroutine pseudo_copy_valence_charge(pseudo)
    type(pseudopotential_new_t), intent(inout) :: pseudo
    call rad_copy(pseudo%chval,pseudo%rho_val)
  end subroutine pseudo_copy_valence_charge

  !-------------------------------------------------------------

  subroutine pseudo_copy_scaled_charge(pseudo)
    type(pseudopotential_new_t), intent(inout) :: pseudo
    call rad_copy(pseudo%rho_val,pseudo%rho_val_scaled)
  end subroutine pseudo_copy_scaled_charge

  !-------------------------------------------------------------

  subroutine pseudo_copy_scaled_ve(pseudo)
    type(pseudopotential_new_t), intent(inout) :: pseudo
    call rad_copy(pseudo%ve_val,pseudo%ve_val_scaled)
  end subroutine pseudo_copy_scaled_ve

  !-------------------------------------------------------------

  subroutine pseudo_convert(old,new)
    use pseudopotential, only : pseudopotential_t
    type(pseudopotential_t), intent(in) :: old
    type(pseudopotential_new_t), intent(inout) ::new

    type(rad_grid_t) :: grid

    integer :: i

    new%name = old%name
    new%zval = old%zval
    new%gen_zval = old%gen_zval
    new%relativistic = old%relativistic
    new%correlation = old%correlation
    new%icorr = old%icorr
    new%irel  = old%irel
    new%nicore = old%nicore
    new%method = old%method
    new%text   = old%text
    new%npotu  = old%npotu
    new%npotd  = old%npotd

    grid = rad_grid_alloc(old%nrval,a=old%a,b=old%b,r=old%r)
    call rad_set_default_length(grid, old%nrval)
    call rad_alloc(new%chcore,old%chcore,grid)
    call rad_alloc(new%chval,old%chval,grid)

    if(old%npotd .gt. 0)then
       new%npotd = old%npotd
       allocate(new%ldown(0:new%npotd-1))
       allocate(new%vdown(0:new%npotd-1))
       do i=1,new%npotd
          new%ldown(i-1) = old%ldown(i)
          call rad_alloc(new%vdown(i-1),old%vdown(i,1:old%nrval),grid)
       end do
    endif

    if(old%npotu .gt. 0)then
       new%npotu = old%npotu
       allocate(new%lup(0:new%npotu-1))
       allocate(new%vup(0:new%npotu-1))
       do i=1,new%npotu
          new%lup(i-1) = old%lup(i)
          call rad_alloc(new%vup(i-1),old%vup(i,1:old%nrval),grid)
       end do
    endif
    
  end subroutine pseudo_convert

  !-------------------------------------------------------------

  subroutine set_rho_val(vps,rho_val)
    type(pseudopotential_new_t), intent(inout) :: vps
    type(rad_func_t), intent(in) :: rho_val

    call rad_copy(rho_val,vps%rho_val)
  end subroutine set_rho_val

  !-------------------------------------------------------------

  subroutine set_rho_val_scaled(vps,rho_val_scaled)
    type(pseudopotential_new_t), intent(inout) :: vps
    type(rad_func_t), intent(in) :: rho_val_scaled

    call rad_copy(rho_val_scaled,vps%rho_val_scaled)
  end subroutine set_rho_val_scaled

  !-------------------------------------------------------------

  subroutine set_ve_val(vps,ve_val)
    type(pseudopotential_new_t), intent(inout) :: vps
    type(rad_func_t), intent(in) :: ve_val

    call rad_copy(ve_val,vps%ve_val)
  end subroutine set_ve_val

  !-------------------------------------------------------------

  subroutine set_ve_val_scaled(vps,ve_val_scaled)
    type(pseudopotential_new_t), intent(inout) :: vps
    type(rad_func_t), intent(in) :: ve_val_scaled

    call rad_copy(ve_val_scaled,vps%ve_val_scaled)
  end subroutine set_ve_val_scaled

  !-------------------------------------------------------------
  

!   subroutine pseudo_read(label,p)
!     character(len=*), intent(in)   :: label
!     type(pseudopotential_new_t)                    :: p

!     !       PS information can be in a .vps file (unformatted)
!     !       or in a .psf file (formatted)

!     character(len=40) :: fname
!     logical           :: found

!     fname  = trim(label) // '.vps'
!     inquire(file=fname, exist=found)
!     if (found) then
!        call pseudo_read_unformatted(fname,p)
!     else
!        fname = trim(label) // '.psf'
!        inquire(file=fname, exist=found)
!        if (found) then
!           call pseudo_read_formatted(fname,p)
!        else
!           write(6,'(/,2a,a20,/)') 'read_pseudo: ERROR: ',&
!                'Pseudopotential file not found: ', fname
!           call die
!        endif
!     endif
!     if (write_ion_plot_files) &
!          call pseudo_dump(trim(p%name) // ".psdump",p)
!   end subroutine pseudo_read
  !
!   subroutine pseudo_read_unformatted(fname,p)
!     character(len=*), intent(in) :: fname
!     type(pseudopotential_new_t)  :: p
!     type(rad_grid_t)             :: grid

!     integer io_ps, i,nr, nrval,j
!     real(dp) :: r2,v,r_int,r3,v2,v3, a, b
!     real(dp), allocatable :: values(:),r(:)

! 8030 format(4(g20.12))

!     call io_assign(io_ps)
!     open(io_ps,file=fname,form='unformatted',status='unknown')
!     write(6,'(3a)') 'Reading pseudopotential information ',&
!          'in unformatted form from ', trim(fname)

!     read(io_ps) p%name, p%icorr, p%irel, p%nicore,(p%method(i),i=1,6),&
!          p%text,p%npotd, p%npotu, nr, b, a, p%zval

!     !
!     !       Old style vps files should have the right info in text.
!     !
!     call read_ps_conf(p%irel,p%npotd-1,p%text,p%gen_zval)

!     nrval = nr + 1
!     allocate(values(1:nrval),r(1:nrval))

!     !Read the grid
!     read(io_ps,8030) (r(j),j=2,nrval)
!     r(1) = 0.d0

!     grid = rad_grid_alloc(nrval,a=a,b=b,r=r)
!     call rad_set_default_length(grid,nrval)

!     if (p%npotd.gt.0) then
!        allocate(p%vdown(0:p%npotd-1)) 
!        allocate(p%ldown(0:p%npotd-1))
!     endif

!     do i=0,p%npotd-1
!        read(io_ps) p%ldown(i)
!        read(io_ps) (values(j),j=2,nrval)
!        values(1)=values(2)
!        call rad_alloc(p%vdown(i),values,Grid)
!     enddo

!     if (p%npotu.gt.0) then
!        allocate(p%vup(0:p%npotu-1))
!        allocate(p%lup(0:p%npotu-1))
!     endif
!     do i=0,p%npotu-1
!        read(io_ps) p%lup(i)
!        read(io_ps) (values(j),j=2,nrval)
!        values(1)=values(2)
!        call rad_alloc(p%vup(i),values,Grid)
!     enddo

!     !Read chcore
!     read(io_ps) (values(j),j=2,nrval)

!     r2 = r(2)
!     r3 = r(3)
!     r_int  = r2/(r2-r3)

!     !Compute the value of chore at the origin
!     v2 = values(2)
!     v3 = values(3)
!     v  = v2-r_int*(v3-v2)
!     values(1) = v
!     call rad_alloc(p%chcore,values,Grid)

!     !The same for chval
!     read(io_ps) (values(j),j=2,nrval)
!     v2 = values(2)
!     v3 = values(3)
!     v  = v2-r_int*(v3-v2)
!     call rad_alloc(p%chval,values,Grid)

!     call io_close(io_ps)
!     deallocate(values,r)
!     call rad_grid_dealloc(grid)
!   end subroutine pseudo_read_unformatted
!   !----
!   subroutine pseudo_read_formatted(fname,p)
!     character(len=*), intent(in) :: fname
!     type(pseudopotential_new_t)  :: p

!     integer io_ps, i, ios, nr, nrval,ir,j
!     character(len=70) dummy
!     type(rad_grid_t) :: grid
!     real(dp) :: gen_zval_inline,v,v2,v3,r_int,r2,r3,a,b
!     real(dp), allocatable :: values(:),r(:)

!     call io_assign(io_ps)
!     open(io_ps,file=fname,form='formatted',status='unknown')
!     write(6,'(3a)') 'Reading pseudopotential information ',&
!          'in formatted form from ', trim(fname)

! 8000 format(1x,i2)
! 8005 format(1x,a2,1x,a2,1x,a3,1x,a4)
! 8010 format(1x,6a10,/,1x,a70)
! 8015 format(1x,2i3,i5,4g20.12)
! 8030 format(4(g20.12))
! 8040 format(1x,a)

!     read(io_ps,8005) p%name, p%icorr, p%irel, p%nicore
!     read(io_ps,8010) (p%method(i),i=1,6), p%text
!     read(io_ps,8015,iostat=ios) p%npotd, p%npotu, nr, b, a, p%zval,&
!          gen_zval_inline
!     if (ios < 0) gen_zval_inline = 0.0_dp

!     call read_ps_conf(p%irel,p%npotd-1,p%text,p%gen_zval)

!     !
!     !       (Some .psf files contain an extra field corresponding
!     !       to the ps valence charge at generation time. If that
!     !       field is not present, the information has to be decoded
!     !       from the "text" variable.
!     !
!     !       "Zero" pseudos have gen_zval = 0, so they need a special case.

!     if (p%gen_zval == 0.0_dp) then
!        if (gen_zval_inline == 0.0_dp) then
!           if (p%method(1) /= "ZEROPSEUDO")  call die("Cannot get gen_zval")
!        else
!           p%gen_zval = gen_zval_inline
!        endif
!     endif

!     nrval = nr + 1
!     allocate(values(1:nrval),r(1:nrval))

!     read(io_ps,8040) dummy

!     read(io_ps,8030) (r(j),j=2,size(values))
!     r(1) = 0.0_dp

!     grid = rad_grid_alloc(nrval,a=a,b=b,r=r)
!     call rad_set_default_length(grid,nrval)

!     if (p%npotd.gt.0) then
!        allocate(p%vdown(0:p%npotd-1))
!        allocate(p%ldown(0:p%npotd-1))
!     endif

!     do i=0,p%npotd-1
!        read(io_ps,8040) dummy 
!        read(io_ps,8000) p%ldown(i)
!        read(io_ps,8030) (values(ir), ir=2,nrval)
!        values(1)=values(2)
!        call rad_alloc(p%vdown(i),values,Grid)                
!     enddo

!     if (p%npotu.gt.0) then
!        allocate(p%vup(0:p%npotu-1))
!        allocate(p%lup(0:p%npotu-1))
!     endif
!     do i=0,p%npotu-1
!        read(io_ps,8040) dummy 
!        read(io_ps,8000) p%lup(i)
!        read(io_ps,8030)(values(ir), ir=2,nrval) 
!        values(1)=values(2)
!        call rad_alloc(p%vup(i),values,Grid) 
!     enddo

!     !read chcore
!     read(io_ps,8040) dummy
!     read(io_ps,8030) (values(ir), ir=2,nrval)

!     r2 = r(2)
!     r3 = r(3)
!     r_int  = r2/(r2-r3)

!     !Compute the value of chore at the origin
!     v2 = values(2)
!     v3 = values(3)
!     v  = v2-r_int*(v3-v2)
!     values(1) = v
!     call rad_alloc(p%chcore,values,Grid)

!     !The same for chval
!     read(io_ps,8040) dummy
!     read(io_ps,8030) (values(ir), ir=2,nrval) 

!     v2 = values(2)
!     v3 = values(3)
!     v  = v2-r_int*(v3-v2)
!     values(1) = v
!     call rad_alloc(p%chval,values,Grid)

!     call io_close(io_ps)

!     deallocate(values,r)
!     call rad_grid_dealloc(grid)

!   end subroutine pseudo_read_formatted
!   !------
!   !----
   subroutine pseudo_write_formatted(fname,p)
     character(len=*), intent(in) :: fname
     type(pseudopotential_new_t), intent(in)     :: p

     integer :: io_ps, i, nr

     real(dp) :: a,b
     type(rad_grid_t) :: grid

     call io_assign(io_ps)
     open(io_ps,file=fname,form='formatted',status='unknown',&
          action="write",position="rewind")
     write(6,'(3a)') 'Writing pseudopotential information ',&
          'in formatted form to ', trim(fname)

 8000 format(1x,i2)
 8005 format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010 format(1x,6a10,/,1x,a70)
 8015 format(1x,2i3,i5,4g20.12)
! !8030 format(4(g20.12))
8040 format(1x,a)

     grid = rad_get_grid(p%vdown(0))
     a = rad_grid_get_a(grid)
     b = rad_grid_get_b(grid)
     nr = rad_grid_get_length(grid)

     write(io_ps,8005) p%name, p%icorr, p%irel, p%nicore
     write(io_ps,8010) (p%method(i),i=1,6), p%text
     write(io_ps,8015) p%npotd, p%npotu, nr, b, a, p%zval,&
          p%gen_zval

     write(io_ps,8040) "Radial grid follows"
     call rad_grid_dump_ascii_formatted(grid,io_ps)


     do i=0,p%npotd-1
        write(io_ps,8040) "Down potential follows (l on next line)"
        write(io_ps,8000) p%ldown(i)
        call rad_write_ascii_formatted(p%vdown(i),io_ps)
     enddo

     do i=1,p%npotu-1
        write(io_ps,8040) "Up potential follows (l on next line)"
        write(io_ps,8000) p%lup(i)
        call rad_write_ascii_formatted(p%vup(i),io_ps)
     enddo

     write(io_ps,8040) "Core charge follows"
     call rad_write_ascii_formatted(p%chcore,io_ps)
     write(io_ps,8040) "Valence charge follows"
     call rad_write_ascii_formatted(p%chval,io_ps)

     call io_close(io_ps)
   end subroutine pseudo_write_formatted
!   !--------
   subroutine pseudo_dump(fname,p)
     !
     !       Column-oriented output
     !
     character(len=*), intent(in) :: fname
     type(pseudopotential_new_t), intent(in)     :: p

     integer io_ps, i

     call io_assign(io_ps)
     open(io_ps,file=fname,form='formatted',status='unknown',&
          action="write",position="rewind")
     write(6,'(3a)') 'Dumping pseudopotential information ',&
          'in formatted form in ', trim(fname)


     do i = 1, p%npotd-1
        call rad_dump_funcs_ascii(p%vdown(i),p%chval,io_ps)
     enddo

     call io_close(io_ps)
   end subroutine pseudo_dump

!   !------

  subroutine pseudo_header_print(lun,p)
    integer, intent(in) :: lun
    type(pseudopotential_new_t)  :: p

    integer :: i

8005 format(1x,a2,1x,a2,1x,a3,1x,a4)
8010 format(1x,6a10,/,1x,a70)

    write(lun,'(a)') '<pseudopotential_header>'
    write(lun,8005) p%name, p%icorr, p%irel, p%nicore
    write(lun,8010) (p%method(i),i=1,6), p%text
    write(lun,'(a)') '</pseudopotential_header>'

  end subroutine pseudo_header_print

!   !
!   subroutine read_ps_conf(irel,lmax,text,chgvps)
!     !
!     !     Attempt to decode the valence configuration used for
!     !     the generation of the pseudopotential
!     !     (At least, the valence charge)

!     character(len=3), intent(in)  :: irel
!     integer, intent(in)           :: lmax
!     character(len=70), intent(in) :: text
!     real(dp), intent(out)         :: chgvps

!     integer  :: l, itext
!     real(dp) :: ztot, zup, zdown, rc_read
!     character(len=2) :: orb

!     chgvps=0.0_dp

!     if(irel.eq.'isp') then
!        write(6,'(/,2a)') 'Pseudopotential generated from an ', &
!             'atomic spin-polarized calculation'

!        write(6,'(/,a)') 'Valence configuration for pseudopotential generation:'

!        do l=0,min(lmax,3)
!           itext=l*17
!           read(text(itext+1:),err=5000,fmt=8080) orb, zdown, zup, rc_read
! 8080      format(a2,f4.2,1x,f4.2,1x,f4.2)
!           chgvps = chgvps + zdown + zup
!           write(6,8085) orb, zdown, zup, rc_read
! 8085      format(a2,'(',f4.2,',',f4.2,') rc: ',f4.2)
!        enddo

!     else
!        if(irel.eq.'rel') then
!           write(6,'(/,2a)')'Pseudopotential generated from a ',&
!                'relativistic atomic calculation'
!           write(6,'(2a)') 'There are spin-orbit pseudopotentials ',&
!                'available'
!           write(6,'(2a)') 'Spin-orbit interaction is not included in ',&
!                'this calculation'
!        endif

!        write(6,'(/,a)') 'Valence configuration for pseudopotential generation:'

!        do l=0,min(lmax,3)
!           itext=l*17
!           read(text(itext+1:),err=5000,fmt=8090)orb, ztot, rc_read
! 8090      format(a2,f5.2,4x,f5.2)
!           chgvps = chgvps + ztot
!           write(6,8095) orb, ztot, rc_read
! 8095      format(a2,'(',f5.2,') rc: ',f4.2)
!        enddo

!     endif
!     return

! 5000 continue       ! Error return: set chgvps to zero

!   end subroutine read_ps_conf
!   !
!   !***********---------------------------------------------
!   !
   subroutine pseudo_reparametrize(p,a,b)
     !
     !        Interpolate values into new grid, given by a and b
     !
     !        Typical new values:  a = 5x10-4, b=10

     type(pseudopotential_new_t)  :: p
     real(dp), intent(in)         :: a, b


     !--------------------------------------------------------------------
     integer :: i
     type(rad_func_t)                :: rad_tmp

     rad_tmp = rad_reparametrize(p%chcore,a,b)
     call rad_dealloc(p%chcore)
     call rad_copy(rad_tmp,p%chcore)
     call rad_dealloc(rad_tmp)
    
     rad_tmp = rad_reparametrize(p%chval,a,b)   
     call rad_dealloc(p%chval)
     call rad_copy(rad_tmp,p%chval)
     call rad_dealloc(rad_tmp)

     do i=0,p%npotd-1
        rad_tmp = rad_reparametrize(p%vdown(i),a,b)
        call rad_dealloc(p%vdown(i))
        call rad_copy(rad_tmp,p%vdown(i))
        call rad_dealloc(rad_tmp)
        print '(a,i2,f10.4)', "pseudo_reparam: l,rc=",i,rad_cutoff(p%vdown(i))
     enddo

     if (p%npotu > 0) then           
        do i=0,p%npotu-1     ! Only executed if npotu > 0 ...
           rad_tmp = rad_reparametrize(p%vup(i),a,b)
           call rad_dealloc(p%vup(i))   
           call rad_copy(rad_tmp,p%vup(i))
           call rad_dealloc(rad_tmp)
        enddo
     endif

     call pseudo_write_formatted(trim(p%name)// ".Reparam.psf",p)
     if(write_ion_plot_files) &
          call pseudo_dump(trim(p%name) // ".Reparam.psdump",p)

   end subroutine pseudo_reparametrize

end module pseudopotential_new



