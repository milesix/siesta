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
module atom_basis_io
  !
  !     Support for dumping and reading PAO and KB information from
  !     ASCII or NetCDF files.
  !
  !     Alberto Garcia, 2000, 2001
  !     E. Anglada 2008
  use chemical, only : read_chemical_types, species_label, number_of_species, atomic_number
  use sys, only: die
  use precision, only: dp
  use atom_types, only: set_number_of_species, get_number_of_species, hilbert_vector_t, &
       label_length, symbol_length, get_l_v, get_n_v, get_zeta_v, get_pop_v, get_energy_v, &
       get_pol_v, get_rad_func_p_v, get_read_from_file, get_label, get_symbol, get_atomic_number, &
       get_valence_charge, get_mass, get_self_energy, get_number_of_orbs, get_lmax_orbs, &
       get_lmax_kb_proj, get_number_of_orbs, get_number_of_kb_projs, get_orb_v, get_rad_func_v,&
       destroy_vector, get_number_of_orbs_non_deg, has_kbs, get_kb_v, get_number_of_kb_non_deg,&
       has_neutral_atom_potential, get_neutral_atom_potential, has_pseudo_local_charge, &
       get_pseudo_local_charge, has_reduced_vlocal, get_reduced_vlocal, has_core_charge,&
       get_core_charge, get_symbol, set_label, set_atomic_number, set_valence_charge, set_mass,&
       set_self_energy, set_lmax_orbs, set_lmax_kb_proj,init_orbs, init_kb_proj, init_vector,&
       set_symbol, set_orb, set_orbs_deg, set_kb_projs_deg, set_neutral_atom_potential,&
       set_pseudo_local_charge, set_core_charge, set_read_from_file, get_kb_proj_l, &
       get_kb_proj_n, get_kb_proj_energy, get_kb_proj_n, get_orb_l, get_orb_n, get_orb_zeta,&
       get_orb_pol, get_orb_pop, set_has_core_charge,set_orb,set_kb_proj

  use atom_generation_types, only: write_basis_specs, basis_parameters
  use pseudopotential_new, only : pseudopotential_new_t, pseudo_header_print
  use radial, only : rad_func_t, rad_grid_t, rad_read_ascii, rad_dealloc, rad_dump_ascii, &
       rad_dump_fft_file, rad_dump_xml, ntbmax, rad_cutoff, rad_get_grid, rad_grid_get_delta,&
       rad_grid_alloc, rad_dump_netcdf, rad_alloc, rad_grid_dealloc
  
  use fdf
  !      use flib_wxml, only: str
  use xml, only: xml_dump_attribute, xml_dump_element, str

  implicit none

  public :: dump_basis_ascii, read_basis_ascii
  public :: dump_basis_xml
  public :: dump_basis_netcdf, read_basis_netcdf
  public :: read_ion_ascii

  private

CONTAINS

  subroutine read_basis_netcdf

#ifndef CDF
    call die('*** You need netCDF to read the new user-defined basis files...')
  end subroutine  read_basis_netcdf
#else 

  use netcdf

  type(rad_func_t), pointer            :: op
  type(rad_func_t), pointer            :: pp

  integer :: ncid, iret, nspecies

  integer :: nkbs, nkbs_id, ntb_id, proj_id, pjnl_l_id, pjnl_n_id, &
       pjnl_ekb_id, kbdelta_id, kbcutoff_id
  integer :: norbs, norbs_id, orbnl_l_id, orbnl_n_id, orbnl_z_id, &
       cutoff_id, delta_id, orb_id, orbnl_pop_id, orbnl_ispol_id
  integer :: vna_id, chlocal_id, core_id, nprojs

  integer :: isp, j, i, l, nrp_tables, core_flag, nor, nk, m,itmp
  real(dp) :: tmp, cutoff, delta
  character(len=40) :: filename
  character(len=40) :: dummy
  character(len=symbol_length) :: symbol
  type(rad_func_t) :: rad_tmp
  type(rad_grid_t) :: grid_tmp

  real(dp), pointer :: r_tmp(:),e(:),pop(:)
  integer, pointer :: l_tmp(:),n_tmp(:),zeta_tmp(:),aux(:)
  logical, pointer :: pol(:)
  type(hilbert_vector_t) :: vector
  !
  call read_chemical_types
  nspecies = number_of_species()

  call set_number_of_species(nspecies)

  do isp = 1, nspecies

     call set_label(isp,species_label(isp))
     call set_read_from_file(isp,.true.)
     
     write(filename,'(a,a)') trim(get_label(isp)), ".ion.nc"
     iret = nf90_open(trim(filename),NF90_NOWRITE,ncid)

     iret = nf90_inq_dimid(ncid,'norbs',norbs_id)
     iret = nf90_inquire_dimension(ncid,norbs_id,len=norbs)
!     if (norbs .gt. maxnorbs) &
!          call die("read_user_basis: Increase maxnorbs in atm_types.f")
!
!     spp%n_orbnl = norbs
!
     iret = nf90_inq_dimid(ncid,'nkbs',nkbs_id)
     iret = nf90_inquire_dimension(ncid,nkbs_id,len=nkbs)
!     spp%n_pjnl = nkbs
!
     !
     !        For now, it is assumed that *all* the radial arrays have
     !        the same length.
     !
     iret = nf90_inq_dimid(ncid,'ntb',ntb_id)
     iret = nf90_inquire_dimension(ncid,ntb_id,len=nrp_tables)
     if (nrp_tables .ne. NTBMAX) call die("NTBMAX mismatch")
!
!     allocate(spp%orbnl(norbs))
!     allocate(spp%pjnl(nkbs))
!
!
     iret = nf90_get_att(ncid,nf90_global,'Element',symbol)
     call set_symbol(isp,symbol)
     iret = nf90_get_att(ncid,nf90_global,'Label',dummy)
!     !! Sanity check here??
!
     iret = nf90_get_att(ncid,nf90_global,'Atomic_number',itmp)
     if (atomic_number(isp) .ne. tmp) call die("Atomic number mismatch")
     call set_atomic_number(isp,itmp)
!
     iret = nf90_get_att(ncid,nf90_global,'Valence_charge',tmp)
     call set_valence_charge(isp,tmp)
     iret = nf90_get_att(ncid,nf90_global,'Mass',tmp)
     call set_mass(isp,tmp)

     iret = nf90_get_att(ncid,nf90_global,'Self_energy', tmp)
     call set_self_energy(isp,tmp)
     iret = nf90_get_att(ncid,nf90_global, 'Number_of_orbitals',norbs)
     iret = nf90_get_att(ncid,nf90_global, 'L_max_basis',itmp)
     call set_lmax_orbs(isp,itmp)

     iret = nf90_get_att(ncid,nf90_global, 'Number_of_projectors',nprojs)

     iret = nf90_get_att(ncid,nf90_global,  'L_max_projs',itmp)
     call set_lmax_kb_proj(isp,itmp)

!
!     !! Orbitals
!

     iret = nf90_inq_varid(ncid,'orbnl_l',orbnl_l_id)
     iret = nf90_inq_varid(ncid,'orbnl_n',orbnl_n_id)
     iret = nf90_inq_varid(ncid,'orbnl_z',orbnl_z_id)
     iret = nf90_inq_varid(ncid,'orbnl_ispol',orbnl_ispol_id)
     iret = nf90_inq_varid(ncid,'orbnl_pop',orbnl_pop_id)
!
     iret = nf90_inq_varid(ncid,'cutoff',cutoff_id)
     iret = nf90_inq_varid(ncid,'delta',delta_id)
!
!     !!      Projectors
!
     iret = nf90_inq_varid(ncid,'pjnl_l',pjnl_l_id)
     call check(iret)
     iret = nf90_inq_varid(ncid,'pjnl_n',pjnl_n_id)
     iret = nf90_inq_varid(ncid,'pjnl_ekb',pjnl_ekb_id)
     iret = nf90_inq_varid(ncid,'kbcutoff',kbcutoff_id)
     iret = nf90_inq_varid(ncid,'kbdelta',kbdelta_id)
     call check(iret)
!
     iret = nf90_inq_varid(ncid,'orb',orb_id)
     call check(iret)
!     !
!     !       Neutral atom potential
!     !
     iret = nf90_inq_varid(ncid,'vna',vna_id)
     iret = nf90_get_att(ncid,vna_id,'Vna_cutoff',cutoff)
     iret = nf90_get_att(ncid,vna_id,'Vna_delta',delta)

     grid_tmp = rad_grid_alloc(length=ntbmax,delta=delta)
     
     allocate(r_tmp(1:ntbmax))

     iret = nf90_get_var(ncid,vna_id,r_tmp(1:),start=(/1/),count=(/NTBMAX/))
     call check(iret)

     call rad_alloc(rad_tmp,r_tmp,grid_tmp)
     call set_neutral_atom_potential(isp,rad_tmp)
     call rad_dealloc(rad_tmp)
     call rad_grid_dealloc(grid_tmp)
     deallocate(r_tmp)

     !
     !       Local potential charge density
     !
     iret = nf90_inq_varid(ncid,'chlocal',chlocal_id)
     iret = nf90_get_att(ncid,chlocal_id, 'Chlocal_cutoff',cutoff)
     iret = nf90_get_att(ncid,chlocal_id, 'Chlocal_delta',delta)

     grid_tmp = rad_grid_alloc(length=ntbmax,delta=delta)
     allocate(r_tmp(1:ntbmax))

     iret = nf90_get_var(ncid,chlocal_id,r_tmp(1:),start=(/1/),count=(/NTBMAX/))
     call check(iret)

     call rad_alloc(rad_tmp,r_tmp,grid_tmp)
     call set_pseudo_local_charge(isp,rad_tmp)
     call rad_dealloc(rad_tmp)
     call rad_grid_dealloc(grid_tmp)
     deallocate(r_tmP)
     
     !
     !       Core charge
     !
     iret = nf90_get_att(ncid,nf90_global,'Core_flag',core_flag)
     if (core_flag .eq. 1) then
        call set_has_core_charge(isp,.true.)
     else
        call set_has_core_charge(isp,.false.)
     endif
!
     if (has_core_charge(isp)) then
        iret = nf90_inq_varid(ncid,'core',core_id)
        iret = nf90_get_att(ncid,core_id,'Core_cutoff',cutoff)
        iret = nf90_get_att(ncid,core_id,'Core_delta',delta)
        grid_tmp = rad_grid_alloc(length=ntbmax,delta=delta)
        allocate(r_tmp(1:ntbmax))
        iret = nf90_get_var(ncid,chlocal_id,r_tmp(1:),start=(/1/),count=(/NTBMAX/))
        call check(iret)

        call rad_alloc(rad_tmp,r_tmp,grid_tmp)
        call set_core_charge(isp,rad_tmp)
        call rad_dealloc(rad_tmp)
        call rad_grid_dealloc(grid_tmp)
        deallocate(r_tmP)
!     else
!        call rad_zero(spp%core)
     endif
!
     call check(iret)
     iret = nf90_inq_varid(ncid,'proj',proj_id)
     call check(iret)

     allocate(l_tmp(1:nkbs),n_tmp(1:nkbs),e(1:nkbs))
!
     iret = nf90_get_var(ncid,pjnl_l_id,l_tmp,count=(/nkbs/))
     call check(iret)
     iret = nf90_get_var(ncid,pjnl_n_id,n_tmp,count=(/nkbs/))
     call check(iret)
     iret = nf90_get_var(ncid,pjnl_ekb_id,e, count=(/nkbs/))
     call check(iret)

     do i = 1, nkbs
        iret = nf90_get_var(ncid,kbcutoff_id,cutoff,start=(/i/))
        call check(iret)
        iret = nf90_get_var(ncid,kbdelta_id,delta,start=(/i/))

        grid_tmp = rad_grid_alloc(length=ntbmax,delta=delta)
        allocate(r_tmp(1:ntbmax))
        iret = nf90_get_var(ncid,proj_id,r_tmp,start=(/1,i/),count=(/NTBMAX,1/))
        call check(iret)

        call rad_alloc(rad_tmp,r_tmp,grid_tmp)
        call init_vector(vector,rad_tmp,0,l_tmp(i),i,0.0_dp,e(i),.false.)
        call set_kb_proj(isp,vector,i)

        call rad_dealloc(rad_tmp)
        call rad_grid_dealloc(grid_tmp)
        call destroy_vector(vector)
     enddo
     call set_kb_projs_deg(isp)

     deallocate(l_tmp,n_tmp,e)

     allocate(l_tmp(1:norbs),n_tmp(1:norbs),zeta_tmp(1:norbs),pop(1:norbs),pol(1:norbs))

     iret=nf90_get_var(ncid,orbnl_l_id,l_tmp,count=(/norbs/))
     iret=nf90_get_var(ncid,orbnl_n_id,n_tmp,count=(/norbs/))
     iret=nf90_get_var(ncid,orbnl_z_id,zeta_tmp,count=(/norbs/))
!
     deallocate(aux)
     allocate(aux(1:norbs))
     iret = nf90_get_var(ncid,orbnl_ispol_id,aux,count=(/norbs/))
     do i = 1, norbs
        pol(i) =  aux(i) .eq. 1
     enddo
     call check(iret)
     iret = nf90_get_var(ncid,orbnl_pop_id,pop,count=(/norbs/))
     call check(iret)
!

     do i = 1, norbs
        iret = nf90_get_var(ncid,cutoff_id,cutoff,start=(/i/))
        call check(iret)
        iret = nf90_get_var(ncid,delta_id,delta,start=(/i/))

        grid_tmp = rad_grid_alloc(length=ntbmax,delta=delta)
        allocate(r_tmp(1:ntbmax))
        iret = nf90_get_var(ncid,orb_id,r_tmp,start=(/1,i/),count=(/NTBMAX,1/))
        call check(iret)

        call rad_alloc(rad_tmp,r_tmp,grid_tmp)
        call init_vector(vector,rad_tmp,n_tmp(i),l_tmp(i),zeta_tmp(i),pop(i),0.0_dp,pol(i))
        call set_orb(isp,vector,i)

        call rad_dealloc(rad_tmp)
        call rad_grid_dealloc(grid_tmp)
        call destroy_vector(vector)
     enddo
     call set_orbs_deg(isp)



!
!
!
!     nor = 0
!     do i = 1, norbs
!        op => spp%orbnl(i)
!        call rad_alloc(op,NTBMAX)
!        iret = nf90_get_var(ncid,orb_id,op%f(1:),start=(/1,i/),count=(/NTBMAX,1/))
!        call check(iret)
!        iret = nf90_get_var(ncid,cutoff_id,op%cutoff,start=(/i/))
!        call check(iret)
!        iret = nf90_get_var(ncid,delta_id,op%delta,start=(/i/))
!        call check(iret)
!        call rad_setup_d2(op)
!        l = spp%orbnl_l(i)
!        do m = -l,l
!           nor = nor+1
!           spp%orb_n(nor) = spp%orbnl_n(i)
!           spp%orb_l(nor) = spp%orbnl_l(i)
!           spp%orb_m(nor) = m
!           spp%orb_pop(nor) = spp%orbnl_pop(i) / (2*l+1)
!           spp%orb_index(nor) = i
!        enddo
!     enddo
!     spp%norbs = nor
!     iret = nf90_close(ncid)
!     call check(iret)
  enddo
  !
CONTAINS

  subroutine check(status)

    integer, intent(in):: status
    if (status .ne. nf90_noerr) then
       print  *, trim(nf90_strerror(status))
       call die()
    endif
  end subroutine check

end subroutine read_basis_netcdf
#endif
 
!=======================================================================

subroutine read_basis_ascii

  integer :: is, nspecies

  call read_chemical_types
  nspecies = get_number_of_species()

  call set_number_of_species(nspecies)

  do is = 1, nspecies
     call set_label(is, species_label(is))
     call set_read_from_file(is, .true.)
     call read_ion_ascii(is)
  enddo

end subroutine read_basis_ascii
!
!----------------------
subroutine read_ion_ascii(is)
  integer, intent(in) :: is

  character(len=20) filename
  integer:: i, l, lun, ispol,izeta, n
  real(dp) :: pop,ekb
  type(rad_func_t) :: rad_func
  type(hilbert_vector_t) :: orb

  write(filename,'(a,a)') trim(get_label(is)), ".ion"
  call io_assign(lun)
  open(lun,file=filename,status='old',form='formatted')
  rewind(lun)

  call read_header(is,lun)
  read(lun,*) 

  do i=1,get_number_of_orbs(is)
     read(lun,*) l, n, izeta,pop,ispol
     call rad_read_ascii(rad_func,lun)
     call init_vector(orb,rad_func,n,l,izeta,pop,0.0_dp,ispol.eq.1)
     call set_orb(is,orb,i)
     call destroy_vector(orb)
     call rad_dealloc(rad_func)
  enddo

  call set_orbs_deg(is) 

  ! KBs

  do i=1,get_number_of_orbs(is)
     read(lun,*) l, n, ekb
     call rad_read_ascii(rad_func,lun)
     call init_vector(orb,rad_func,n,l,1,0.0_dp,ekb,.false.)
     call set_orb(is,orb,i)
     call destroy_vector(orb)
     call rad_dealloc(rad_func)
  enddo

  call set_kb_projs_deg(is)
  !
  !Vna
  read(lun,*)
  call rad_read_ascii(rad_func,lun)
  call set_neutral_atom_potential(is,rad_func)
  call rad_dealloc(rad_func)

  !
  !Chlocal
  read(lun,*)
  call rad_read_ascii(rad_func,lun)
  call set_pseudo_local_charge(is,rad_func)
  call rad_dealloc(rad_func)

  !
  !Core charge
  read(lun,*,end=9999)
  call rad_read_ascii(rad_func,lun)
  call set_core_charge(is,rad_func)
  call rad_dealloc(rad_func)

9999 continue
  call io_close(lun)

CONTAINS

  subroutine read_header(is,unit)
    integer, intent(in)   :: is
    integer, intent(in)   :: unit

    character(len=78) :: line
    character(len=symbol_length) :: symbol
    character(len=label_length) :: label
    real(dp) :: zval,mass,self_energy
    integer :: lmax_basis, lmax_projs, norbs, nkbs, z

    read(unit,'(a)') line
    if (trim(line) .eq. '<preamble>') then
1      continue
       read(unit,'(a)') line
       if (trim(line) .ne. '</preamble>') goto 1
    endif

    read(unit,'(a2)') symbol
    read(unit,'(a20)') label
    read(unit,*) z
    read(unit,*) zval
    read(unit,*) mass
    read(unit,*) self_energy
    read(unit,*) lmax_basis, norbs
    read(unit,*) lmax_projs, nkbs
    
    call set_symbol(is,symbol)
    call set_label(is,label)
    call set_atomic_number(is,z)
    call set_valence_charge(is,zval)
    call set_mass(is,mass)
    call set_self_energy(is,self_energy)
    call set_lmax_orbs(is,lmax_basis)
    call set_lmax_kb_proj(is,lmax_projs)

    call init_orbs(is,norbs)
    call init_kb_proj(is,nkbs)

  end subroutine read_header

end subroutine read_ion_ascii


#ifndef CDF
subroutine dump_basis_netcdf
  !     Do nothing
end subroutine dump_basis_netcdf
#else
subroutine dump_basis_netcdf

  use netcdf

  type(rad_func_t), pointer            :: pp
  type(rad_func_t), pointer            :: op

  integer :: ncid, iret, nspecies

  integer :: nkbs, nkbs_id, ntb_id, proj_id, pjnl_l_id, pjnl_n_id,&
     pjnl_ekb_id, kbdelta_id, kbcutoff_id
  integer :: norbs, norbs_id, orbnl_l_id, orbnl_n_id, orbnl_z_id, &
     cutoff_id, delta_id, orb_id, orbnl_pop_id, orbnl_ispol_id
  integer :: vna_id, chlocal_id, reduced_vlocal_id, core_id

  !integer aux(maxnorbs)
  type(rad_func_t) :: rad_tmp
  type(rad_grid_t) :: grid_tmp
  type(hilbert_vector_t) :: vec_tmp
  real(dp) :: delta

  integer :: isp, j, i, l, nfuncs
  integer, pointer :: i_tmp(:)
  real(dp), pointer :: r_tmp(:)
  character*40 filename

  nspecies = get_number_of_species()
  do isp = 1, nspecies
!     spp => species(is)
     if (get_read_from_file(isp)) cycle  !! Do not dump


     write(filename,'(a,a)') trim(get_label(isp)), ".ion.nc"
     write(6,'(2a)') 'Dumping basis to NetCDF file ',trim(filename)

     iret = nf90_create(trim(filename),NF90_CLOBBER,ncid)

     nkbs = get_number_of_kb_non_deg(isp) !spp%n_pjnl
     norbs = get_number_of_orbs_non_deg(isp) !spp%n_orbnl

     iret = nf90_def_dim(ncid,'norbs',norbs,norbs_id)
     call check(iret)
     iret = nf90_def_dim(ncid,'nkbs',nkbs,nkbs_id)
     call check(iret)
     iret = nf90_def_dim(ncid,'ntb',NTBMAX,ntb_id)
     call check(iret)

!     !! Orbitals

     iret = nf90_put_att(ncid,nf90_global,'Element',get_symbol(isp))
     iret = nf90_put_att(ncid,nf90_global,'Label',get_label(isp))
     iret = nf90_put_att(ncid,nf90_global,'Atomic_number',get_atomic_number(isp))
     iret = nf90_put_att(ncid,nf90_global,'Valence_charge',get_valence_charge(isp))
     iret = nf90_put_att(ncid,nf90_global,'Mass',get_mass(isp))
     iret = nf90_put_att(ncid,nf90_global,'Self_energy',get_self_energy(isp))
     iret = nf90_put_att(ncid,nf90_global,'Number_of_orbitals',get_number_of_orbs(isp))
     iret = nf90_put_att(ncid,nf90_global,'L_max_basis',get_lmax_orbs(isp))
     iret = nf90_put_att(ncid,nf90_global,'Number_of_projectors',&
          get_number_of_kb_projs(isp))
     iret = nf90_put_att(ncid,nf90_global,'L_max_projs',get_number_of_kb_projs(isp))

     iret = nf90_def_var(ncid,'orbnl_l',nf90_int,norbs_id,orbnl_l_id)
     iret = nf90_def_var(ncid,'orbnl_n',nf90_int,norbs_id,orbnl_n_id)
     iret = nf90_def_var(ncid,'orbnl_z',nf90_int,norbs_id,orbnl_z_id)
     iret = nf90_def_var(ncid,'orbnl_ispol',nf90_int,norbs_id,orbnl_ispol_id)
     iret = nf90_def_var(ncid,'orbnl_pop',nf90_double,norbs_id,orbnl_pop_id)
!
     iret = nf90_def_var(ncid,'cutoff',nf90_double, norbs_id,cutoff_id)
     iret = nf90_def_var(ncid,'delta',nf90_double, norbs_id,delta_id)
!
!     !!      Projectors
!
     iret = nf90_def_var(ncid,'pjnl_l',nf90_int,nkbs_id,pjnl_l_id)
     call check(iret)
     iret = nf90_def_var(ncid,'pjnl_n',nf90_int,nkbs_id,pjnl_n_id)
     iret = nf90_def_var(ncid,'pjnl_ekb',nf90_double, nkbs_id,pjnl_ekb_id)
     iret = nf90_def_var(ncid,'kbcutoff',nf90_double, nkbs_id,kbcutoff_id)
     iret = nf90_def_var(ncid,'kbdelta',nf90_double, nkbs_id,kbdelta_id)
     call check(iret)
!
     iret = nf90_def_var(ncid,'orb',nf90_double, (/ntb_id,norbs_id/),orb_id)
     call check(iret)
     !
     !       Local potential
     !
     iret = nf90_def_var(ncid,'vna',nf90_double, (/ntb_id/),vna_id)
     rad_tmp = get_neutral_atom_potential(isp)
     grid_tmp = rad_get_grid(rad_tmp)
     iret = nf90_put_att(ncid,vna_id, 'Vna_cutoff',rad_cutoff(rad_tmp))
     delta = rad_grid_get_delta(grid_tmp)
     iret = nf90_put_att(ncid,vna_id,'Vna_delta',delta)
     call check(iret)
!!AG     call rad_dump_netcdf(rad_tmp,ncid,vna_id)
     call rad_dealloc(rad_tmp)
     call rad_grid_dealloc(grid_tmp)

     !
     !       Local potential charge density
     !
     iret = nf90_def_var(ncid,'chlocal',nf90_double, (/ntb_id/),chlocal_id)
     rad_tmp = get_pseudo_local_charge(isp)
     iret = nf90_put_att(ncid,chlocal_id, 'Chlocal_cutoff',rad_cutoff(rad_tmp))
     grid_tmp = rad_get_grid(rad_tmp)
     delta = rad_grid_get_delta(grid_tmp)
     iret = nf90_put_att(ncid,chlocal_id, 'Chlocal_delta',delta)
     call check(iret)
!!AG     call rad_dump_netcdf(rad_tmp,ncid,chlocal_id)
     call rad_dealloc(rad_tmp)
     call rad_grid_dealloc(grid_tmp)
     !
     !       Reduced Local potential (rV+2*Zval)
     !
     iret = nf90_def_var(ncid,'reduced_vlocal',nf90_double, (/ntb_id/),reduced_vlocal_id)
     rad_tmp = get_reduced_vlocal(isp)
     iret = nf90_put_att(ncid,reduced_vlocal_id, 'Reduced_vlocal_cutoff', rad_cutoff(rad_tmp))
     grid_tmp = rad_get_grid(rad_tmp)
     delta = rad_grid_get_delta(grid_tmp)
     iret = nf90_put_att(ncid,reduced_vlocal_id,'Reduced_vlocal_delta', delta)
     call check(iret)
!!AG     call rad_dump_netcdf(rad_tmp,ncid,reduced_vlocal_id)
     call rad_dealloc(rad_tmp)
     call rad_grid_dealloc(grid_tmp)

      !
      !       Core charge
      !
      if (has_core_charge(isp)) then
         iret = nf90_put_att(ncid,nf90_global, 'Core_flag',1)
         iret = nf90_def_var(ncid,'core',nf90_double,(/ntb_id/),core_id)
         rad_tmp = get_core_charge(isp)
         grid_tmp = rad_get_grid(rad_tmp)
         delta = rad_grid_get_delta(grid_tmp)
         iret = nf90_put_att(ncid,core_id,'Core_cutoff',rad_cutoff(rad_tmp))
         iret = nf90_put_att(ncid,core_id,'Core_delta',delta)
         call check(iret)

!!AG         call rad_dump_netcdf(rad_tmp,ncid,core_id)
         call rad_dealloc(rad_tmp)
         call rad_grid_dealloc(grid_tmp)
      else
         iret = nf90_put_att(ncid,nf90_global,'Core_flag',0)
      endif

      call check(iret)
      iret = nf90_def_var(ncid,'proj',nf90_double, (/ntb_id,nkbs_id/),proj_id)
      call check(iret)
! !
!!AG  End of definition phase.

      iret = nf90_enddef(ncid)
      call check(iret)
!!AG
      rad_tmp = get_neutral_atom_potential(isp)
      call rad_dump_netcdf(rad_tmp,ncid,vna_id)
      call rad_dealloc(rad_tmp)

      rad_tmp = get_pseudo_local_charge(isp)
      call rad_dump_netcdf(rad_tmp,ncid,chlocal_id)
      call rad_dealloc(rad_tmp)
      
      rad_tmp = get_reduced_vlocal(isp)
      call rad_dump_netcdf(rad_tmp,ncid,reduced_vlocal_id)
      call rad_dealloc(rad_tmp)

      if (has_core_charge(isp)) then
         rad_tmp = get_core_charge(isp)
         call rad_dump_netcdf(rad_tmp,ncid,core_id)
         call rad_dealloc(rad_tmp)
      endif

!!AG
      nfuncs = get_number_of_kb_non_deg(isp)
      allocate (i_tmp(1:nfuncs))
      do i=1,nfuncs
         i_tmp(i)=get_kb_proj_l(isp,i)
      enddo
      iret = nf90_put_var(ncid,pjnl_l_id,i_tmp,count=(/nfuncs/))
      call check(iret)

      do i=1,nfuncs
         i_tmp(i)=get_kb_proj_n(isp,i)
      enddo
      iret = nf90_put_var(ncid,pjnl_n_id,i_tmp,count=(/nfuncs/))
      call check(iret)

      allocate(r_tmp(1:nfuncs))
      do i=1,nfuncs
         r_tmp(i)=get_kb_proj_energy(isp,i)
      enddo
      iret = nf90_put_var(ncid,pjnl_ekb_id,r_tmp,count=(/nfuncs/))
      call check(iret)

      do i=1,nfuncs
         vec_tmp = get_kb_v(isp,i)
         rad_tmp = get_rad_func_v(vec_tmp)
         call rad_dump_netcdf(rad_tmp,ncid,proj_id)

         iret = nf90_put_var(ncid,kbcutoff_id,rad_cutoff(rad_tmp),start=(/i/))
         call check(iret)

         grid_tmp = rad_get_grid(rad_tmp)
         delta = rad_grid_get_delta(grid_tmp)
         iret = nf90_put_var(ncid,kbdelta_id,delta,start=(/i/))
         call check(iret)

         call rad_dealloc(rad_tmp)
         call destroy_vector(vec_tmp)
         call rad_grid_dealloc(grid_tmp)
      end do
      deallocate(i_tmp,r_tmp)

      nfuncs = get_number_of_orbs_non_deg(isp)

      allocate(i_tmp(1:nfuncs))
      do i=1,nfuncs
         i_tmp(i)=get_orb_l(isp,i)
      enddo
      iret = nf90_put_var(ncid,orbnl_l_id,i_tmp,count=(/norbs/))

      i_tmp(1:nfuncs) = 0
      do i=1,nfuncs
         i_tmp(i)=get_orb_n(isp,i)
      enddo
      iret = nf90_put_var(ncid,orbnl_n_id,i_tmp,count=(/norbs/))

      i_tmp(1:nfuncs) = 0
      do i=1,nfuncs
         i_tmp(i)=get_orb_zeta(isp,i)
      enddo
      iret = nf90_put_var(ncid,orbnl_z_id,i_tmp,count=(/norbs/))

!      if (norbs .gt. maxnorbs) call die("dump_basis_netcdf: Increase maxnorbs")

      i_tmp(1:nfuncs) = 0
      do i = 1, nfuncs
         if (get_orb_pol(isp,i)) i_tmp(i)=1
      enddo
      iret = nf90_put_var(ncid,orbnl_ispol_id,i_tmp,count=(/norbs/))
      call check(iret)

      allocate(r_tmp(1:nfuncs))
      r_tmp = 0.0_dp
      do i=1,nfuncs
         r_tmp(i)=get_orb_pop(isp,i)
      enddo
      iret = nf90_put_var(ncid,orbnl_pop_id,r_tmp,count=(/norbs/))
      call check(iret)


      do i=1,nfuncs
         vec_tmp = get_orb_v(isp,i)
         rad_tmp = get_rad_func_v(vec_tmp)
         call rad_dump_netcdf(rad_tmp,ncid,orb_id)

         iret = nf90_put_var(ncid,cutoff_id,rad_cutoff(rad_tmp),start=(/i/))
         call check(iret)

         grid_tmp = rad_get_grid(rad_tmp)
         delta = rad_grid_get_delta(grid_tmp)
         iret = nf90_put_var(ncid,delta_id,delta,start=(/i/))
         call check(iret)

         call rad_dealloc(rad_tmp)
         call destroy_vector(vec_tmp)
         call rad_grid_dealloc(grid_tmp)
      end do
      deallocate(i_tmp,r_tmp)

  enddo

contains
  subroutine check(status)

    integer, intent(in):: status
    if (status .ne. nf90_noerr) then
       print  *, trim(nf90_strerror(status))
       call pxfabort()
    endif
  end subroutine check

end subroutine dump_basis_netcdf

#endif

!------------------------------------

subroutine dump_basis_ascii

  integer :: is, nspecies

  nspecies = get_number_of_species()

  do is = 1, nspecies
     call dump_ion_ascii(is)
  enddo
end subroutine dump_basis_ascii

!------------------------------------

subroutine dump_ion_ascii(is)
  integer, intent(in)      :: is

  character*40 filename
  character*30 fileid
  integer i, lun, lun2, n_series
  type(hilbert_vector_t) :: vector
  type(rad_func_t) :: rfunc

  if (get_read_from_file(is)) return !! Do not dump

  write(filename,'(a,a)') trim(get_label(is)), ".ion"
  call io_assign(lun)
  open(lun,file=filename,status='replace',form='formatted')

  write(lun,'(a)') '<preamble>'
  call write_basis_specs(lun,is)
  call pseudo_header_print(lun,basis_parameters(is)%pseudopotential)
  write(lun,'(a)') '</preamble>'
  call write_header(is,lun)
  write(lun,'(a)') "# PAOs:__________________________"
  n_series = 0
  do i=1,get_number_of_orbs_non_deg(is)
     vector = get_orb_v(is,i)
     call dump_vector(vector,"orb",trim(get_label(is)),"ascii",lun)
     call destroy_vector(vector)
  enddo

  write(lun,'(a)') "# KBs:__________________________"
  if (has_kbs(is)) then
     n_series = 0
     do i=1,get_number_of_kb_non_deg(is)
        vector = get_kb_v(is,i)
        call dump_vector(vector,"kb",trim(get_label(is)),"ascii",lun)
        call destroy_vector(vector)
     enddo
  endif

  if(has_neutral_atom_potential(is)) then
     write(lun,'(a)') "# Vna:__________________________"
     write(fileid,'(a,a)') trim(get_label(is)),"-VNA.dat"
     call io_assign(lun2)
     open(unit=lun2,file=fileid,status='replace', form='formatted')
     write(lun2,'(3a)') "# ",trim(get_label(is)), " Vna"
     rfunc=get_neutral_atom_potential(is)
     call rad_dump_ascii(rfunc,lun2,header=.false.)
     call io_close(lun2)
     call rad_dump_ascii(rfunc,lun)
     !write(fileid,'(a,a)') trim(get_label(is)),"-VNA-fft.dat"
     !call rad_dump_fft_file(rfunc,fileid)
     call rad_dealloc(rfunc)
  endif

  if (has_pseudo_local_charge(is))then
     write(lun,'(a)') "# Chlocal:__________________________"
     write(fileid,'(a,a)') trim(get_label(is)),"-CHLOCAL.dat"
     call io_assign(lun2)
     open(unit=lun2,file=fileid,status='replace', form='formatted')
     write(lun2,'(3a)') "# ",trim(get_label(is)), " ChLocal"
     rfunc = get_pseudo_local_charge(is)
     call rad_dump_ascii(rfunc,lun2,header=.false.)
     call io_close(lun2)
     call rad_dump_ascii(rfunc,lun)
     call rad_dealloc(rfunc)
  endif

  if (has_reduced_vlocal(is))then
  !
  !        Vlocal does not go to the .ion file
  !
     write(fileid,'(a,a)') trim(get_label(is)),"-RED_VLOCAL.dat" 
     call io_assign(lun2)
     open(unit=lun2,file=fileid,status='replace',  form='formatted')
     write(lun2,'(3a)') "# ",trim(get_label(is))," Red_Vlocal"
     rfunc = get_reduced_vlocal(is)
     call rad_dump_ascii(rfunc,lun2,header=.false.)
     call rad_dealloc(rfunc)
     call io_close(lun2)
  endif

  !
  if (has_core_charge(is)) then
     write(lun,'(a)') "# Core:__________________________"
     write(fileid,'(a,a)') trim(get_label(is)),"-CHCORE.dat"
     call io_assign(lun2)
     open(unit=lun2,file=fileid,status='replace', form='formatted')
     write(lun2,'(3a)') "# ",trim(get_label(is)), " ChCore"
     rfunc = get_core_charge(is)
     call rad_dump_ascii(rfunc,lun2,header=.false.)
     call io_close(lun2)
     call rad_dump_ascii(rfunc,lun)
     !FFT
     !write(fileid,'(a,a)') trim(get_label(is)),"-CHCORE-fft.dat"
     !call io_assign(lun2)
     !open(unit=lun2,file=fileid,status='replace', form='formatted')
     !write(lun2,'(3a)') "# ",trim(get_label(is)), " ChCore"
     !call rad_dump_fft_ascii(rfunc,lun2)
     !call rad_dump_fft_file(rfunc,fileid) 
     !call io_close(lun2)
     
     call rad_dealloc(rfunc)
  endif
  call io_close(lun)

CONTAINS

  subroutine write_header(is,unit)
    integer, intent(in) :: is
    integer, intent(in) :: unit

    write(unit,'(a2,28x,a)') get_symbol(is), "# Symbol"
    write(unit,'(a20,10x,a)') get_label(is), "# Label"
    write(unit,'(i5,25x,a)') get_atomic_number(is), "# Atomic number"
    write(unit,'(g22.12,25x,a)') get_valence_charge(is),   "# Valence charge"
    write(unit,'(g22.12,4x,a)') get_mass(is), "# Mass "
    write(unit,'(g22.12,4x,a)') get_self_energy(is),      "# Self energy "
    write(unit,'(2i4,22x,a)') get_lmax_orbs(is), get_number_of_orbs_non_deg(is), &
         "# Lmax for basis, no. of nl orbitals "
    write(unit,'(2i4,22x,a)') get_lmax_kb_proj(is),get_number_of_kb_non_deg(is),&
         "# Lmax for projectors, no. of nl KB projectors "

  end subroutine write_header

end subroutine dump_ion_ascii
!
!-----------------------------------------------------------------
subroutine dump_basis_xml

  integer :: is, nspecies
  nspecies = get_number_of_species()

  do is = 1, nspecies
     call dump_ion_xml(is)
  enddo
end subroutine dump_basis_xml
!---
subroutine dump_ion_xml(is)

  integer, intent(in)      :: is

  type(hilbert_vector_t)   :: vector
  type(rad_func_t)         :: rfunc
  character*40 filename
  integer i, lun

  if (get_read_from_file(is)) return    !! Do not dump

  write(filename,'(a,a)') trim(get_label(is)), ".ion.xml"
  call io_assign(lun)
  open(lun,file=filename,status='replace',form='formatted')

  write(lun,'(a)') '<ion version="0.1">'
  call xml_dump_element(lun,'symbol',str(get_symbol(is)))
  call xml_dump_element(lun,'label',str(get_label(is)))
  call xml_dump_element(lun,'z',str(get_atomic_number(is)))
  call xml_dump_element(lun,'valence',str(get_valence_charge(is)))
  call xml_dump_element(lun,'mass',str(get_mass(is)))
  call xml_dump_element(lun,'self_energy',str(get_self_energy(is)))
  call xml_dump_element(lun,'lmax_basis',str(get_lmax_orbs(is)))
  call xml_dump_element(lun,'norbs_nl',str(get_number_of_orbs(is)))
  call xml_dump_element(lun,'lmax_projs',str(get_lmax_kb_proj(is))) 
  call xml_dump_element(lun,'nprojs_nl',str(get_number_of_kb_projs(is)))

  write(lun,'(a)') '<preamble>'
  call write_basis_specs(lun,is)
  call pseudo_header_print(lun,basis_parameters(is)%pseudopotential)
  write(lun,'(a)') '</preamble>'

  write(lun,'(a)') "<paos>"
  do i=1,get_number_of_orbs_non_deg(is)
     vector = get_orb_v(is,i)
     call dump_vector(vector,"orb",trim(get_label(is)),"xml",lun)
     call destroy_vector(vector)
  enddo
  write(lun,'(a)') "</paos>"

  write(lun,'(a)') "<kbs>"
  if(has_kbs(is))then
     do i=1,get_number_of_kb_non_deg(is)
        vector=get_kb_v(is,i)
        call dump_vector(vector,"kb",trim(get_label(is)),"xml",lun)
        call destroy_vector(vector)
     enddo
  else
     write(lun,'(a)') "No kbs."
  endif
  write(lun,'(a)') "</kbs>"

  write(lun,'(a)') "<vna>"
  if (has_neutral_atom_potential(is)) then
     rfunc=get_neutral_atom_potential(is)
     call rad_dump_xml(rfunc,lun)
     call rad_dealloc(rfunc)
  endif
  write(lun,'(a)') "</vna>"

  write(lun,'(a)') "<chlocal>"
  if (has_pseudo_local_charge(is)) then
     rfunc=get_pseudo_local_charge(is)
     call rad_dump_xml(rfunc,lun)
     call rad_dealloc(rfunc)
  else
     write(lun,'(a)') "No chlocal"
  endif
  write(lun,'(a)') "</chlocal>"

  write(lun,'(a)') "<reduced_vlocal>"
  if(has_reduced_vlocal(is)) then
     rfunc=get_reduced_vlocal(is)
     call rad_dump_xml(rfunc,lun)
     call rad_dealloc(rfunc)
  endif
  write(lun,'(a)') "</reduced_vlocal>"

  if (has_core_charge(is)) then
     write(lun,'(a)') "<core>"
     rfunc=get_core_charge(is)
     call rad_dump_xml(rfunc,lun)
     call rad_dealloc(rfunc)
     write(lun,'(a)') "</core>"
  endif

  write(lun,'(a)') "</ion>"

  call io_close(lun)

end subroutine dump_ion_xml

!---------------------------------------------------------------------

subroutine dump_vector(vector,kind,label,output_kind,lun)
  type(hilbert_vector_t), intent(in) :: vector
  character(len=*), intent(in)     :: kind,output_kind,label
  integer, intent(in)              :: lun
 
  type(rad_func_t), pointer        :: func 
  
  integer                          :: n,l,zeta,ispol=0
  integer                          :: lun2
  real(dp)                         :: pop,energy
  character(len=60)                :: filename

  l=get_l_v(vector)
  n=get_n_v(vector)
  zeta=get_zeta_v(vector)
  pop=get_pop_v(vector)
  energy=get_energy_v(vector)

  if( get_pol_v(vector) ) ispol=1

  if (output_kind == "ascii") then    
     if (kind == "orb") then
        !lun2: individual files
        write(filename,'(a,a,a,a,i1,a,i1,a,i1,a)') trim(label),"-",&
             trim(kind),"-n=",n,"-l=",l,"-zeta=",zeta,".dat"

        call io_assign(lun2)
        open(unit=lun2,file=filename,status='replace',form='formatted')
        func => get_rad_func_p_v(vector)
        call rad_dump_ascii(func,lun2,header=.false.)
        call io_close(lun2)

        write(lun,'(4i3,f10.6,2x,a)') l,n,zeta,ispol,pop, &
             " #orbital l, n, z, is_polarized, population"
        call rad_dump_ascii(func,lun,header=.true.)

        !FFT
        !write(filename,'(a,a,a,a,i1,a,i1,a,i1,a)') trim(label),"-",&
        !     trim(kind),"-n=",n,"-l=",l,"-zeta=",zeta,"-fft.dat"
        !call io_assign(lun2)
        !open(unit=lun2,file=filename,status='replace',form='formatted')
        !write(lun2,'(a)') "#(species label, l, n, z, is_polarized, popul)"
        !call rad_dump_fft_file(func,filename)
        !call io_close(lun2)

     elseif(kind=="kb")then
         !lun2: individual files
        write(filename,'(a,a,a,a,i1,a,i1,a)') trim(label),"-",&
             trim(kind),"-n=",n,"-l=",l,".dat"

        call io_assign(lun2)
        open(unit=lun2,file=filename,status='replace',form='formatted')
        func => get_rad_func_p_v(vector)
        !write(lun,'(2i3,f10.6,2x,a)') l,n,energy, &
        !     " #kb l, n, Reference energy"
        call rad_dump_ascii(func,lun2,header=.true.)
        call io_close(lun2)

        write(lun,'(2i3,4x,f10.6,2x,a)') l,n,energy, &
             " #kb l, n, Reference energy"
        call rad_dump_ascii(func,lun,header=.true.)
     endif
  elseif(output_kind == "xml")then
     if(kind=="orb")then
        write(lun,'(a)') "<orbital "
        call xml_dump_attribute(lun,'l',str(l))
        call xml_dump_attribute(lun,'n',str(n))
        call xml_dump_attribute(lun,'z',str(zeta))
        call xml_dump_attribute(lun,'ispol',str(ispol))
        call xml_dump_attribute(lun,'population',str(pop))
        write(lun,'(a)') " >"
        call rad_dump_xml(get_rad_func_v(vector),lun)
        write(lun,'(a)') "</orbitals>"
     elseif(kind=="kb")then
         write(lun,'(a)') "<projector "
         call xml_dump_attribute(lun,'l',str(l))
         call xml_dump_attribute(lun,'n',str(n))
         call xml_dump_attribute(lun,'ref_energy',str(energy))
         write(lun,'(a)') " >"
         call rad_dump_xml(get_rad_func_v(vector),lun)
         write(lun,'(a)') "</projector>"
     else
        call die("basis_io: dump_xml: unknown function!")
     end if
endif
end subroutine dump_vector

end module atom_basis_io












