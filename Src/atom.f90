module atom
  !Common modules
  use precision, only: dp
  use sys, only: die
  use fdf
  use units, only: eV
  
  !Generic module for radial (1D) functions
  use radial
  
  !Main module where all the atomic info is stored.
  use atm_types, only : get_atomic_number, get_symbol,species_info_t, &
       species,nspecies, symbol_length, set_no_reduced_vlocal, set_no_kb, &
       set_no_neutral_atom_potential
  
  !Auxiliary module where all the intermediate info is stored.
  use atom_generation_types, only : basis_def_t, basis_parameters,lshell_t, shell_t, &
       write_basis_specs
  
  !All the auxiliary modules
  use atmparams
  use periodic_table, only : symbol
  use pseudopotential, only : pseudopotential_t
  
  !All the experte subroutines which generate the corresponding functions.
  use corechg, only: coreChargeSetup
  
  use bessel_m, only : bessel
  use parallel, only:ionode
  use na, only: gen_vna
  
  implicit none      

  private

  !All the atomic info is generated when this subroutine is called.
  public:: generate_all_atomic_info 

CONTAINS

  subroutine atom_main(isp)
	!Main, expert subroutine which calls the all the subroutines
	!which generate ALL the atomic info.
	
	!Arguments
	!integer, intent(in) :: isp !Species index. All 

    use kb,    only: kbgen
    use vloc,  only: gen_vlocal 
    use basis, only: basis_gen

    integer , intent(in) :: isp

    write(6,'(/,a)') 'ATOM: Species begin__________________________ '

    call species_init(isp)

    if (basis_parameters(isp)%bessel) then
       call set_no_reduced_vlocal(species(isp))
       call set_no_kb(species(isp))
       call set_no_neutral_atom_potential(species(isp))
       call bessel(isp)
    elseif (basis_parameters(isp)%floating) then
       call species_pseudopotential_init(isp)
       call species_charge(isp)
       call set_no_reduced_vlocal(species(isp))
       call set_no_kb(species(isp))
       call set_no_neutral_atom_potential(species(isp))
       call basis_gen(isp)
    else
       call species_pseudopotential_init(isp)       
       call species_charge(isp)
       call gen_vlocal(isp)
       call kbgen(isp)
       call basis_gen(isp)
       call gen_vna(isp)
    endif

    write(6,'(/,a)') 'ATOM: Species end_____________________________ '

  end subroutine atom_main

  !----------------------------------------------------

  subroutine calculate_rho_potential(charge,rho,potential,irel,rho_core)
    real(dp), intent(in)                   :: charge
    type(rad_func_t), intent(inout)        :: rho
    type(rad_func_t), intent(out)          :: potential
    integer, intent(in)                    :: irel
    type(rad_func_t), optional, intent(in) :: rho_core

    type(rad_func_t)                :: rho4pir2!rho/4pir2
    type(rad_func_t)                :: rho_core4pir2
    type(rad_func_t)                :: rho_scaled !rho multiplied by charge
    type(rad_func_t)                :: pot_xc     !XC potential due to rho
    type(rad_func_t)                :: pot_hartree !Hartree pot. due to rho
    type(rad_func_t)                :: rho_tot    !rho/4pir2+core_charge
    real(dp)                        :: rmax, exc

    rmax = rad_cutoff(rho)

    !------Hartree-potential-----------------
    rho_scaled = rad_multiply_each_value(rho,charge)
    pot_hartree = rad_vhartree(rho_scaled)
    !call rad_dump_file(pot_hartree,"pot_hartree.dat")

    ! NOTE that atomxc expects true rho(r), not 4pir^2*rho(r)
    rho4pir2=rad_divide_by_4pir2(rho_scaled,.false.)


    !------XC-potential----------------------
    !If exists add the core charge due to nonlinear core correcctions
    call rad_dealloc(rho)
    if (present(rho_core)) then
       rho_core4pir2 = rad_divide_by_4pir2(rho_core,.false.)
       !call rad_dump_file(rho_core4pir2,"rho_core4pir2.dat")
       rho_tot = rad_sum(rho4pir2,rho_core4pir2,rmax)
       !call rad_dump_file(rho_tot,"rho_tot.dat")
       call rad_dealloc(rho_core4pir2)
    else
       call rad_copy(rho4pir2,rho_tot)
    endif

    !    Compute xc energy and potential and add to ve
    pot_xc = rad_vxc(rho_tot,irel,exc)
    !call rad_dump_file(pot_xc,"pot_xc.dat")

    !Update total potential
    potential =  rad_sum(pot_hartree,pot_xc)
    !Store rho
    call rad_copy(rho_scaled,rho)

    !Memory release
    call rad_dealloc(rho_scaled)
    call rad_dealloc(pot_xc)
    call rad_dealloc(pot_hartree)
    call rad_dealloc(rho4pir2)
    call rad_dealloc(rho_tot)
  endsubroutine calculate_rho_potential

  !----------------------------------------

  subroutine generate_all_atomic_info()

    ! Routine to initialize the Pseudopotentials and Atomic Orbitals.
    ! Substantially modified by Alberto Garcia (2000)
    !
    ! The PAO and KB information can optionally be read from ASCII files
    ! (those produced by a standard run of Siesta or Base, but with the extension 
    ! renamed to '.input' instead of '.dump'), or from NetCDF files (if NetCDF
    ! is available). Note that there must be files for *all* species.
    !
    ! This behavior is controlled by the FDF logicals 'user-basis' and
    ! 'user-basis-netcdf'.
    !
    ! The old 'USER' basis type has been removed.
    !
    ! This routine also outputs information about the basis specification
    ! determined by the routines in the 'basis_specs' modules. 
    !

    use basis_specs, only: read_basis_specs
    use basis_io, only: read_basis_ascii, dump_basis_ascii,dump_basis_xml, dump_basis_netcdf, read_basis_netcdf
    use elec_correction, only: elec_corr_setup
    !use atmfuncs, only: check_atmfuncs

    implicit none

    ! Internal variables ...................................................
    integer is

    logical user_basis, user_basis_netcdf, read_from_file

    write(6,'(/2a)')'initatom: Reading input for the', &
         'pseudopotentials and atomic orbitals'

    user_basis = fdf_boolean('user-basis',.false.)
    user_basis_netcdf = fdf_boolean('user-basis-netcdf',.false.)

    read_from_file = user_basis .or. user_basis_netcdf

    if (user_basis_netcdf) then

       write(6,'(a)') 'Reading PAOs and KBs from NetCDF files...'
       call read_basis_netcdf()

    else if (user_basis) then

       write(6,'(a)') 'Reading PAOs and KBs from ascii files...'
       call read_basis_ascii()

    else

       !     Read basis specification and generate PAOs and KB projectors
       !     Note : In this section, the number of species is nspecies
       call read_basis_specs()

       allocate(species(1:nspecies))

       do is = 1,nspecies
          call write_basis_specs(6,is)
          call atom_main(is) 
       enddo

       !     Compute the electrostatic correction tables
       call elec_corr_setup()

    endif

    !     Dump all the atomic info
    call dump_basis_ascii()

    call dump_basis_netcdf()

    call dump_basis_xml()

    !call check_atmfuncs()

  end subroutine generate_all_atomic_info

  !----------------------------------------------------

  subroutine species_charge(isp)
    use atm_types, only : get_valence_charge
    integer, intent(in) :: isp

    !Internal vars
    real(dp)                         :: charge 
    type(pseudopotential_t),pointer  :: vps 
    type(species_info_t), pointer    :: spp
    type(basis_def_t), pointer       :: basp
    integer                          :: irel
    type(rad_func_t), pointer        :: chcore
    real(dp)                         :: rescale_charge
    character(len=3)                 :: xcfunc    !Xc used
    character(len=4)                 :: xcauth    !Authors of xc


    basp => basis_parameters(isp)
    spp => species(isp)
    vps => basp%pseudopotential


    if (abs(get_atomic_number(spp)-vps%gen_zval).gt.1.0d-3) then 
       write(6,'(/,a,/,a,f5.2)') &
            'ATOM: Pseudopotential generated from an ionic configuration',&
            'ATOM: with net charge', get_valence_charge(spp)-vps%gen_zval
    endif

    !
    !           AG: Note zval-chgvps = (Znuc-Zcore)-gen_zval
    !           Example: Ba with 5s and 5p in the valence:
    !           zval from atom (scaled with zratio): 10 (5s2 5p6 6s2)
    !           chgvps (only 5s2 and 5p6 occupied in ref config): 8
    !           gen_zval = chgvps
    !           ==> net charge (zval-chgvps) = 2
    !           Znuc = 56,  Zcore= Z(Xe)-2-6 = 54-2-6 = 46
    !           Znuc-Zcore-true_val = 56 - 46 - 8 = 2
    !
    !           If we stop the practice of scaling the valence charge
    !           in the ps file, we would need info about Zcore to setup
    !           things correctly (both for Siesta and for the PW program...)
    !           BUT actually, we have that information in the ps file!!
    !           The potentials are stored as r*V, so at large r they all
    !           go to -2*(Znuc-Zcore).....

    !    
    !    Set 'charge':
    !    1. If 'basp%ionic_charge' is zero (that is, not set in the fdf file)
    !    then set it to basp%zval-vps%chg.
    !    2. If 'charge' is equal to basp%zval-vps%chg, set it to that.
    !    
    charge = basp%ionic_charge
    if ( abs(charge) == 0.d0  .or. &
         (abs(charge-vps%zval+vps%gen_zval) .lt. 1.0d-3) ) then   
       ! We can't overwrite basp%ionic_charge because it may be used 
       ! if there are semicore states 
       charge = get_valence_charge(spp) - vps%gen_zval
    endif

    !    Fill species core density used for the nonlinear core corrections
    call coreChargeSetup(isp)

    !    Relativistic pseudo
    if (vps%irel.eq.'rel') irel=1
    if (vps%irel.ne.'rel') irel=0

    xcfunc = fdf_string('xc.functional','LDA')
    xcauth = fdf_string('xc.authors','PZ')
    call xc_check(xcfunc,xcauth,vps%icorr)

    !-----------------------Rho, hartree potential for KBs-------------------
    !    CALCULATION OF THE VALENCE SCREENING POTENTIAL FROM THE READ CHARGE
    !    DENSITY
    !    
    !    For Kleinman-Bylander projectors calculation, use
    !    the true valence charge used for pseudopotential generation

    call rad_copy(vps%chval,vps%rho_val)

    !    Rho has to be updated as rho in the VPS/PSF file is rescaled
    !    to the charge of a neutral atom.

    
    rescale_charge = vps%gen_zval/get_valence_charge(spp)

    if(basp%pseudopotential%nicore.ne.'nc ') then
       chcore => vps%chcore
       call calculate_rho_potential(rescale_charge,vps%rho_val,vps%ve_val,irel,chcore)
    else
       call calculate_rho_potential(rescale_charge,vps%rho_val,vps%ve_val,irel)
    endif
    !    vps%Rho now contains the 'true' charge used in the pseudopotential
    !    calculation,
    
    !----------------------------------------------------------


    !---------------------Rho and Hartree pot for PAOs-----------
    !     For PAO basis functions calculations 
    ! We use the "scaled" charge density of an ion of total charge "charge"
    ! As seen above, this ion could be the one involved in ps generation,
    ! or another specified at the time of basis generation.
    ! Example: Ba: ps ionic charge: +2
    !              basis gen specified charge: +0.7
   
    !The rho for pao generation
    call rad_copy(vps%rho_val,vps%rho_val_scaled)


    rescale_charge = (get_valence_charge(spp)-charge)/vps%gen_zval

    if(basp%pseudopotential%nicore.ne.'nc ') then
      
       call calculate_rho_potential(rescale_charge,vps%rho_val_scaled,vps%ve_val_scaled,irel,chcore)
       nullify(chcore)
    else
       call calculate_rho_potential(rescale_charge,vps%rho_val_scaled,vps%ve_val_scaled,irel)
    endif

    if (charge <= 0.0_dp) then

       if (charge < 0.0_dp) then
          write(6,'(/,a)')  'ATOM: basis set generated (by rescaling the valence charge)'
          write(6,'(a,f8.4)') 'ATOM: for an anion of charge ',charge 
       endif


    else if (charge > 0.0_dp) then
       ! charge > 0
       if (abs(charge-basp%pseudopotential%zval+vps%gen_zval).gt.1.0d-3) then 
          write(6,'(/,a)') 'ATOM: basis set generated (by rescaling the valence charge)'
          write(6,'(a,f8.4)') 'ATOM: for a cation of charge ',charge 
       else
          write(6,'(/,a)') 'ATOM: basis set generated from the ionic configuration used'
          write(6,'(a)')   'ATOM: to generate the pseudopotential'
       endif

       call rad_copy(vps%ve_val, vps%ve_val_scaled)

    Endif

  end subroutine species_charge

!------------------------------------------------------------------------

  subroutine species_init(isp)
    use atm_types, only : set_symbol, set_label,set_atomic_number, &
         set_valence_charge,set_read_from_file,set_mass,get_symbol, set_lmax_orbs,&
         get_lmax_orbs, get_label
    use chemical
    
    integer , intent(in) :: isp

    integer :: iz
    type(basis_def_t), pointer :: basp
    character(len=symbol_length) :: sym

    basp => basis_parameters(isp)

    iz = basp%z
    sym = symbol(abs(iz))
    call set_atomic_number(species(isp),basp%z)
    call set_symbol(species(isp),sym)
    call set_label(species(isp),basp%label)
    call set_mass(species(isp),basp%mass)
    call set_valence_charge(species(isp),basp%pseudopotential%zval)
    call set_read_from_file(species(isp),.false.)
    call set_lmax_orbs(species(isp),basp%lmxo)

    if (basp%floating) then  
       write(6,'(3a,i4,a,a)') 'ATOM: Called for ', get_label(species(isp)), '  (Z =',iz,')',' ( Floating basis ) '

    elseif (basp%bessel) then
       write(6,'(a,i4,a)')'ATOM: Called for Z=',iz,'( Floating Bessel functions)'  
    elseif (basp%synthetic) then
       write(6,'(3a,i4,a)') 'ATOM: Called for (synthetic) ',  get_label(species(isp)),'  (Z =', iz,')'
    else
        write(6,'(3a,i4,a)') 'ATOM: Called for ', get_label(species(isp)), '  (Z =',iz,')' 
    endif

  end subroutine species_init

  !--------------------------------------------------------------------------

  subroutine species_pseudopotential_init (isp)
    !    Fills vps with info read from the psf/vps file 
    implicit none 

    integer,      intent (in)  :: isp

    !    Internal vars.

    integer :: lmax,linput,nodd,l,nrval,i,ndown,exponent

    type(pseudopotential_t), pointer :: vp
    type(basis_def_t),       pointer :: basp
    type(rad_func_t)                 :: ve
    type(rad_func_t)                 :: rad_tmp      
    type(rad_grid_t)                 :: grid
    !---

    basp => basis_parameters(isp)
    vp   => basp%pseudopotential

    linput=max(basp%lmxo,basp%lmxkb)
    lmax=min(vp%npotd-1,linput)

    if (lmax.lt.linput) then
       write(6,'(a)')  'read_vps: ERROR: You must generate a pseudopotential'
       write(6,'(a,i4)') 'read_vps: ERROR: for each L up to ',linput
       call die
    endif

    nrval = rad_get_ir_from_r(vp%vdown(0),rad_cutoff(vp%vdown(0)))
    if (rmax_radial_grid /= 0.0_dp) then
       nrval = rad_get_ir_from_r(vp%vdown(0),rmax_radial_grid)
       write(6,"(a,f10.5,i5)") &
            "Maximum radius (at nrval) set to ", &
            rmax_radial_grid, nrval
       grid = rad_get_grid(vp%vdown(0))
       call rad_set_default_length(grid,nrval)
    endif


    if (restricted_grid) then
       nodd=mod(nrval,2)
       nrval=nrval-1+nodd ! Will be less than or equal to vp%nrval
    endif

    if(nrval.gt.nrmax) then
       write(6,'(a,i4)') 'read_vps: ERROR: Nrmax must be increased to at least',nrval
       call die
    endif


    write(6,'(/,a)')'read_Read: Pseudopotential generation method:'
    write(6,'(7a)') 'read_vps: ',vp%method(1),(vp%method(i),i=3,6)

    !We are going to find the charge configuration
    !used for the pseudopotential generation using the information given in
    !the 'text' variable.

    write(6,'(a,f10.5)') 'Total valence charge: ', vp%gen_zval

    if (vp%nicore.ne.'nc  ') then
       write(6,'(/,a)') &
            'read_vps: Pseudopotential includes a core correction:'

       if(vp%nicore.eq.'pcec') then
          write(6,'(a)') 'read_vps: Pseudo-core for xc-correction'
       elseif(vp%nicore.eq.'pche') then
          write(6,'(a)')  'read_vps: Pseudo-core for hartree and xc-correction'
          write(6,'(a)') 'Siesta cannot use this pseudopotential'
          write(6,'(a)') 'Use option pe instead of ph in ATOM program'
          call die()
       elseif(vp%nicore.eq.'fcec') then
          write(6,'(a)') 'read_vps: Full-core for xc-correction'
       elseif(vp%nicore.eq.'fche') then
          write(6,'(a)') 'read_vps: Full-core for hartree and xc-correction'
          write(6,'(a)') 'Siesta cannot use this pseudopotential'
          write(6,'(a)') 'Use option pe instead of ph in ATOM program'
          call die()
       endif

    endif

    !Ionic pseudopotentials (Only 'down' used)
    do ndown=0,lmax
       l = vp%ldown(ndown) 
       if(l.ne.ndown) then
          write(6,'(a)') 'atom: Unexpected angular momentum  for pseudopotential'
          write(6,'(a)') 'atom: Pseudopotential should be ordered by increasing l'
       endif
       exponent = -1
       rad_tmp = rad_multiply_by_rl(vp%vdown(l),exponent)
       call rad_dealloc(vp%vdown(l))
       call rad_copy(rad_tmp,vp%vdown(l))
       call rad_dealloc(rad_tmp)       

    enddo

    !***  OBTAIN AN IONIC-PSEUDOPOTENTIAL IF CORE CORRECTION FOR HARTREE****
    !    POTENTIAL
    if((vp%nicore.eq.'pche').or.(vp%nicore.eq.'fche')) then

       ve = rad_vhartree(vp%chcore)


       do l=0,lmax
          rad_tmp = rad_sum(vp%vdown(l),ve)
          call rad_dealloc(vp%vdown(l))
          call rad_copy(rad_tmp,vp%vdown(l))
          call rad_dealloc(rad_tmp)
       enddo
       call rad_dealloc(ve)
    endif

    nullify(vp,basp)

  end subroutine species_pseudopotential_init

  !---------------------------------------------------------

  subroutine xc_check(xcfunc,xcauth,icorr)

    !  Checking the functional used for exchange-correlation energy.
    !Written by D. Sanchez-Portal, Aug. 1998

    character xcfunc*3, xcauth*4,icorr*2

    write(6,'(/a)') 'xc_check: Exchange-correlation functional:'
    if(((xcauth.eq.'CA').or.(xcauth.eq.'PZ')).and.&
         ((xcfunc.eq.'LDA').or.(xcfunc.eq.'LSD'))) then

       write(6,'(a)') 'xc_check: Ceperley-Alder'
       if(icorr.ne.'ca') then
          write(6,'(a)') 'xc_check: WARNING: Pseudopotential generated with'
          if(icorr.eq.'pw') write(6,'(a)')'xc_check: WARNING: Perdew-Wang 1992 functional'
          if(icorr.eq.'pb') write(6,'(a)')'xc_check: WARNING: GGA Perdew, Burke & Ernzerhof 1996 functional'
       endif

    elseif((xcauth.eq.'PW92').and. &
         ((xcfunc.eq.'LDA').or.(xcfunc.eq.'LSD'))) then

       write(6,'(a)') 'xc_check: Perdew-Wang 1992'
       if(icorr.ne.'pw') then
          write(6,'(a)')'xc_check: WARNING: Pseudopotential generated with'
          if(icorr.eq.'ca') write(6,'(a)') 'xc_check: WARNING: Ceperly-Alder functional'
          if(icorr.eq.'pb') write(6,'(a)')'xc_check:WARNING: GGA Perdew, Burke & Ernzerhof 1996 functional'
       endif
    elseif((xcauth.eq.'PBE').and.(xcfunc.eq.'GGA')) then

       write(6,'(a)') 'xc_check: GGA Perdew, Burke & Ernzerhof 1996'
       if(icorr.ne.'pb') then
          write(6,'(a)') 'xc_check: WARNING: Pseudopotential generated with'
          if(icorr.eq.'ca') write(6,'(a)') 'xc_check: WARNING: Ceperly-Alder functional'
          if(icorr.eq.'pw') write(6,'(a)') 'xc_check: WARNING: Perdew-Wang 1992 functional'
       endif

    else

       write(6,'(a)') 'xc_check: ERROR: Exchange-correlation functional not allowed'
       write(6,'(a)') 'xc_check: ERROR: xc.functional= ',xcfunc
       write(6,'(a)') 'xc_check: ERROR: xc.authors= ',xcauth
       call die

    endif

  end subroutine xc_check

end module atom

