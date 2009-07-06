!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

!     This module is the main interface to the atomic part of siesta.
!     The only public subroutine is:
!     subroutine generate_all_atomic_info()
!     No arguments.

module atom
  !Common modules
  use precision, only: dp
  use sys, only: die
  use fdf
  use units, only: eV
  
  !Generic module for radial (1D) functions
  use radial, only : rad_func_t, restricted_grid, rmax_radial_grid, rad_cutoff, &
       rad_get_ir_from_r, rad_get_grid, rad_multiply_by_rl, rad_dealloc, rad_sum,&
       rad_set_default_length, rad_vhartree, rad_multiply_each_value, &
       rad_divide_by_4pir2, rad_vxc, rad_copy
  
  !Main module where all the atomic info is stored.
  use atom_types, only : get_atomic_number, get_symbol, &
       symbol_length, set_no_reduced_vlocal, set_no_kb, &
       set_no_neutral_atom_potential, has_kbs, set_floating, &
       set_has_core_charge, filter_orbs, get_number_of_species, rad_grid_t
  
  !Auxiliary module where all the intermediate info is stored.
  use atom_generation_types, only : basis_def_t, basis_parameters,lshell_t, shell_t, &
       write_basis_specs
  
  !All the auxiliary modules
  use atomparams, only : 
  use periodic_table, only : symbol
  use pseudopotential_new, only : pseudopotential_new_t, get_pseudo_valence_charge,&
       get_pseudo_gen_valence_charge, pseudo_has_core_charge, pseudo_copy_scaled_ve,&
       get_pseudo_method, get_pseudo_npotd, get_pseudo_nicore, pseudo_dealloc_vdown,&
       pseudo_set_vdown, get_pseudo_l, get_pseudo_down, set_rho_val,&
       pseudo_dealloc_vdown, pseudo_set_vdown, get_core_charge, get_irel, get_icorr,&
       get_rho_val, pseudo_copy_valence_charge, set_ve_val_scaled, set_ve_val,&
       get_rho_val_scaled, pseudo_copy_scaled_charge, set_rho_val_scaled
  
  !All the expert subroutines which generate the corresponding functions.
  use atom_corecharge, only: core_charge_Setup
  
  use atom_bessel, only : bessel
  use parallel, only : ionode
  use atom_na, only : gen_vna
  use atom_options, only : get_atom_options, write_ion_plot_files
  use atomparams, only : nzetmx, nrmax
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

    use atom_kb,    only: kbgen
    use atom_vloc,  only: gen_vlocal 
    use atom_basis_gen, only: basis_gen
    use ldauproj, only: ldauproj_gen

    integer , intent(in) :: isp
    logical :: floating = .false.

    write(6,'(/,a)') 'ATOM: Species begin__________________________ '

    if (basis_parameters(isp)%bessel .or. basis_parameters(isp)%floating) floating = .true.
    call species_init(isp,floating)
   

    if (basis_parameters(isp)%bessel) then
       call set_has_core_charge(isp, .false.)
       call bessel(isp)
    elseif (basis_parameters(isp)%floating) then
       call species_pseudopotential_init(isp)
       call species_charge(isp)
       call basis_gen(isp)
    else
       call species_pseudopotential_init(isp)
       call species_charge(isp)    
       call gen_vlocal(isp)
       call kbgen(isp)
       call basis_gen(isp)
       call ldauproj_gen(isp)
    endif
    write(6,'(/,a)') 'ATOM: Species end_____________________________ '

  end subroutine atom_main

  !----------------------------------------------------

  subroutine calculate_rho_potential(charge,rho,potential,irel,rho_core)
    !Scale rho with charge
    !Calculate Vhartree and Vxc
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

    !-------Scale rho with the new charge
    rho_scaled = rad_multiply_each_value(rho,charge)
    call rad_dealloc(rho)
     
     !------Hartree-potential-----------------
    pot_hartree = rad_vhartree(rho_scaled)
    !call rad_dump_file(pot_hartree,"pot_hartree.dat")

    ! NOTE that atomxc expects true rho(r), not 4pir^2*rho(r)
    rho4pir2=rad_divide_by_4pir2(rho_scaled,.false.)

    !------XC-potential----------------------
    !If exists add the core charge due to nonlinear core correcctions
    if (present(rho_core)) then
       rho_core4pir2 = rad_divide_by_4pir2(rho_core,.false.)
       !call rad_dump_file(rho_core4pir2,"rho_core4pir2.dat")
       rho_tot = rad_sum(rho4pir2,rho_core4pir2,rmax)
       !call rad_dump_file(rho_tot,"rho_tot.dat")
       call rad_dealloc(rho_core4pir2)
    else
       call rad_copy(rho4pir2,rho_tot)
    endif

    !------Compute xc energy and potential and add to ve
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

    use atom_basis_specs, only: read_basis_specs
    use atom_basis_io, only: read_basis_ascii, dump_basis_ascii,dump_basis_xml, dump_basis_netcdf, read_basis_netcdf
    use atom_elec_correction, only: elec_corr_setup

    implicit none

    ! Internal variables ...................................................
    integer  :: nspecies, is
    real(dp) :: filter_factor

    logical :: user_basis, user_basis_netcdf, read_from_file, filter_Orbitals

    write(6,'(/2a)')'initatom: Reading input for the', &
         'pseudopotentials and atomic orbitals'

    user_basis = fdf_boolean('user-basis',.false.)
    user_basis_netcdf = fdf_boolean('user-basis-netcdf',.false.)

    read_from_file = user_basis .or. user_basis_netcdf

    call get_atom_options()

    if (user_basis_netcdf) then

       write(6,'(a)') 'Reading PAOs and KBs from NetCDF files...'
       call read_basis_netcdf()

    else if (user_basis) then

       write(6,'(a)') 'Reading PAOs and KBs from ascii files...'
       call read_basis_ascii()

    else

       !     Read basis specification and generate PAOs and 
       !     KB and LDAU projectors
       !     Note : In this section, the number of species is nspecies
       call read_basis_specs()

       nspecies = get_number_of_species()

       do is = 1,nspecies
          call write_basis_specs(6,is)
          call atom_main(is)
       enddo

       filter_Orbitals = fdf_boolean("PAO.Filter",.false.)
       if (filter_orbitals) then
          filter_factor = fdf_double("PAO.FilterFactor",0.7_dp)
          call filter_orbs(filter_factor)
       endif

       !Nonlinear core corrections, neutral atom potential 
       do is = 1, nspecies
          ! If exists, fill species core density (for nonlinear core corrections)
          call core_charge_setup(is)
          call gen_vna(is)
       enddo

       !     Compute the electrostatic correction tables
       call elec_corr_setup()

       call prinput(nspecies)
    endif
    
    
       ! Dump all the atomic information
       call dump_basis_ascii()
       call dump_basis_netcdf()

       if (write_ion_plot_files) then         
          call dump_basis_xml()
       endif

    !Subroutine which dumps all the atomic info. Do not remove, please.
    !call check_atmfuncs()
  end subroutine generate_all_atomic_info

  !----------------------------------------------------

  subroutine prinput(ntotsp)

!**
! Prints the values of the parameter which have been actually 
! used in the generation of the basis (cut-off radius, contraction
! factors, and BasisType option to augment the basis set).
! The information is written in the same format as required for the 
! input file.
! Written by D. Sanchez-Portal, Oct. 1998.
!**

    implicit none 
    integer ntotsp

        
    !***Internal variables
    integer is, nshells, l, izeta
          
    character(len=11) rcchar(nzetmx), lambdachar(nzetmx)

    type (basis_def_t), pointer :: basp
    type (shell_t)    , pointer :: shell
    integer                     :: n
    
    write(6,'(/a,58("-"))') 'prinput: Basis input '


    write(6,'(/2a)')'PAO.BasisType split'

    write(6,'(/a)') '%block ChemicalSpeciesLabel' 
    do is=1,ntotsp
       basp => basis_parameters(is)
       write(6,'(2(1x,i4),1x,2a)') is,basp%z,basp%label, &
            '    # Species index, atomic number, species label' 
       nullify(basp)
    enddo
    write(6,'(a)') '%endblock ChemicalSpeciesLabel' 

    write(6,'(/a)') '%block PAO.Basis                 # Define Basis set'
    do is=1, ntotsp  
       basp => basis_parameters(is)
       nshells=0 

       do l=0,basp%lmxo !Find number of shells (including semicore)
          do n=1,basp%lshell(l)%nn
             shell => basp%lshell(l)%shell(n)
             if(shell%nzeta > 0) nshells=nshells+1
             nullify(shell)
          enddo
       enddo

       if(basp%basis_type .eq. 'split') then  

          if(abs(basp%ionic_charge).lt.1.0d-4) then 
             write(6,'(a10,1x,i2,20x,a)') basp%label, nshells, &
                  '# Species label, number of l-shells'
          else
             write(6,'(a10,1x,i2,1x,f7.3,12x,2a)') basp%label, nshells, &
                  basp%ionic_charge, '# Label, l-shells,',' ionic net charge'
          endif
         
       else 

          if(abs(basp%ionic_charge).lt.1.0d-4) then
             write(6,'(a10,1x,i2,1x,a,10x,2a)') &
                  basp%label, nshells,basp%basis_type, &
                  '# Species label, l-shells,',' basis type ' 
          else
             write(6,'(a10,1x,i2,1x,a,1x,f7.3,1x,2a)') &
                  basp%label, nshells, basp%basis_type, &
                  basp%ionic_charge, '# Label, l-shells, type,', &
                  ' ionic net charge'
          endif
           
       endif
 
       do l=0,basp%lmxo
          do n=1,basp%lshell(l)%nn
             shell => basp%lshell(l)%shell(n)

             if(shell%auto_polarizes .or. shell%nzeta == 0) then
                exit
             else
                if (shell%auto_polarized) then
                   write(6,'(1x,a,i1,2(1x,i3),a,i3,19x,2a)') &
                        'n=',shell%n,l,shell%nzeta, ' P ',shell%nzeta_pol,&
                        '# n, l, Nzeta, ','Polarization, NzetaPol'
                else
                   write(6,'(1x,a,i1,2(1x,i3),25x,a)') &
                        'n=',shell%n,l, shell%nzeta, '# n, l, Nzeta '
                endif
                do izeta=1, shell%nzeta
                   write(rcchar(izeta),'(1x,f7.3)') shell%rc(izeta)
                   write(lambdachar(izeta),'(1x,f7.3)') shell%lambda(izeta)
                enddo
                write(6,'(20a)') (rcchar(izeta), izeta=1,shell%nzeta)
!                    ,'        # rc(izeta=1,Nzeta)(Bohr)'
                write(6,'(20a)') (lambdachar(izeta), izeta=1,shell%nzeta)
!                    ,'        # scaleFactor(izeta=1,Nzeta)'
             endif
          enddo
       enddo
    enddo
    write(6,'(a)') '%endblock PAO.Basis' 

    write(6,'(/a,70("-")/)') 'prinput: '

  end subroutine prinput



  !----------------------------------------------------
  subroutine species_charge(isp)
    use atom_types, only : get_valence_charge
    integer, intent(in) :: isp

    !Internal vars
    real(dp)                         :: charge 
    type(pseudopotential_new_t),pointer  :: vps => Null()
    type(basis_def_t), pointer       :: basp => Null()
    integer                          :: irel
    type(rad_func_t)                 :: chcore, rho_val, ve_val, rho_val_scaled, ve_val_scaled
    real(dp)                         :: rescale_charge, zval, gen_zval
    character(len=3)                 :: xcfunc, rel_label    !XC used, relativistic?
    character(len=4)                 :: xcauth               !Authors of xc
    character(len=2)                 :: pseudo_XC_label      !XC used in the pseudo
    
    basp => basis_parameters(isp)
    vps => basp%pseudopotential

    zval = get_pseudo_valence_charge(vps)
    gen_zval = get_pseudo_gen_valence_charge(vps)

    if (abs(get_atomic_number(isp)-gen_zval).gt.1.0d-3) then 
       write(6,'(/,a,/,a,f5.2)') &
            'ATOM: Pseudopotential generated from an ionic configuration',&
            'ATOM: with net charge', get_valence_charge(isp)-gen_zval
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
         (abs(charge-zval+gen_zval) .lt. 1.0d-3) ) then   
       ! We can't overwrite basp%ionic_charge because it may be used 
       ! if there are semicore states 
       charge = get_valence_charge(isp) - gen_zval
    endif

    !    Relativistic pseudo
    rel_label = get_irel(vps)
    pseudo_XC_label = get_icorr(vps)
    if (rel_label.eq.'rel') irel=1
    if (rel_label.ne.'rel') irel=0

    xcfunc = fdf_string('xc.functional','LDA')
    xcauth = fdf_string('xc.authors','PZ')

    call xc_check(xcfunc,xcauth,pseudo_XC_label)

    !-----------------------Rho, hartree potential for KBs-------------------
    !    CALCULATION OF THE VALENCE SCREENING POTENTIAL FROM THE READ CHARGE
    !    DENSITY
    !    
    !    For Kleinman-Bylander projectors calculation, use
    !    the true valence charge used for pseudopotential generation

    !call rad_copy(vps%chval,vps%rho_val)
    call pseudo_copy_valence_charge(vps)

    !    Rho has to be updated as rho in the VPS/PSF file is rescaled
    !    to the charge of a neutral atom.

    
    rescale_charge = gen_zval/get_valence_charge(isp)

    !if(basp%pseudopotential%nicore.ne.'nc ') then
    rho_val = get_rho_val(vps)

    if(pseudo_has_core_charge(vps))then
       chcore = get_core_charge(vps)
       call calculate_rho_potential(rescale_charge,rho_val,ve_val,irel,chcore)
    else
       call calculate_rho_potential(rescale_charge,rho_val,ve_val,irel)
    endif

    call set_rho_val(vps, rho_val)
    call set_ve_val(vps, ve_val)

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
    call pseudo_copy_scaled_charge(vps)

    rescale_charge = (get_valence_charge(isp)-charge)/gen_zval
    rho_val_scaled = get_rho_val_scaled(vps)

    if(pseudo_has_core_charge(vps)) then      
       call calculate_rho_potential(rescale_charge,rho_val_scaled,ve_val_scaled,irel,chcore)
       call rad_dealloc(chcore) !no longer needed.
    else
       call calculate_rho_potential(rescale_charge,rho_val_scaled,ve_val_scaled,irel)
    endif

    call set_rho_val_scaled(vps, rho_val_scaled)
    call set_ve_val_scaled(vps, ve_val_scaled)

    if (charge <= 0.0_dp) then

       if (charge < 0.0_dp) then
          write(6,'(/,a)')  'ATOM: basis set generated (by rescaling the valence charge)'
          write(6,'(a,f8.4)') 'ATOM: for an anion of charge ',charge 
       endif


    else if (charge > 0.0_dp) then
       
       if (abs(charge-zval+gen_zval).gt.1.0d-3) then 
          write(6,'(/,a)') 'ATOM: basis set generated (by rescaling the valence charge)'
          write(6,'(a,f8.4)') 'ATOM: for a cation of charge ',charge 
       else
          write(6,'(/,a)') 'ATOM: basis set generated from the ionic configuration used'
          write(6,'(a)')   'ATOM: to generate the pseudopotential'
       endif

       !call rad_copy(vps%ve_val, vps%ve_val_scaled)
       !call pseudo_copy_scaled_ve(vps)
    Endif

    call rad_dealloc(rho_val)
    call rad_dealloc(ve_val)
    call rad_dealloc(rho_val_scaled)
    call rad_dealloc(ve_val_scaled)
  end subroutine species_charge

!------------------------------------------------------------------------

  subroutine species_init(isp,floating)
    use atom_types, only : set_symbol, set_label,set_atomic_number, &
         set_valence_charge,set_read_from_file,set_mass,get_symbol, set_lmax_orbs,&
         get_lmax_orbs, get_label, set_lmax_ldau_proj

    use chemical
    
    integer , intent(in) :: isp
    logical , intent(in) :: floating

    integer :: iz
    real(dp) :: zval
    type(basis_def_t), pointer :: basp
    character(len=symbol_length) :: sym

    basp => basis_parameters(isp)

    iz = basp%z
    sym = symbol(abs(iz))
    call set_atomic_number(isp,basp%z)
    call set_symbol(isp,sym)
    call set_label(isp,basp%label)
    call set_mass(isp,basp%mass)
    zval = get_pseudo_valence_charge(basp%pseudopotential)
    call set_valence_charge(isp,zval)
    call set_read_from_file(isp,.false.)
    call set_lmax_orbs(isp,basp%lmxo)

    if (basp%floating) then  
       write(6,'(3a,i4,a,a)') 'ATOM: Called for ', get_label(isp), '  (Z =',iz,')',' ( Floating basis ) '

    elseif (basp%bessel) then
       write(6,'(a,i4,a)')'ATOM: Called for Z=',iz,'( Floating Bessel functions)'  
    elseif (basp%synthetic) then
       write(6,'(3a,i4,a)') 'ATOM: Called for (synthetic) ',  get_label(isp),'  (Z =', iz,')'
    else
        write(6,'(3a,i4,a)') 'ATOM: Called for ', get_label(isp), '  (Z =',iz,')' 
    endif

    if (floating) then
       call set_floating(isp,.true.)
       call set_no_reduced_vlocal(isp)
       call set_no_kb(isp)
       call set_no_neutral_atom_potential(isp)
    endif
  end subroutine species_init

  !--------------------------------------------------------------------------

  subroutine species_pseudopotential_init (isp)
    !    Fills vps with info read from the psf/vps file 
    implicit none 

    integer,      intent (in)  :: isp

    !    Internal vars.

    integer :: lmax,linput,nodd,l,nrval,i,ndown,exponent, npotd
    real(dp) ::gen_zval

    type(pseudopotential_new_t), pointer :: vp
    type(basis_def_t),       pointer :: basp
    type(rad_func_t)                 :: ve
    type(rad_func_t)                 :: rad_tmp,vdown, chcore      
    type(rad_grid_t)                 :: grid
    character(len=10)                :: method(6)
    character(len=4)                 :: nicore
    !---

    basp => basis_parameters(isp)
    vp   => basp%pseudopotential

    linput=max(basp%lmxo,basp%lmxkb)
    npotd = get_pseudo_npotd(vp)
    lmax=min(npotd-1,linput)

    if (lmax.lt.linput) then
       write(6,'(a)')  'read_vps: ERROR: You must generate a pseudopotential'
       write(6,'(a,i4)') 'read_vps: ERROR: for each L up to ',linput
       call die
    endif

    vdown = get_pseudo_down(vp,0)

    nrval = rad_get_ir_from_r(vdown,rad_cutoff(vdown))
    if (rmax_radial_grid /= 0.0_dp) then
       nrval = rad_get_ir_from_r(vdown,rmax_radial_grid)
       write(6,"(a,f10.5,i5)") &
            "Maximum radius (at nrval) set to ", &
            rmax_radial_grid, nrval
       grid = rad_get_grid(vdown)
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

    method = get_pseudo_method(vp)
    write(6,'(/,a)')'read_Read: Pseudopotential generation method:'
    write(6,'(7a)') 'read_vps: ',method(1),(method(i),i=3,6)

    !We are going to find the charge configuration
    !used for the pseudopotential generation using the information given in
    !the 'text' variable.

    gen_zval = get_pseudo_gen_valence_charge(vp)
    write(6,'(a,f10.5)') 'Total valence charge: ', gen_zval

    if (pseudo_has_core_charge(vp)) then
       write(6,'(/,a)') &
            'read_vps: Pseudopotential includes a core correction:'

       nicore = get_pseudo_nicore(vp)
       if(nicore.eq.'pcec') then
          write(6,'(a)') 'read_vps: Pseudo-core for xc-correction'
       elseif(nicore.eq.'pche') then
          write(6,'(a)')  'read_vps: Pseudo-core for hartree and xc-correction'
          write(6,'(a)') 'Siesta cannot use this pseudopotential'
          write(6,'(a)') 'Use option pe instead of ph in ATOM program'
          call die()
       elseif(nicore.eq.'fcec') then
          write(6,'(a)') 'read_vps: Full-core for xc-correction'
       elseif(nicore.eq.'fche') then
          write(6,'(a)') 'read_vps: Full-core for hartree and xc-correction'
          write(6,'(a)') 'Siesta cannot use this pseudopotential'
          write(6,'(a)') 'Use option pe instead of ph in ATOM program'
          call die()
       endif

    endif

    !Ionic pseudopotentials (Only 'down' used)
    do ndown=0,lmax
       l = get_pseudo_l(vp,ndown)
       !l = vp%ldown(ndown) 
       if(l.ne.ndown) then
          write(6,'(a)') 'atom: Unexpected angular momentum  for pseudopotential'
          write(6,'(a)') 'atom: Pseudopotential should be ordered by increasing l'
       endif
       exponent = -1
       vdown = get_pseudo_down(vp,l)
       rad_tmp = rad_multiply_by_rl(vdown,exponent)
       !call rad_dealloc(vp%vdown(l))
       !call rad_copy(rad_tmp,vp%vdown(l))
       call pseudo_dealloc_vdown(vp,l)
       call pseudo_set_vdown(vp,l,rad_tmp)
       call rad_dealloc(rad_tmp)       
    enddo

    !***  OBTAIN AN IONIC-PSEUDOPOTENTIAL IF CORE CORRECTION FOR HARTREE****
    !    POTENTIAL
    if((nicore.eq.'pche').or.(nicore.eq.'fche')) then

       chcore = get_core_charge(vp)
       ve = rad_vhartree(chcore)

       do l=0,lmax
          vdown = get_pseudo_down(vp,l)
          rad_tmp = rad_sum(vdown,ve)
          call pseudo_dealloc_vdown(vp,l)
          !call rad_dealloc(vp%vdown(l))
          call pseudo_set_vdown(vp,l,rad_tmp)
          !call rad_copy(rad_tmp,vp%vdown(l))
          call rad_dealloc(rad_tmp)
          call rad_dealloc(vdown)
       enddo
       call rad_dealloc(ve)
       call rad_dealloc(chcore)

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
     elseif ((xcauth.eq.'RPBE').and.(xcfunc.eq.'GGA')) then

       write(6,'(a)')  &
          'xc_check: GGA  Hammer, Hansen and Norskov, PRB 59, 7413 (1999)'

     elseif ((xcauth.eq.'REVPBE').and.(xcfunc.eq.'GGA')) then

       write(6,'(a)') 'xc_check: GGA Zhang and  Yang, PRL 80,890(1998)'

     elseif ((xcauth.eq.'LYP').and.(xcfunc.eq.'GGA')) then

       write(6,'(a)')  &
        'xc_check:  GGA Becke-Lee-Yang-Parr'

       elseif ((xcauth.eq.'WC').and.(xcfunc.eq.'GGA')) then

       write(6,'(a)') 'xc_check: GGA Wu-Cohen'

    else

       write(6,'(a)') 'xc_check: ERROR: Exchange-correlation functional not allowed'
       write(6,'(a)') 'xc_check: ERROR: xc.functional= ',xcfunc
       write(6,'(a)') 'xc_check: ERROR: xc.authors= ',xcauth
       call die

    endif

  end subroutine xc_check

end module atom

