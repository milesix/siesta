!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     


module basis
  use precision,        only: dp
  use fdf
  use hilbert_vector_collection
  use hilbert_vector_m, only:hilbert_vector_t
  use radial,           only: rad_func_t, rad_copy
  use atm_types,        only: species_info_t, get_lmax_orbs, set_orbs_deg, &
       set_lmax_orbs, species, init_orbs, get_atomic_number, set_orb
  use pseudopotential,  only:pseudopotential_t
  use atom_generation_types,      only:basis_def_t,basis_parameters,shell_t,lshell_t,energies_t
  use sys,              only:die
  use multiple_z,       only: generate_multiple_zeta
  use pol_orb,          only: generate_polarization_orbital
  use pao_util,         only: calculate_energies, normalize_orbital, &
       generate_vsoft, sym, number_of_orbs, total_charge
  use periodic_table,   only: qvlofz,cnfig
  use atom_options, only: write_ion_plot_files
  use units, only: eV
  implicit none

  private
  public basis_gen

  real(dp), parameter             :: deltmax=0.05d0

  !     Default energy-shift to define the cut off radius of orbitals
  ! In Rydbergs
  real(dp), parameter             :: eshift_default=0.02d0

  !    Default norm-percentage for the automatic definition of
  !    multiple-zeta orbitals with the 'SPLIT' option
  real(dp), parameter             :: splnorm_default=0.15d0

contains



  subroutine basis_gen(isp)
    !Generation of all the pseudo atomic orbitals

    !  There are different flavours of atomic orbitals:
    !  **PAO's orbitals (Pseudo-atomicorbitals). We can form several
    !   types of basis sets depending on how we double the basis set:
    !    1) 'NODES': Double-Z, Triple-Z,... are orbitals of PAO 
    !        type with one, two,... nodes
    !    2) 'SPLIT': Double-Z, Triple-Z,... are generated from the 
    !        original orbital, being equal to the PAO outside a 
    !        radius Rm and given by (a*r^2+b)*r^l inside Rm.
    !    3) 'NONODES': Double-Z, Triple-Z,... are orbitals of PAO
    !        type with no-nodes, but different radii or different
    !        scale-factors.
    !    4) 'SPLITGAUSS':  Double-Z, Triple-Z,... are orbitals of GTO
    !        (Gaussian Type Orbitals).
    !    5) 'USER': Read from a file, provided by the user.
    !


    !    Calculates the atomic orbitals basis set, using the option SPLIT 
    !    for the generation of the augmentation orbitals.
    !    Written by D. Sanchez-Portal, Aug. 1998
    !    Modified by DSP, July 1999



    integer, intent(in)        :: isp
    

    !***  Internal variables**

    integer :: l,iorb,norbs
    integer :: izeta, nsm, nvalence

    type(basis_def_t), pointer :: basp     !Parameters corresponding to this basis set.
    type(pseudopotential_t),pointer :: vps !Psuedopotential info.
    type(species_info_t),pointer :: spp    !Pointer to a species.
    type(shell_t),     pointer :: shell    !Pointer to a shell.
    type(lshell_t),    pointer :: lshell   !Pointer to a l-shell.
    
    type(hilbert_vector_t) :: orb_vector
    type(rad_func_t)   :: rad_tmp
    
    integer            :: atomic_number
    logical            :: filterOrbitals
    real(dp)           :: filterFactor
    real(dp),dimension(0:3) :: qatm
    integer, dimension(0:3) :: config
    logical                 :: is_pol

    basp => basis_parameters(isp)
    spp  => species(isp)
    vps  => basp%pseudopotential

    !**   Filter the orbitals
    filterOrbitals = fdf_boolean("PAO.Filter",.false.)
    filterFactor = fdf_double("PAO.FilterFactor", 0.7_dp)

    qatm(0:3)=0.0d0

    if(basp%floating) then
       atomic_number = -1*get_atomic_number(spp)
       qatm = 0.0_dp
    elseif(basp%synthetic) then
       qatm = basp%ground_state%occupation
    else
       atomic_number = get_atomic_number(spp)
       call qvlofz(atomic_number,qatm)
    endif
    
    call cnfig(atomic_number,config)

    norbs = number_of_orbs(isp)
  
    call init_orbs(species(isp),norbs)
    
    print ('(/,a)'), "BASISgen begin"

    !LOOP over angular momenta    
    iorb = 1

    do l=0,get_lmax_orbs(spp) ! species(isp)%lmax_basis
       lshell => basp%lshell(l)
      
       !loop over all the semicorestates.
       do nsm=1,lshell%nn
          
          shell => lshell%shell(nsm)
          shell%i_sm = nsm
          shell%population = 0.0_dp

          if(shell%nzeta.le.0) exit
          
          if(.not. shell%polarizes) &
               write(6,'(/A,I2)')'SPLIT: Orbitals with angular momentum L=',l

          !IF THE COMPRESSION FACTOR IS NEGATIVE OR ZERO THE ORBITALS ARE
          !LEFT UNTOUCHED
          if(shell%lambda(1).le.0.0d0) shell%lambda(1)=1.0d0
          call rad_copy(vps%vdown(l), shell%pseudo)

          shell%z_valence  = basp%ground_state%z_valence
           !Set pseudo for schroed. eq.
          call pseudopotential_setup(shell,vps)

          do izeta=1, shell%nzeta
             is_pol = .false.
             !COMPRESSION FACTOR IS ONLY ACTIVE FOR THE INITIAL PAO WHEN USING****
             !SPLIT OPTION FOR THE GENERATION OF THE BASIS SET****
             shell%lambda(izeta)=shell%lambda(1)            
             nvalence = config(shell%l)

             if(izeta.eq.1) then                
                                                             
                !Polarization
                if (shell%polarizes)then 
                   call generate_polarization_orbital(isp,shell) 
                   is_pol = .true.
                else

                   !Population analysis
                   ! floating atom?
                   if (get_atomic_number(spp) > 0) then !Not a floating atom
                      if(shell%n == nvalence)then
                         shell%population(izeta) = qatm(l)/dble(2*l+1)
                      elseif(shell%n < nvalence) then
                         shell%population(izeta) = 2.0_dp
                      endif
                   endif

                   !Automatic rc?
                   if (shell%rc(1) .eq. 0.0_dp) call auto_rc(shell,vps,l,nsm)
                   !Generate the orbital                   
                   call generate_first_zeta(shell,basp%label)                                     
                  
                endif

             else                    
                shell%population(izeta) = 0.0_dp

                if(basp%basis_type.eq.'split')then
                   shell%multiple_z_kind = 'split'                   
                elseif(basp%basis_type.eq. 'splitgauss')then
                   shell%multiple_z_kind = 'splitgauss'
                elseif(basp%basis_type.eq. 'nonodes')then
                   shell%multiple_z_kind = 'nonodes'
                elseif(basp%basis_type.eq. 'nodes')then
                   shell%multiple_z_kind = 'nodes'
                else
                   call die('basis gen: unknown multiple-z scheme')
                endif
                call generate_multiple_zeta(lshell,izeta,nsm) 

             endif

             call print_shell_info(shell,izeta)

             !filter
             if (filterorbitals) then                
                call rad_copy(shell%orb(izeta),rad_tmp)
                call rad_dealloc(shell%orb(izeta))
                shell%orb(izeta) = rad_filter(rad_tmp,l,filterFactor,2) 
                shell%rc(izeta) = rad_cutoff(shell%orb(izeta))
                call rad_dealloc(rad_tmp)
                call rad_copy(shell%orb(izeta),rad_tmp)
                call rad_dealloc(shell%rphi(izeta))
                shell%rphi(izeta) = rad_multiply_by_rl(rad_tmp,l+1)
                call rad_dealloc(rad_tmp)
                call calculate_energies(shell, izeta)
                print *, "----Filtered-orb--energies"
                call print_shell_info(shell, izeta)
             endif            

             call init_vector(orb_vector,shell%orb(izeta),shell%n,l,izeta,&
                     shell%population(izeta),shell%energies(izeta)%total,is_pol)
             !Store the orb and realese the memory
             call set_orb(spp,orb_vector,iorb)
             call destroy_vector(orb_vector)

             iorb = iorb + 1
          enddo
          
       enddo

    enddo
    
    call set_orbs_deg(spp)

    if (.not. basp%floating) then
       call species_vxc(isp)
    endif
      
    print ('(/,a)'), "BASISgen end"

  end subroutine basis_gen

  !-------------------------------------------------------------------------

  subroutine pseudopotential_setup(shell,vps)
    type(shell_t), intent(inout)         :: shell
    type(pseudopotential_t),intent(in)   :: vps    !Pseudo information

    call rad_copy(vps%ve_val,shell%ve)
    call rad_copy(vps%ve_val_scaled,shell%ve_PAO)

  end subroutine pseudopotential_setup

  !-------------------------------------------------------------------------
 
  subroutine print_shell_info(shell,izeta,message)
    type(shell_t), intent(in)            :: shell
    integer, intent(in)                  :: izeta
    character(len=*),optional,intent(in) :: message

    if(present(message)) print ('(/,a)'), message
    if(izeta.eq.1) then  

       write(6,'(/,(3x,a,i2),3(/,a25,f12.6))') 'izeta =',izeta, &
            'lambda =',shell%lambda(izeta), 'rc =',shell%rc(izeta), 'Total energy =',shell%energies(iZeta)%total  

    elseif(izeta.gt.1) then 

       write(6,'(/,(3x,a,i2),3(/,a25,f12.6))') 'izeta =',izeta,&
            'rmatch =',shell%rc(izeta),'splitnorm =',shell%spln(izeta),'Total energy =',shell%energies(iZeta)%total

       if (shell%spln(izeta) < 0.05_dp) then
          write(6,"(a)")"WARNING: effective split_norm is quite small." 
          write(6,"(a)")"         Orbitals will be very similar."
       endif

    endif

    write(6,'(a25,f12.6)') 'kinetic =',shell%energies(izeta)%ekin
    write(6,'(a25,f12.6)') 'potential(screened) =',shell%energies(izeta)%screened
    write(6,'(a25,f12.6)') 'potential(ionic) =',shell%energies(izeta)%ionic

  end subroutine print_shell_info

  !-------------------------------------------------------

  subroutine auto_rc(shell,vps,l,nsm)
    type(shell_t) :: shell
    Type(pseudopotential_t) :: vps
    integer, intent(in) :: nsm !which semicore shell we're generating
    integer, intent(in) :: l

    integer  :: nnodes, nprin
    real(dp) :: eigen, eshift, el,rc
    !real(dp), pointer, dimension(:) :: rphi

    type(rad_func_t) :: rad_orb
    !Automatic determination of the cut off radius for the PAOs
    !***  READING THE ENERGY-SHIFT TO DEFINE THE CUT-OFF RADIUS OF ORBITALS***
    eshift=fdf_physical('PAO.EnergyShift',eshift_default,'Ry')

    nnodes=nsm
    nprin=l+nsm


    rc = rad_cutoff(vps%ve_val)
    rad_orb = rad_schro(shell%pseudo,vps%ve_val,rc,l,nnodes,nprin,shell%z_valence,eigen)

    !Rc given by eshift   
    if(eigen.gt.0.0d0) then 
       write(6,'(/2A,I2,A)') 'SPLIT: ERROR ', &
            'Orbital with angular momentum L=',l,' not bound in the atom'
       write(6,'(2A)') 'SPLIT: ERROR a cut off radius must be explicitely given' 
       call die
    endif

    if(abs(eshift).gt.1.0d-5) then
       el=eigen+eshift
       shell%rc(1) = rad_rc_vs_e(rad_orb,shell%pseudo,vps%ve_val,l,el,nnodes)
    else
       shell%rc(1) = rad_default_length(rad_orb) - 0.2_dp
       !rad_get_r_from_ir(shell%pseudo, rad_default_length(shell%pseudo)-2)
    endif
    write(6,'(/,A,/,A,f10.6,A)') &
         'SPLIT: PAO cut-off radius determinated from an', &
         'SPLIT: energy shift=',eshift,' Ry'
    call rad_dealloc(rad_orb)
  end subroutine auto_rc

  !--------------------------------------------------------------------
  subroutine generate_first_zeta(shell,label) 

    type(shell_t),intent(inout)    :: shell
    character(len=*), intent(in)   :: label

    !Internal vars
    integer                        :: nprin,l
    real(dp)                       :: eorb !, dnrm
    real(dp)                       :: rc
    logical                        :: new_split_code, fix_split_table, split_tail_norm
    type(rad_func_t)               :: vtot

    split_tail_norm=fdf_boolean('PAO.SplitTailNorm',.false.)
    fix_split_table=fdf_boolean('PAO.FixSplitTable',.false.)
    if (split_tail_norm .or. fix_split_table) then
       new_split_code = .true.
    else
       new_split_code=fdf_boolean('PAO.NewSplitCode',.false.)
    endif

    l=shell%l
    write(6,'(/A,I1,A)') 'SPLIT: Basis orbitals for state ',shell%n, sym(l)

    rc=shell%rc(1)/shell%lambda(1)    
    if (rc < 0.0_dp)  call die("rc < 0 for first-zeta orbital")

    nprin=l+shell%i_sm

    !SOFT-CONFINEMENT            
    !Calculate the soft confinement potential for this shell                
    vtot = generate_vsoft(shell,label,write_ion_plot_files)
  
    shell%rphi(1) = rad_schro(shell%pseudo,vtot,rc,l,shell%i_sm,nprin,shell%z_valence,eorb)
    
    !The orb is multiplied by r and the spherical harmonics routines
    !expect the orb divided by r**l, so we divide it by r**(l+1)
    shell%orb(1) = rad_divide_by_r_l_1(shell%rphi(1),l,shell%lambda(1))

    if (split_tail_norm) then
       write(*,"(a)") "Split based on tail norm"
       ! Follows the criterion of the JPC paper, but
       ! with more contrast 
       !(use the actual math norm (sqrt(int(f^2)))
       shell%split_table = rad_split_scan_tail(shell%rphi(1))
    else
       !  Do a full scan using the old method
       !  (scan norm of tail+parabola)
       shell%split_table = rad_split_scan_tail_parabola(shell%rphi(1),l,fix_split_table)

    endif
    
    call normalize_orbital(shell,1)
    call calculate_energies(shell,1)
    shell%rc(1) = rad_cutoff(shell%orb(1))
    call rad_dealloc(vtot)

  end subroutine generate_first_zeta

  !----------------------------------------------------------------------------------

  subroutine species_vxc(is)
    !Calculates exchange correlation energy of a given atomic species
    !The result can be compared with SIESTA's result, as long as there is only
    !one  atom in the system.
    !E. Anglada 2008

    integer, intent(in) :: is !Which specie we are interested in.

    type(rad_func_t) :: rho, vxc, rho4pi
    real(dp)         :: charge, rc, exc

    rho = total_charge(is,charge,rc)
    rho4pi = rad_divide_by_4pir2(rho,.false.)
    vxc = rad_vxc(rho4pi,0,exc)

    write(6,'(/,a)') " Check species Exc"
    write(6,'(a,f14.4)') "   Species valence charge =", charge
    write(6,'(a,f14.4)') "   Species Exc (eV) =", exc/eV

    call rad_dealloc(vxc)
    call rad_dealloc(rho)
    call rad_dealloc(rho4pi)
  end subroutine species_vxc

end module basis



