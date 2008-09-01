!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

! Generation of polarization orbitals
module pol_orb
  use fdf
  use precision
  use atom_generation_types, only:basis_def_t,basis_parameters,shell_t,lshell_t,energies_t
  use pseudopotential, only:pseudopotential_t
  use atm_types, only:species_info_t,species
  use pao_util,  only:calculate_energies, normalize_orbital, sym, generate_vsoft
  use atom_options, only: write_ion_plot_files
  use radial
  use radial_logGrid
  !use schro  
  implicit none

  contains


  subroutine generate_polarization_orbital(isp,shell)

    !    Calculates the polarization  orbitals for the basis set augmentation.
    !    Written by D. Sanchez-Portal, Aug. 1998.
    !    Modify by DSP, July 1999
    !    Filter added by E. Anglada 2006

    integer, intent(in)              ::  isp  !Specie index
    type(shell_t),intent(inout)      :: shell !The shell being filled with orbitals.

    !Internal vars.
    integer                          :: lpol !l of shell being polarized.
    real(dp)                         :: rc,epao
    logical                          :: new_split_code, fix_split_table, split_tail_norm

    type(species_info_t),    pointer :: spp
    type(basis_def_t),       pointer :: basp
    type(shell_t),           pointer :: shell_pol !The shell being polarized
    type(rad_func_t)                 :: vtot

    logical                          :: old_pol = .false.

    split_tail_norm=fdf_boolean('PAO.SplitTailNorm',.false.)
    fix_split_table=fdf_boolean('PAO.FixSplitTable',.false.)
    if (split_tail_norm .or. fix_split_table) then
       new_split_code = .true.
    else
       new_split_code=fdf_boolean('PAO.NewSplitCode',.false.)
    endif

    basp => basis_parameters(isp)
    spp  => species(isp)
    
    lpol = shell%l_shell_polarized
    shell_pol => basp%lshell(lpol)%shell(1)

    rc=shell_pol%rc(1)/shell_pol%lambda(1)

    if (rc < 0.0_dp)  call die("rc < 0 for first-zeta orbital")

    write(6,'(/A,I1,A)') 'POLgen: Polarization orbital for state ',shell_pol%n , &
         sym(shell_pol%l)

    !SOFT-CONFINEMENT            
    !Calculate the soft confinement potential for this shell                
    vtot = generate_vsoft(shell_pol,basp%label,write_ion_plot_files)

    !Generate the polarization function perturbatively from the original PAO**
    ePAO = shell_pol%energies(1)%total

    old_pol = fdf_boolean('PAO.OLD_POL',.false.)

    if ( old_pol ) then

       shell%rphi(1) = rad_polarization(shell_pol%rphi(1),shell_pol%pseudo,vtot,rc,lpol,ePao)

    else
       shell%rphi(1) = rad_polarization(shell_pol%rphi(1),shell%pseudo,vtot,rc,lpol,ePao)

    endif
    shell%orb(1) = rad_divide_by_r_l_1(shell%rphi(1),lpol+1,shell%lambda(1))
    shell%rc(1) = rad_cutoff(shell%orb(1))
    
    !allocate(shell%split_table(1:nrc))
    if (split_tail_norm) then
       write(*,"(a)") "Split based on tail norm"
       ! Follows the criterion of the JPC paper, but
       ! with more contrast 
       !(use the actual math norm (sqrt(int(f^2)))
       shell%split_table = rad_split_scan_tail(shell%rphi(1)) !sqrt(max(1.0_dp-rnrm(1:nrc),0.0_dp))

    else
       !  Do a full scan of the old method
       !  (norm of tail+parabola)
       
       shell%split_table = rad_split_scan_tail_parabola(shell%rphi(1),shell%l,fix_split_table)
    endif

    !Normalization of basis functions***
    call normalize_orbital(shell,1)

    !Calculation of the mean value of kinetic and potential energy**
    call calculate_energies(shell,1)

    call rad_dealloc(vtot)
  end subroutine generate_polarization_orbital

  !-----------------------------------------------------
end module pol_orb
