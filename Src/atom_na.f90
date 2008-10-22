!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

!Generation of neutral atom potential, used to screen the hartree potential.
module atom_na
  !
  !     Generates the neutral-atom potential Vna
  !     using only the orbitals and reduced_vlocal
  !     Needs also the orbital populations
  !
  !     Needs to create an *artificial* log grid
  !
  !c     The routine gen_Vna assumes that the data structures in atm_types
  !c     (the "species" type) have been initialized and contain
  !c     all the relevant information. Thus, this routine should be
  !c     called after reading the information from file, or after
  !c     calling "atom".

  use precision
  use fdf
  use atom_types, only: get_atomic_number,get_lmax_orbs,&
       get_valence_charge, set_neutral_atom_potential,set_self_energy, &
       set_no_neutral_atom_potential, is_floating, orbs_kc_max
  use atom_generation_types, only: basis_parameters,lshell_t,shell_t,basis_def_t
  use atom_pao_util, only: total_charge
  use schro, only: vhrtre
  use radial, only : rad_func_t, rad_cutoff, rad_filter, rad_vhartree, rad_multiply_each_value, &
       rad_sum, rad_dealloc, rad_self_energy, eps
  use pseudopotential_new, only : pseudopotential_new_t, get_vlocal

  implicit none

  private
  public gen_vna

contains

  subroutine gen_Vna(is) 

    !Generates the neutral-atom pseudopotential.

    integer, intent(in)    ::  is
    
    logical                         :: filterVna
    real(dp)                        :: chval, rVna, rcocc 
    real(dp)                        :: self_energy, reference_chval,ch_temp
    type(pseudopotential_new_t),pointer :: vps
    type(rad_func_t)                :: rho, vna,vna_tmp,filt_vna, vlocal
    type(basis_def_t),  pointer     :: basp
 
    !-----------------------------------------

    !Only compute vna if it's a real species (not floating, bessel, etc)
    if (is_floating(is)) return

    Write(6,*) 'na: Computing Vna for species ' , is

    basp=>basis_parameters(is)
    vps => basis_parameters(is)%pseudopotential

    if (get_atomic_number(is) .le.0) then
       Write(6,*) "na: No Vna for species: ", is
       call set_self_energy(is,0.0_dp)
       call set_no_neutral_atom_potential(is)
       !call rad_zero(Vna)
       !call set_neutral_atom_potential(spp,vna)
       RETURN
    endif

    !     Compute total charge, using the occupied basis orbitals.
    rho =  total_charge(is,chval,rcocc)

    reference_chval = dble(get_valence_charge(is))
    write(6,'(a,2f10.5)') '     Vna: chval, zval: ', chval, reference_chval
    
    !     Make sure that it adds up to the total valence charge
    !
    !eps=1.0d-4
    if(abs(chval-reference_chval ).gt.eps) then
       ch_temp = reference_chval/chval
       rho = rad_multiply_each_value(rho,ch_temp)      
    endif

    !C CALCULATION OF THE HARTREE POTENTIAL 
    !C DUE TO THE NEW VALENCE CHARGE 
    vna_tmp = rad_vhartree(rho)

    !Calculation of real
    vlocal = get_vlocal(vps)
    vna=rad_sum(vlocal,vna_tmp,rcocc)

    rVna = rad_cutoff(vna)

    write(6,'(/,a,f10.6)') &
         'Vna:  Cut-off radius for the neutral-atom potential: ', rVna

    if(rVna.gt.(rcocc+0.5d0)) then
       write(6,"(2a,f12.5)")'Vna: WARNING: ', &
            'Cut-off radius for the neutral-atom potential, rVna =', & 
            rVna
       write(6,"(2a,f12.5)")'Vna: WARNING: ', &
            'Cut-off radius for charge density =', rcocc
       write(6,"(2a)")'Vna: WARNING: ', &
            'Check ATOM: Look for the sentence:'
       write(6,"(2a)")'Vna: WARNING: ', &
            'LOCAL NEUTRAL-ATOM PSEUDOPOTENTIAL'
       write(6,"(2a)")'Vna: WARNING: ', &
            'Increasing the tolerance parameter EPS'
       write(6,"(2a)")'Vna: WARNING: ', &
            'might be a good idea'
    endif

   
    !Filter the high-k components of Vna
    filterVna = fdf_boolean("Vna.Filter",.false.)

    if (filterVna)then
       filt_vna = rad_filter(vna,0,1.0_dp,0,orbs_kc_max)
       call set_neutral_atom_potential(is,filt_vna)
       call rad_dealloc(filt_vna)
    else
       !Store the neutral atom potential
       call set_neutral_atom_potential(is,vna)
    endif

    !Self energy
    self_energy = rad_self_energy(vlocal)
    call set_self_energy(is,self_energy)

    !Release memory
    call rad_dealloc(rho)
    call rad_dealloc(vna)
    call rad_dealloc(vna_tmp)
    call rad_dealloc(vlocal)
    
    Write(6,*) 'na: Finished computing Vna for species ' , is
    
  end subroutine gen_Vna
 
  !=================================================================

end module atom_na


