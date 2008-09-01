!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

module corechg

  use precision
  use pseudopotential, only : pseudopotential_t
  use atm_types, only: species, set_has_core_charge, set_core_charge
  use atom_generation_types , only :basis_parameters
  use radial
  use fdf
  implicit none

contains

  subroutine corechargeSetup(isp)
    !    Fill species(isp)%core with the pseudo-core information contained in core.
    !    This charge density will be used in the calculation of the nonlinear core corrections.
    !    D. Sanchez-Portal, Aug. 1998
    !    A. Garcia             
    !    E. Anglada              2007-08

    implicit none

    integer, intent(in)             :: isp
   
    !    Internal variables     
    real(dp)                     :: rcore
    logical                      :: filterPCC
    type(rad_func_t), pointer    :: core_charge
    type(rad_func_t)             :: core_charge_tmp,core_charge_filtered
  
    if (basis_parameters(isp)%pseudopotential%nicore == 'nc ') then

       call set_has_core_charge(species(isp), .false.)
       
    else

       core_charge => basis_parameters(isp)%pseudopotential%chcore
       call set_has_core_charge(species(isp), .true.)
       core_charge_tmp = rad_divide_by_4pir2(core_charge,.true.)
       rcore = rad_cutoff(core_charge_tmp)
       write(6,'(a,f10.6)') 'corechg: Radius of charge density used in the'
       write(6,'(a,f10.6)') '         nonlinear core corrections =',Rcore
       
       !Filter
       filterPCC = fdf_boolean("PCC.Filter",.false.)
         
       if (filterPCC) then
          core_charge_filtered = rad_filter(core_charge_tmp,0,1.0_dp,0)   
          call set_core_charge(species(isp),core_charge_filtered)
          call rad_dealloc(core_charge_filtered)         
       else
          !Store the core charge
          call set_core_charge(species(isp),core_charge_tmp)          
       endif
       call rad_dealloc(core_charge_tmp)
       nullify(core_charge)

    endif

  end subroutine corechargeSetup


end module corechg



