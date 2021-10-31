module thermal_flux
  !! Core virtual_step_procs are responsible for saving/restoring
  !! the state of the system between `BASE` and `VIRTUAL` steps
  !! in VMD scheme, and for calling subroutines that implement
  !! VMD logic for corresponding substep.
  !! For the setup phase, they also aim to replicate some functionality
  !! of reciprocal lattice-related routines similar to what is realized
  !! for QE in `gvect`, `recvec` etc.
  !! Since some of these entities are common for separate components
  !! calculated in VMD-scheme, but not used in general in Siesta,
  !! and since they are bound to calculation of `charge_g`,
  !! I place them here.
  !!
  !!NOTE: wherever the comments address `VMD' (Virtual-step MD)
  !! it is the former codename of the substeps scheme for taking
  !! derivatives of functions within SIESTA DFT method.

  use precision, only: dp, grid_p
  use siesta_geom, only: xa, va, na_u
  ! use siesta_options, only: virtual_dt

  use thermal_flux_data
  use derivation_routines, only: apply_derivation
  ! use m_siesta_forces
  ! external :: die

contains

  subroutine init_thermal_flux()
    !! Reads coordinates and velocities of the snapshot
    !! (formely: `BASE` step), calls GK settings object init.
    allocate(xa_in(3, na_u))
    allocate(va_in(3, na_u))

    call read_xvs()             ! will fail if coordinates not specified

    call gk_setup%init_thermal_flux_settings()
    call gk_results%init_thermal_flux_results()
  end subroutine init_thermal_flux


  subroutine init_data_base_step()
    !! Logic to be executed at the end of the 0-`base' substep.
    use thermal_flux_jion,  only: init_ion_flux_data
    use thermal_flux_jzero, only: init_zero_flux_data

    call init_ion_flux_data()
    call init_zero_flux_data()
    ! allocate(thtr_Rho_deriv, source=thtr_Rho) <- no need; done in `dhscf'
  end subroutine init_data_base_step


  subroutine init_data_each_step()
    !! Logic to be executed at the end of each substep.
    use thermal_flux_jhart, only: store_vhart_substep

    call store_vhart_substep()
  end subroutine init_data_each_step


  subroutine read_xvs()
    integer        :: fu, rc
    integer        :: ia

    character(len=*), parameter   :: FILE_NAME = 'XVS.dat'

    open (action='read', file=FILE_NAME, iostat=rc, newunit=fu)

    if (rc /= 0) then
       call die('Thermal Transport requires "' // FILE_NAME // '" file, aborting')
    end if

    do ia=1,na_u
       read (fu, *, iostat=rc) xa_in(1:3, ia), va_in(1:3, ia)
       if (rc /= 0) call die('Error reading parameters from ' // FILE_NAME)
    end do

    close (fu)
  end subroutine read_xvs


  subroutine make_small_move()
    integer  :: step_interval

    step_interval = gk_setup%func_step_intervals(substep)

    xa(:,:) = xa_in(:,:) + step_interval * va_in(:,:) * gk_setup%virtual_dt

  end subroutine make_small_move


  subroutine compute_derivatives()

    call apply_derivation(DM_save, gk_setup%virtual_dt, gk_setup%dpoints)
    call apply_derivation(Rho_save, gk_setup%virtual_dt, gk_setup%dpoints)
    call apply_derivation(Vhart_save, gk_setup%virtual_dt, gk_setup%dpoints)

  end subroutine compute_derivatives


  subroutine compute_flux()
    use thermal_flux_jks
    use thermal_flux_jhart
    use thermal_flux_jion
    use thermal_flux_jzero

    call compute_jks()
    call compute_jhart()

    compute_jxc: do i=1,size(Rho_save, 1)
       gk_results%Jxc(1:3) = gk_results%Jxc(1:3) &
            & - Rho_save(i,2) * thtr_dexcdGD(i,1:3,1) * xc_flux_dvol
    end do compute_jxc

    call compute_jion()
    call compute_jzero()

  end subroutine compute_flux


  subroutine thermal_flux_tesults()

    call gk_results%write_thermal_flux_results()

  end subroutine thermal_flux_tesults


  subroutine reset_thermal_flux()
    use alloc,       only : de_alloc

    deallocate(xa_in)
    deallocate(va_in)
    deallocate(DM_save)

    if (allocated(thtr_dexcdGD)) deallocate(thtr_dexcdGD)
    deallocate(Rho_save)

    deallocate(Vhart_save)
    call de_alloc(charge_g, routine="reset_thermal_flux")
    deallocate(charge_g_base)

    deallocate(I_first_g, I_second_g)

  end subroutine reset_thermal_flux

end module thermal_flux
