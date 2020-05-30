module xc_flux_data
  use precision, only: dp, grid_p

  implicit none

  real(grid_p), allocatable, save :: thtr_Rho(:)
  !! Local charge Rho. Used for Hartree and XC-flux components.
  real(grid_p), allocatable, save :: thtr_Rho_deriv(:)
  !! Array to store derivative of local charge Rho
  !! computed within virtual-step MD program flow.
  !!
  !! NOTE: It is so far one-dimensional:
  !! only first spin component.
  real(grid_p), allocatable, save :: thtr_dexcdGD(:,:,:)
  !! `dexcdGD` from the call to cellXC (LibGridXC)
  real(dp), save  :: xc_flux_Jxc(3)
  !! XC-component of the heat flux.
  real(dp), save :: xc_flux_dvol
  !! Element of volume size saved from `dhscf`

contains

  subroutine init_xc_flux_data ()

    xc_flux_Jxc(:) = 0.0        ! re-init resulting Jxc-component
    allocate(thtr_Rho_deriv, source=thtr_Rho)

  end subroutine init_xc_flux_data


  subroutine reset_xc_flux_data ()

    ! In case we only use thtr_Rho for J_Hartree flux component
    ! I do the following check:
    if (allocated(thtr_Rho_deriv)) deallocate(thtr_Rho_deriv)
    ! FIXME: do not allocate/process dexcdGD on every step in `proc(dhscf)`
    if (allocated(thtr_dexcdGD)) deallocate(thtr_dexcdGD)
    deallocate(thtr_Rho)

  end subroutine reset_xc_flux_data

end module xc_flux_data
