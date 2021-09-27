module m_thermal_flux_settings

  use precision, only: dp, grid_p
  use fdf
  implicit none

  public

  type :: thermal_flux_settings_type
     logical :: init = .false.  !! Flag to show if settings object was initialized
     integer :: dpoints !! Derivatiion scheme indicator (number of points): 2, 3 or 5
     integer, allocatable :: func_step_intervals(:) !! Array of 0-centered indexes of dimension `dpoints`
     !! Used to get system coordinates for corresponding 'virtual' steps
     real(dp) :: virtual_dt   !! Virtual time step for projection of positions for time derivatives
   contains
     procedure :: init_self_deriv_scheme
     procedure :: init_thermal_flux_settings
  end type thermal_flux_settings_type

contains

  subroutine init_self_deriv_scheme(this, dpts)
    class(thermal_flux_settings_type), intent(inout) :: this
    integer, intent(in) :: dpts

    select case (dpts)
    case (2)
       write(*,*) "derivation scheme: 2-step derivative init"
       allocate(this%func_step_intervals(2))
       this%func_step_intervals(1:2) = [0, 1]
    case (3)
       write(*,*) "derivation scheme: 3-step midpoint derivative init"
       allocate(this%func_step_intervals(3))
       this%func_step_intervals(1:3) = [0, -1, 1]
    case (5)
       write(*,*) "derivation scheme: 5-step midpoint derivative init"
       allocate(this%func_step_intervals(5))
       this%func_step_intervals(1:5) = [0, -2, -1, 1, 2]
    case default
       call die('derivation scheme: ERROR points number must be 2, 3 or 5')
    end select

    this%dpoints = dpts

  end subroutine init_self_deriv_scheme


  subroutine init_thermal_flux_settings(this)
    class(thermal_flux_settings_type), intent(inout) :: this
    integer :: dpts_in

    if (.not.(this%init)) then
       write(*,*) "======= Green-Kubo ThermalFlux Calculation ======="
       dpts_in = fdf_get('ThermalFlux.NumDerivPoints', 3)
       call this%init_self_deriv_scheme(dpts_in)

       this%virtual_dt = fdf_get('ThermalFlux.Virtual.dt',0.1_dp,'fs')
       write(*,*) "substep delta_t: ", this%virtual_dt, " fs"

       this%init = .true.
       write(*,*) "done with thermal flux init"
       write(*,*) "=================================================="
    end if
  end subroutine init_thermal_flux_settings

end module m_thermal_flux_settings


module thermal_flux_data

  use precision, only: dp, grid_p
  use m_thermal_flux_settings

  type(thermal_flux_settings_type) :: gk_setup
  !! Short distinct name of the thermal
  !! transport settings' singleton example

  integer :: substep !! Current substep number
  !! (will go from 1 to `gk_setup%dpoints` in runtime)

  real(dp), allocatable :: xa_in(:,:) !! system coordinates input for the snapshot
  real(dp), allocatable :: va_in(:,:) !! system velocities input for the snapshot

  real(dp), allocatable :: DM_save(:,:,:,:) !! Full density matrix storage +1 dim:
  !! Last dimension goes over substep points
  !! (the derivation scheme parameter)

end module thermal_flux_data
