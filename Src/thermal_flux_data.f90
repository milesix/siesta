module m_thermal_flux_settings
  !! Variables and data structures for the Green-Kubo Thermal Transport
  !! calculations within SIESTA.
  !!
  !!NOTE: I keep reciprocal space G-vectors storage variables here and not
  !! in the dedicated module (`thermal_gvecs') also here so far. These data
  !! should be renamed/moved if used elsewhere outside Thermal routines.

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
     real(dp) :: meshcutoff !! MeshCutoff value copy.
     real(dp) :: alat !! LatticeConstant value copy.
     !!NOTE: at the moment it must be set in the .fdf-file for
     !! the cubic cell parameter used in Ewald scheme for Jion,
     !! even though positions and velocities are then read from `XVS.dat'-file.
     real(dp) :: eta_ewald = 0.1_dp
     !! Ewald factor for convergence.
     !! Read from .fdf input under `ThermalFlux.Jion.Eta`
     integer  :: n_max_ewald = 5
     !! Number of periodic cell images for Ewald scheme
     !! Read from .fdf input under `ThermalFlux.Jion.Nmax`
     logical :: use_sternheimer =.false.  !! Use Sternheimer equation method
     !! for `psi_hat_c`-s in computation of Jks
     logical :: qeheat_units   = .false.  !! Use QuantumEspresso Heat suite thermal flux units for comparison
     logical :: verbose_output = .false.  !! Eval and print extra debug test info
   contains
     procedure :: init_self_deriv_scheme
     procedure :: init_thermal_flux_settings
  end type thermal_flux_settings_type

  type :: thermal_flux_results_type
     real(dp) :: Jtotal(3), Jhart(3), Jxc(3)
     real(dp) :: Jks(3), Jks_A(3), Jks_B(3), Jele(3)
     real(dp) :: Jion(3), Jion_A(3), Jion_B(3), Jion_C(3), Jion_D(3), Jion_E(3)
     real(dp) :: Jzero(3), Jzloc(3), Jznl(3)
   contains
     procedure :: init_thermal_flux_results
     procedure :: write_thermal_flux_results
  end type thermal_flux_results_type

contains

  subroutine init_self_deriv_scheme(this, dpts)
    class(thermal_flux_settings_type), intent(inout) :: this
    integer, intent(in) :: dpts

    select case (dpts)
    case (2)
       write(*,"(2X,A)") "[gk:.init] Derivation scheme: 2-step derivative init"
       allocate(this%func_step_intervals(2))
       this%func_step_intervals(1:2) = [0, 1]
    case (3)
       write(*,"(2X,A)") "[gk:.init] Derivation scheme: 3-step midpoint derivative init"
       allocate(this%func_step_intervals(3))
       this%func_step_intervals(1:3) = [0, -1, 1]
    case (5)
       write(*,"(2X,A)") "[gk:.init] Derivation scheme: 5-step midpoint derivative init"
       allocate(this%func_step_intervals(5))
       this%func_step_intervals(1:5) = [0, -2, -1, 1, 2]
    case default
       call die('[gk:.init] Derivation scheme: ERROR points number must be 2, 3 or 5')
    end select

    this%dpoints = dpts

  end subroutine init_self_deriv_scheme


  subroutine init_thermal_flux_settings(this)
    class(thermal_flux_settings_type), intent(inout) :: this
    integer :: dpts_in

    if (.not.(this%init)) then
       write(*,"(2X,A)") "============= Green-Kubo ThermalFlux Calculation ============="
       dpts_in = fdf_get('ThermalFlux.NumDerivPoints', 3)
       call this%init_self_deriv_scheme(dpts_in)


       this%alat = fdf_physical('LatticeConstant',0.0_dp,'Bohr')
       if (this%alat==0.0_dp) call die('init_thermal_flux', &
            &'ThermalFlux requires LatticeConstant set!')
       write(*,"(2X,A,F10.4,A)") "[gk:.init] Lattice Constant parameter: ", this%alat, " Bohr"

       this%meshcutoff = fdf_physical('MeshCutoff',0.0_dp,'Ry')
       if (this%meshcutoff==0.0_dp) call die('zero_flux:init', &
            &'ThermalFlux requires meshcutoff set')
       write(*,"(2X,A,F10.4,A)") "[gk:.init] Mesh Cutoff parameter: ", this%meshcutoff, " Ry"

       this%virtual_dt = fdf_get('ThermalFlux.Virtual.dt',0.1_dp,'fs')
       write(*,"(2X,A,F10.4,A)") "[gk:.init] Substep delta_t: ", this%virtual_dt, " fs"

       this%eta_ewald = fdf_double("ThermalFlux.Jion.Eta", 0.1_dp)
       write(*,"(2X,A,F10.4,A)") "[gk:.init] Ewald factor for convergence (eta): ", this%eta_ewald

       this%n_max_ewald = fdf_integer("ThermalFlux.Jion.Nmax", 5)
       write(*,"(2X,A,I2)") "[gk:.init] Number of periodic cell images for Ewald scheme: ", this%n_max_ewald

       if(fdf_get('ThermalFlux.QEHeat.Units', .false.)) then
          this%qeheat_units = .true.
          write(*,"(2X,A)") "[gk:.init] Using QEHeat units for thermal flux results output."
       end if

       if(fdf_get('ThermalFlux.UseSternheimer', .false.)) then
          this%use_sternheimer = .true.
          write(*,"(2X,A)") "[gk:.init] Sternheimer Equation method for computation of Jks."
       end if

       if(fdf_get('ThermalFlux.VerboseOutput', .false.)) then
          this%verbose_output = .true.
          write(*,"(2X,A)") "[gk:.init] Verbose output requested."
       end if

       this%init = .true.
       write(*,"(2X,A)") "[gk:.init] Done with thermal flux init."
       write(*,"(2X,A)") "=============================================================="
    end if
  end subroutine init_thermal_flux_settings


  subroutine init_thermal_flux_results(this)
    !! Init all results to zero
    class(thermal_flux_results_type), intent(inout) :: this

    this%Jtotal(:) = 0.0_dp

    this%Jks(:)    = 0.0_dp
    this%Jks_A(:)  = 0.0_dp
    this%Jks_B(:)  = 0.0_dp
    this%Jele(:)   = 0.0_dp

    this%Jxc(:)    = 0.0_dp

    this%Jhart(:)  = 0.0_dp

    this%Jion(:)   = 0.0_dp
    this%Jion_A(:) = 0.0_dp
    this%Jion_B(:) = 0.0_dp
    this%Jion_C(:) = 0.0_dp
    this%Jion_D(:) = 0.0_dp
    this%Jion_E(:) = 0.0_dp

    this%Jzero(:)  = 0.0_dp
    this%Jzloc(:)  = 0.0_dp
    this%Jznl(:)   = 0.0_dp
  end subroutine init_thermal_flux_results


  subroutine write_thermal_flux_results(this, setup)
    class(thermal_flux_results_type),  intent(inout) :: this
    class(thermal_flux_settings_type), intent(inout) :: setup

    real(dp), parameter :: qeheat_factor = 0.0483776900146_dp
    real(dp)            :: scale_factor  = 1.0_dp

    if (setup%qeheat_units) scale_factor = qeheat_factor

    write(*,*)
    write(*,*) "#           ================ Green-Kubo ThermalFlux Results ================"

    write(*,*) "[gk:.Jks...]", this%Jks * scale_factor
    if (setup%verbose_output) then
       write(*,*) "[gk:.Jks-A.]", this%Jks_A * scale_factor
       write(*,*) "[gk:.Jks-B.]", this%Jks_B * scale_factor
       write(*,*) "[gk:.Jele..]", this%Jele * scale_factor
    end if

    write(*,*) "[gk:.Jxc...]", this%Jxc * scale_factor

    write(*,*) "[gk:.Jhart.]", this%Jhart * scale_factor

    write(*,*) "[gk:.Jion..]", this%Jion * scale_factor
    if (setup%verbose_output) then
       write(*,*) "[gk:.Jion-A]", this%Jion_A * scale_factor
       write(*,*) "[gk:.Jion-B]", this%Jion_B * scale_factor
       write(*,*) "[gk:.Jion-C]", this%Jion_C * scale_factor
       write(*,*) "[gk:.Jion-D]", this%Jion_D * scale_factor
       write(*,*) "[gk:.Jion-E]", this%Jion_E * scale_factor
    end if

    write(*,*) "[gk:.Jzero.]", this%Jzero * scale_factor
    if (setup%verbose_output) then
       write(*,*) "[gk:.Jzloc.]", this%Jzloc * scale_factor
       write(*,*) "[gk:.Jznl..]", this%Jznl * scale_factor
    end if

    write(*,*) "#-----------"

    write(*,*) "[gk:.Jtotal]", this%Jtotal * scale_factor

    if (setup%qeheat_units) then
       write(*,*) "#           ======================    QEHeat Units    ======================"
    else
       write(*,*) "#           ======================    Siesta Units    ======================"
    end if
    write(*,*)
  end subroutine write_thermal_flux_results

end module m_thermal_flux_settings


module thermal_flux_data
  !! NOTE: Not memory-wise optimal;
  !! look out for `^mem' - probably reduce-able entities.

  use precision, only: dp, grid_p
  use m_thermal_flux_settings

  type(thermal_flux_settings_type) :: gk_setup
  !! Short distinct name of the thermal
  !! transport settings' singleton example.

  type(thermal_flux_results_type) :: gk_results
  !! Object to store the results of the
  !! Thermal Transport calculation.

  integer :: substep !! Current substep number
  !! (will go from 1 to `gk_setup%dpoints` in runtime)

  real(dp), allocatable, save :: xa_in(:,:) !! system coordinates input for the snapshot
  real(dp), allocatable, save :: va_in(:,:) !! system velocities input for the snapshot
  real(dp), allocatable, save :: v_cm(:,:)  !! center of mass velocities for each atomic kind
  integer, allocatable, save  :: nasp(:) !! number of atoms of a certain kind

  real(dp), allocatable, save :: DM_save(:,:,:,:) !! Full density matrix storage +1 dim:
  !! Last dimension goes over substep points
  !! (the derivation scheme parameter)
  !!
  !! As for all functions meant to be differentiated,
  !! after call to `compute_derivatives()'
  !! the last rank index `1' stores the "base" step value
  !! and the last rank index `2' will hold the derivative.

  real(dp), allocatable, save :: gradS_base(:,:)
  !! Gradient of overlap copy from "base" step.
  real(dp), allocatable, save :: Rmat_base(:,:)
  !! Rmatrix copy from "base" step.
  real(dp), allocatable, save   :: hr_commutator(:,:)
  !! [H,r] computed at base step

  real(dp), pointer, save :: Sinv(:,:) => null()
  !! Inverse of the overlap matrix.
  !! Computed in diagg at the beginning of the scf cycle of
  !! the base step. Possibility of memory optimization (^mem).
  real(dp), pointer, save :: S_base(:,:) => null() ! ^mem
  real(dp), pointer, save :: H_base(:,:) => null() ! ^mem
  real(dp), pointer, save :: eo_base(:) => null()
  real(dp), pointer, save :: psi_base(:,:) => null()

  ! Common auxilliary state data regarding mesh, gvectors etc:
  real(dp), dimension(3,3), save :: cell_vmd !! Unit cell vectors
  integer, dimension(3), save    :: mesh_vmd
  !! Number of global mesh divisions.
  !! Used in particular for `ngm_vmd` limits estimation.
  integer, dimension(3), save    :: ntml_vmd
  !! Total number of mesh points stored locally
  integer, save                  :: nsm_vmd
  !! Number of sub-mesh points per mesh point along each axis
  real(dp), save :: alat_vmd !! Lattice constant, retrieved in Bohr
  real(dp), save :: at_vmd(3,3) !! Lattice vectors scaled by alat
  real(dp), save :: Vna_G0
  !! Vna(G=0) obtained as a sum(Vna)*dvol/volume
  !! on a real space grid.

  real(grid_p), pointer, save :: charge_g(:,:) => null()
  !! charge density in the reciprocal space.
  !! Used to get Hartree potential and for computation of
  !! local part of the zero flux during execution of VMD logic.
  real(grid_p), allocatable, target, save :: charge_g_base(:,:) ! ^mem
  !! Charge density in the reciprocal space at `base' substep.

  integer, save :: gstart_vmd = 2
  !! index of the first G vector whose module is > 0
  !! In QE it's used for parallel execution, but not only:
  !! it also omits G=0 for loops in special cases where the component
  !! in question diverges at G=0.
  !! (when parallelized should work like this:
  !!  gstart=2 for the proc that holds G=0, gstart=1 for all others.)
  integer, save :: ngm_plus_vmd = 0 !! number of G vectors in G>
  !!TODO: (when parallelized: local, on this processor)
  integer, save :: ngm_max_vmd = 0 !! Maximum possible number of G vectors
  integer, save :: ngl_vmd = 0  !! number of G-vector shells
  real(dp), allocatable, target :: g_vmd(:,:)
  !! G-vectors cartesian components
  !! (for QE it's in units tpiba =(2pi/a)  )
  integer, allocatable :: igsrt_vmd(:)
  !! `igsrt_vmd(ngm_max_vmd)` : fft indexes for g-vectors,
  !! sorted in ascending order.
  integer, allocatable :: igplus_vmd(:)
  !! `igplus_vmd(ngm_plus_vmd)` : fft indexes for g-vectors in G>,
  !! sorted in ascending order.
  !! (should correspond to QE's `nl` array for G>)
  real(dp), allocatable, target :: gg_vmd(:)
  !! G^2 (for QE it's in units of tpiba2=(2pi/a)^2)
  real(dp), allocatable, target :: gl_vmd(:)
  !! gl(i) = i-th shell of G^2 (in units of tpiba2)
  integer, allocatable, target :: igtongl_vmd(:)
  !! shell index for n-th G-vector
  ! integer, allocatable, target :: mill_vmd(:,:) !DEBUG
  !! mill = miller index of G vectors

  real(grid_p), allocatable, save :: Rho_save(:,:) ! ^mem
  !! Array to store derivative of local charge Rho.
  !!
  !! NOTE: It is so far one-dimensional:
  !! only first spin component.
  real(grid_p), allocatable, save :: thtr_dexcdGD(:,:,:)
  !! `dexcdGD` from the call to cellXC (LibGridXC)
  real(dp), save :: xc_flux_dvol
  !! Element of volume size saved from `dhscf`

  real(grid_p), allocatable, target, save :: Vhart_save(:,:,:)
  !! Obtained from charge density in the reciprocal space for
  !! computation of its derivative for the Hartree flux component.

  real(dp) :: I_prime, I_prime_rec
  real(dp), allocatable :: I_first_g(:,:,:)
  real(dp), allocatable :: I_second_g(:)

contains

    subroutine calculate_cm_velocities()
      !! Calculate center of mass velocities for each atomic type.
      !!
      !! Subtraction is not performed here (as is not in QEHeat by default),
      !! rather the J_cm[sp][1-3] outputted for post-processing (decorrelation).

      use atm_types,   only: species, species_info
      use basis_types, only: nsp
      use siesta_geom, only: isa, na_u

      integer :: iatom, itype
      real(dp) :: delta(3), mean(3)

      allocate (v_cm(3, nsp))
      allocate (nasp(nsp))
      nasp = 0
      v_cm = 0.0_dp

      do iatom = 1, na_u
         itype = isa(iatom)
         nasp(itype) = nasp(itype) + 1
         delta = (va_in(:, iatom) - v_cm(:, itype))/real(nasp(itype), dp)
         v_cm(:, itype) = v_cm(:, itype) + delta
      end do

  end subroutine calculate_cm_velocities


  subroutine write_cm_velocities()
    use basis_types, only: nsp
    integer :: itype
    real(dp), parameter :: qeheat_factor = 0.0483776900146_dp
    real(dp)            :: scale_factor  = 1.0_dp

    if (gk_setup%qeheat_units) scale_factor = qeheat_factor

    write (*,*) "#           ====================== Center of mass currents/ atom type:"
    do itype = 1, nsp
       write (*, '(A,I0.3,A,3E20.12)') "[gk:.Jcm", itype, "]", real(nasp(itype), dp) * v_cm(:, itype) * scale_factor
    end do

  end subroutine write_cm_velocities

end module thermal_flux_data
