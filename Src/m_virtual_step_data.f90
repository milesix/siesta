module m_virtual_step_data
  !! Module that holds state data in between `Base` and `Virtual` MD steps.
  !! State of the system includes atom positions, velocities and forces
  !! on atoms, as well as mesh information and G-vectors in a way similar
  !! to gvecs handling in QE.

  use precision, only: dp, grid_p

  implicit none

  public

  integer, parameter :: BASE_STEP = 0 !! Designates base step in an SCF cycle
  integer, parameter :: VIRTUAL_STEP = 1
  !! Designates `virtual` step: a displacement over `virtual_dt` time step
  !! that we use to get derivatives as components to e.g. heat flux.
  integer :: substep = 0 !! Current substep within SCF cycle
  integer, save :: istep_vmd   !! Global step-index copy

  real(dp), save  :: Jtotal(3)
  !! Total value of the heat flux calculated in terms of the G-K formalism.

  real(dp), pointer, save :: Dfull(:,:,:) => null()
  !! Full density matrix where we need it full (those cases
  !! where we cannot account for sparsity in strides over orbitals).
  !! Processed in the following modules/subroutines:
  !! - `diagg()`
  !! - `siesta_forces()`
  real(dp), pointer, save :: Dderiv(:,:,:) => null()
  !! Time derivative of the `Dfull`.
  !! Serves as accumulator, storing full DM value in each VMD substep
  !! (NOTE: now, in every scf cycle. That should be optimized) in `diagg()`.
  !! Then in the upper `siesta_forces`, when `SCFconverged`:
  !! - on "Base" step it is copied to Dfull;
  !! - on "Virtual" it is turned into its derivative.
  !!
  !! Processed in the following modules/subroutines:
  !! - `diagg()`
  !! - `siesta_forces()`

  real(dp), pointer, save :: Sinv(:,:) => null()
  !! Inverse of the overlap matrix.
  !! Computed in diagg at the beginning of the scf cycle of
  !! the base step
  real(dp), pointer, save :: S_base(:,:) => null()
  real(dp), pointer, save :: H_base(:,:) => null()
  real(dp), pointer, save :: eo_base(:) => null()

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

  real(grid_p), pointer, save :: charge_g(:,:) => null()
  !! charge density in the reciprocal space.
  !! Used to get Hartree potential and for computation of
  !! local part of the zero flux during execution of VMD logic.

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

  !NOTE: Should thtr_Rho be defined here, since it's sort of global too?
  ! Atm it resides in `xc_flux_data`.

  ! Forces, coordinates and velocities before and after `base` move:
  real(dp), pointer, save :: fa_before_move(:, :) => null()
  real(dp), pointer, save :: xa_before_move(:, :) => null()
  real(dp), pointer, save :: va_before_move(:, :) => null()

  real(dp), pointer, save :: fa_after_move(:, :) => null()
  real(dp), pointer, save :: xa_after_move(:, :) => null()
  real(dp), pointer, save :: va_after_move(:, :) => null()

end module m_virtual_step_data
