module m_virtual_step_data
  !! Module that holds atom positions, velocities and forces
  !! on atoms in between `Base` and `Virtual` MD steps.

  use precision, only: dp

  implicit none

  public

  integer, parameter :: BASE_STEP = 0 !! Designates base step in an SCF cycle
  integer, parameter :: VIRTUAL_STEP = 1
  !! Designates `virtual` step: a displacement over `virtual_dt` time step
  !! that we use to get derivatives as components to e.g. heat flux.
  integer :: substep = 0 !! Current substep within SCF cycle
  integer, save :: istep_vmd   !! Global step-index copy

  ! Forces, coordinates and velocities before and after `base` move:
  real(dp), pointer, save :: fa_before_move(:, :) => null()
  real(dp), pointer, save :: xa_before_move(:, :) => null()
  real(dp), pointer, save :: va_before_move(:, :) => null()

  real(dp), pointer, save :: fa_after_move(:, :) => null()
  real(dp), pointer, save :: xa_after_move(:, :) => null()
  real(dp), pointer, save :: va_after_move(:, :) => null()

end module m_virtual_step_data
