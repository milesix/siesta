module thermal_flux_data

  use precision, only: dp, grid_p

  real(dp), allocatable :: xa_in(:,:) !! system coordinates input for the snapshot
  real(dp), allocatable :: va_in(:,:) !! system velocities input for the snapshot

  integer, save :: vs_direction !! direction of the "virtual" step

  real(dp), allocatable :: Dfull_1st(:,:,:) !! Full density matrix for derivation at 1st point
  real(dp), allocatable :: Dfull_2nd(:,:,:) !! Full density matrix for derivation at 2nd point
  !! Will contain full DM derivative after finite difference

end module thermal_flux_data
