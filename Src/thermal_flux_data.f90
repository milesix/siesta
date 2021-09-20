module thermal_flux_data

  use precision, only: dp, grid_p

  real(dp), allocatable :: xa_in(:,:) !! system coordinates input for the snapshot
  real(dp), allocatable :: va_in(:,:) !! system velocities input for the snapshot

end module thermal_flux_data
