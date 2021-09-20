module thermal_flux

  use precision, only: dp, grid_p
  use siesta_geom, only: xa, va, na_u
  use siesta_options, only: virtual_dt

  use thermal_flux_data
  use m_siesta_forces
  external :: die

contains

  subroutine init_thermal_flux()
    allocate(xa_in(3, na_u))
    allocate(va_in(3, na_u))

    call read_xvs()             ! will fail if coordinates not specified
  end subroutine init_thermal_flux


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

    xa(:,:) = xa_in(:,:) + vs_direction * 0.5_dp * va_in(:,:) * virtual_dt

  end subroutine make_small_move


  subroutine reset_thermal_flux()
    deallocate(xa_in)
    deallocate(va_in)
    deallocate(Dfull_1st)
    deallocate(Dfull_2nd)
  end subroutine reset_thermal_flux

end module thermal_flux
