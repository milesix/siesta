module thermal_flux

  use precision, only: dp, grid_p
  use siesta_geom, only: xa, va, na_u

  use thermal_flux_data
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

    do ia=1,na_u
       print *, ia, xa_in(1:3, ia), va_in(1:3, ia)
    end do

    close (fu)
  end subroutine read_xvs


  subroutine reset_thermal_flux()
    deallocate(xa_in)
    deallocate(va_in)
  end subroutine reset_thermal_flux

end module thermal_flux
