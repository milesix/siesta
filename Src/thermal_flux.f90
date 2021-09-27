module thermal_flux

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
    integer  :: step_interval

    step_interval = gk_setup%func_step_intervals(substep)

    xa(:,:) = xa_in(:,:) + step_interval * va_in(:,:) * gk_setup%virtual_dt

  end subroutine make_small_move


  subroutine compute_derivatives()
    ! integer :: k, j, i          !debug

    call apply_derivation(DM_save, gk_setup%virtual_dt, gk_setup%dpoints)

    ! debug
    ! do k = 1,size(Dfull, 3)
    !    do j = 1,size(Dfull,2)
    !       do i = 1,size(Dfull,1)
    !          write((100+gk_setup%dpoints),*) Dfull(i,j,k,2)
    !       end do
    !    end do
    ! end do
    ! debug

  end subroutine compute_derivatives


  subroutine reset_thermal_flux()
    deallocate(xa_in)
    deallocate(va_in)
    deallocate(DM_save)
  end subroutine reset_thermal_flux

end module thermal_flux


module derivation_routines
  use precision, only: dp, grid_p
  implicit none

contains

  subroutine apply_derivation(F, h, scheme)
    real(dp), intent(inout) :: F(:,:,:,:)
    real(dp), intent(in) :: h
    integer, intent(in)  :: scheme

    select case(scheme)
    case(2)
       call derivation_2pts(F, h)
    case(3)
       call derivation_3pts_mid(F, h)
    case(5)
       call derivation_5pts_mid(F, h)
    case default
       call die('Non-existent derivation scheme requested, aborrting.')
    end select
  end subroutine apply_derivation

  subroutine derivation_2pts(F, h)
    real(dp), intent(inout) :: F(:,:,:,:)
    real(dp), intent(in) :: h

    F(:,:,:,2) = (F(:,:,:,2) - F(:,:,:,1)) / h
  end subroutine derivation_2pts

  subroutine derivation_3pts_mid(F, h)
    real(dp), intent(inout) :: F(:,:,:,:)
    real(dp), intent(in) :: h

    F(:,:,:,2) = (F(:,:,:,3) - F(:,:,:,2)) / (2.0_dp * h)
  end subroutine derivation_3pts_mid

  subroutine derivation_5pts_mid(F, h)
    real(dp), intent(inout) :: F(:,:,:,:)
    real(dp), intent(in) :: h

    F(:,:,:,2) = (F(:,:,:,2) - 8.0_dp * F(:,:,:,3) &
         & + 8.0_dp * F(:,:,:,4) - F(:,:,:,5)) / (12.0_dp * h)
  end subroutine derivation_5pts_mid

end module derivation_routines
