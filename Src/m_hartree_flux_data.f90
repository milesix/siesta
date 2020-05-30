module hartree_flux_data
  use precision, only: dp, grid_p

  real(grid_p), pointer, save :: Vhart_deriv(:,:) => null()
  real(grid_p), public, save  :: h_flux_Jhart(3) = 0.0

contains

  ! subroutine init_hartree_flux_data ()
  !   call re_alloc(Vhart_deriv, 1, size(Vhart)&
  !        & name='Vhart_deriv', routine='init_hartree_flux_data')

  !   Vhart_deriv(:) = Vhart(:)
  ! end subroutine init_hartree_flux_data

  !TODO: put in the teardown section of SIESTA
  subroutine reset_hartree_flux_data ()
    use alloc, only : de_alloc

    ! call de_alloc(Vhart_rec, "Vhart_rec", "virtual_step_md")
    ! nullify(Vhart)
    call de_alloc(Vhart_deriv, "Vhart_deriv", "virtual_step_md")
    nullify(Vhart_deriv)
  end subroutine reset_hartree_flux_data
end module hartree_flux_data
