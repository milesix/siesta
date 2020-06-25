module zero_flux_data

  use precision, only: dp

  real(dp), save  :: zero_flux_Jzero(3)
  !! Zero-flux component of the heat flux.

contains

  ! subroutine init_zero_flux_data ()
  !   use sparse_matrices, only: S, gradS, Dscf, Rmat
  !   use alloc, only : re_alloc

  !   call re_alloc( Dscf_deriv,&
  !        & 1, size(Dscf,1), 1, size(Dscf,2),&
  !        & name='Dscf_deriv', routine='init_zero_flux_data' )

  !   Dscf_deriv(:,:) = Dscf(:,:)

  ! end subroutine init_zero_flux_data

  ! subroutine reset_zero_flux_data ()
  !   use alloc, only : de_alloc

  ! end subroutine reset_zero_flux_data

end module zero_flux_data
