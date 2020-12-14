module ks_flux_data

  use class_Sparsity
  use precision, only: dp, grid_p

  real(dp), allocatable       :: psi_hat_c(:,:,:)
  real(dp), allocatable       :: psi_dot_c(:,:)

  type(Sparsity), public      :: base_step_sparse_pattern
  real(dp), pointer, save     :: Dscf_deriv(:,:)
  !! Array to store derivative of density matrix
  !! computed within virtual-step MD program flow.
  !!
  !!NOTE: It is so far one-dimensional:
  !! only first spin component.
  real(dp), dimension(:), target, allocatable :: ks_flux_S
  !! Storage for overlap from the `BASE` step.
  real(dp), dimension(:,:), target, allocatable :: ks_flux_gS
  !! Storage for overlap gradient from the `BASE` step.
  real(dp), dimension(:,:), target, allocatable :: ks_flux_H
  !! Storage for the Hamiltonian from the `BASE` step.
  real(dp), dimension(:,:), target, allocatable :: ks_flux_D
  !! Storage for density matrix from the `BASE` step.
  real(dp), dimension(:,:), target, allocatable :: ks_flux_Rmat
  !! Storage for R-matrix from the `BASE` step.

  real(dp), save  :: ks_flux_Jks(3)
  !! Kohn-Sham component of the heat flux.

  real(dp), save  :: ks_flux_Jele(3)
  !! Electronic flux.

contains

  subroutine init_ks_flux_data ()
    use sparse_matrices, only: H, S, gradS, Dscf, Rmat
    use alloc, only : re_alloc

    call re_alloc( Dscf_deriv,&
         & 1, size(Dscf,1), 1, size(Dscf,2),&
         & name='Dscf_deriv', routine='init_ks_flux_data' )

    Dscf_deriv(:,:) = Dscf(:,:)

    allocate(ks_flux_S, source=S)
    allocate(ks_flux_gS, source=gradS)
    allocate(ks_flux_H, source=H)
    allocate(ks_flux_D, source=Dscf)
    allocate(ks_flux_Rmat, source=Rmat)
  end subroutine init_ks_flux_data

  subroutine reset_ks_flux_data ()
    use alloc, only : de_alloc

    call de_alloc(Dscf_deriv, "Dscf_deriv", "virtual_step_md")
    nullify(Dscf_deriv)

    deallocate(ks_flux_S)
    deallocate(ks_flux_gS)
    deallocate(ks_flux_H)
    deallocate(ks_flux_D)
    deallocate(ks_flux_Rmat)
  end subroutine reset_ks_flux_data

end module ks_flux_data
