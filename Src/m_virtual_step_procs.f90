module m_virtual_step_procs
  use precision, only: dp, grid_p

  use siesta_options, only: virtual_dt
  use siesta_options, only: virtual_md_Jhart
  use siesta_options, only: virtual_md_Jks
  use siesta_options, only: virtual_md_Jxc
  use siesta_options, only: virtual_md_Jion
  use siesta_options, only: virtual_md_Jzero
  use siesta_options, only: want_virtual_step_md
  use siesta_options, only: want_vmd_in_dhscf

  use gvecs, only: gvecs_init, gvecs_teardown, gvecs_gen, gshells

  implicit none

contains

  subroutine base_md_step_logic ()
    !! Wrapper over logic performed in the `Base` MD step,
    !! e.g. output integrated energy components etc.

    use sparse_matrices, only: Dscf

    use ks_flux_data, only: Dscf_deriv
    use ks_flux_data, only: init_ks_flux_data
    use xc_flux_data, only: init_xc_flux_data
    use ion_flux_data,only: init_ion_flux_data
    use ion_flux_data, only: ion_flux_mesh, ion_flux_cell

    if ( virtual_md_Jks ) then
       call init_ks_flux_data()

       Dscf_deriv(:,:) = Dscf(:,:)
    end if

    if ( virtual_md_Jxc ) call init_xc_flux_data()

    if ( virtual_md_Jion .or. virtual_md_Jzero ) then
       ! Initialize and generate G-vectors data
       call gvecs_init(ion_flux_mesh)
       call gvecs_gen(ion_flux_cell, ion_flux_mesh)
       call gshells(.false.)    ! init G-vector shells assuming non-variable cell
    end if

    if ( virtual_md_Jion ) call init_ion_flux_data()

  end subroutine base_md_step_logic


  subroutine virtual_md_step_logic ()
    !! Wrapper over logic performed in the `Virtual` MD step,
    !! e.g. computation of derivatives of heat flux components etc.

    use sparse_matrices, only: Dscf, gradS

    use hartree_flux_data, only: h_flux_Jhart

    use ks_flux_data, only: Dscf_deriv
    use ks_flux_procs, only: compute_Jks

    use xc_flux_data, only: thtr_Rho, thtr_Rho_deriv, thtr_dexcdGD
    use xc_flux_data, only: xc_flux_Jxc, xc_flux_dvol

    use ion_flux_data, only: ion_flux_a, ion_flux_b
    use ion_flux_data, only: ion_flux_c, ion_flux_d, ion_flux_e
    use ion_flux_data, only: ion_flux_Jion
    use ion_flux_procs, only: compute_Jion_a, compute_Jion_b
    use ion_flux_procs, only: compute_Jion_cde

    use zero_flux_data,  only: zero_flux_Jzero
    use zero_flux_procs, only: compute_Jzero

    use ks_flux_data,      only: reset_ks_flux_data
    use xc_flux_data,      only: reset_xc_flux_data
    use ion_flux_data,     only: reset_ion_flux_data
    use hartree_flux_data, only: reset_hartree_flux_data

    ! use ion_flux_data, only: ion_flux_mesh, ion_flux_cell

    integer :: i

    if ( virtual_md_Jhart ) then
       print*, "[Jhart] ", h_flux_Jhart(:)
    end if

    if ( virtual_md_Jks ) then
       Dscf_deriv(:,:) = (Dscf_deriv(:,:)-Dscf(:,:))/virtual_dt

       call compute_Jks()

       ! call reset_ks_flux_data()
    end if

    if ( virtual_md_Jxc ) then
       ! compute Rho derivative on a virtual step:
       thtr_Rho_deriv(:) = (thtr_Rho_deriv(:)-thtr_Rho(:))/virtual_dt

       ! xc_flux_Jxc(1:3) = 0.0_dp <- CHECK: should be reset to 0 already from `init_xc_flux_data`

       do i=1,size(thtr_Rho)
          xc_flux_Jxc(1:3) = xc_flux_Jxc(1:3) &
               & - thtr_Rho_deriv(i) * thtr_dexcdGD(i,1:3,1) &
               & * xc_flux_dvol
       end do

       print*, "[Jxc] ", xc_flux_Jxc(:)

    end if

    if ( virtual_md_Jion ) then

       call compute_Jion_a()

       print*, "[Jion] flux A: ", ion_flux_a(:)

       call compute_Jion_b()

       print*, "[Jion] flux B: ", ion_flux_b(:)

       call compute_Jion_cde()

       print*, "[Jion] flux C: ", ion_flux_c(:)
       print*, "[Jion] flux D: ", ion_flux_d(:)
       print*, "[Jion] flux E: ", ion_flux_e(:)

       ! call reset_ion_flux_data()
    end if

    if ( virtual_md_Jzero ) then
       call compute_Jzero()
    end if

    if ( virtual_md_Jion .or. virtual_md_Jzero ) then
       ! tearing down G-vectors data
       call gvecs_teardown()
    end if

    !FIXME: regroup resetters in corresponding subroutine
    if ( virtual_md_Jion) call reset_ion_flux_data()
    if ( want_vmd_in_dhscf ) call reset_xc_flux_data()
    if ( virtual_md_Jks ) call reset_ks_flux_data()

    ! if ( want_vmd_in_dhscf ) call reset_xc_flux_data()

  end subroutine virtual_md_step_logic


  subroutine reset_virtual_step_md ()
    !! Deallocates auxiliary arrays for before- and after-Base
    !! MD move system state.

    use class_Sparsity
    use alloc, only : de_alloc
    use m_virtual_step_data
    use ks_flux_data, only: base_step_sparse_pattern

    use ks_flux_data,      only: reset_ks_flux_data
    use xc_flux_data,      only: reset_xc_flux_data
    use ion_flux_data,     only: reset_ion_flux_data
    use hartree_flux_data, only: reset_hartree_flux_data

    if ( want_virtual_step_md ) then !TODO: make this check in `siesta_end`

       if (virtual_md_Jhart) call reset_hartree_flux_data()
       ! if ( virtual_md_Jion) call reset_ion_flux_data()
       ! if ( want_vmd_in_dhscf ) call reset_xc_flux_data()
       ! if ( virtual_md_Jks ) call reset_ks_flux_data()

       ! clean the aux `Base`-step sparse pattern
       call delete( base_step_sparse_pattern )

       call de_alloc(fa_before_move, "fa_before_move", "virtual_step_md")
       call de_alloc(xa_before_move, "xa_before_move", "virtual_step_md")
       call de_alloc(va_before_move, "va_before_move", "virtual_step_md")

       call de_alloc(fa_after_move, "fa_after_move", "virtual_step_md")
       call de_alloc(xa_after_move, "xa_after_move", "virtual_step_md")
       call de_alloc(va_after_move, "va_after_move", "virtual_step_md")

       nullify(fa_before_move)
       nullify(xa_before_move)
       nullify(va_before_move)

       nullify(fa_after_move)
       nullify(xa_after_move)
       nullify(va_after_move)
    end if

  end subroutine reset_virtual_step_md

end module m_virtual_step_procs
