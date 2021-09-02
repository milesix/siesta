module zero_flux_data

  use precision, only: dp

  real(dp), save  :: zero_flux_Jzero(3)
  !! Zero-flux component of the heat flux.
  real(DP), allocatable, save :: tab_local(:,:,:)
  !! One-dimensional table for radial part.
  real(dp), allocatable, save :: H_g(:,:,:,:)
  !! Table for the values of `h_ab` used in computation
  !! of the local, long-range part of the Zero heat current.
  complex(DP), allocatable, save :: u_g(:,:)

  real(dp), parameter :: dq = 0.01_dp   !! space between points in the pseudopotential tab.
  real(dp), parameter :: cell_factor = 1.0_dp
  !! maximum expected (linear) cell contraction
  !! during relaxation/MD

  real(dp)  :: pi   = 4.0_dp * atan(1.0_dp)
  real(dp)  :: alat, tpiba, meshcutoff, pref
  real(dp)  :: gcutm  !! potential cut-off
  integer   :: nqxq   !! size of interpolation table


contains

  subroutine init_zero_flux_data ()
    use fdf, only: fdf_physical
    use utils, only: die
    use cellsubs, only : volcel
    use m_virtual_step_data, only: cell_vmd
    !NOTE: This should be present even with `ion_flux` calculation switched off,
    ! though ofc it should be moved into a proper scope.

    real(dp) :: omega

    alat = fdf_physical('LatticeConstant',0.0_dp,'Bohr')
    if (alat==0.0_dp) call die('zero_flux:init', 'VMD Jzero requires alat set')

    omega = volcel(cell_vmd)

    tpiba = 2.0_dp * pi / alat
    pref  = 4.0_dp * pi / omega

    meshcutoff = fdf_physical('MeshCutoff',0.0_dp,'Ry')
    if (meshcutoff==0.0_dp) call die('zero_flux:init', 'VMD Jzero requires meshcutoff set')

    ! gcutm = meshcutoff / tpiba**2
    nqxq  = int((sqrt(meshcutoff)/dq + 4) * cell_factor)

    zero_flux_Jzero(:) = 0.0_dp ! reset Jzero for this step

  end subroutine init_zero_flux_data

  ! subroutine reset_zero_flux_data ()
  !   use alloc, only : de_alloc

  ! end subroutine reset_zero_flux_data

end module zero_flux_data
