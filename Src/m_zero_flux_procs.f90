module zero_flux_procs
  !! Contains procedures for computation of J_zero.
  !! See notes on components of the Kohn-Sham flux for thermal transport.
  use precision, only: dp
  ! use atomlist, only: no_u, no_l, indxuo, iaorb, qtot
  use siesta_geom, only: na_u
  use atomlist, only: no_u
  use sparse_matrices, only: listh, listhptr, numh

  use zero_flux_data, only: zero_flux_Jzero

  !FIXME: not all of them are needed:
  use matel_mod, only: init_matel_main_tables
  use matel_mod, only: init_matel_thermal_transport
  use matel_mod, only: init_matel_orb_XYZ_orb
  use matel_mod, only: init_matel_optical_P
  use matel_mod, only: new_matel

  implicit none

  character(len=1), dimension(3), parameter :: coord_table = [ 'X', 'Y', 'Z' ]
  integer :: mu, nu, I, alp

  private
  public :: compute_Jzero

contains

subroutine compute_Jzero ()
  use ks_flux_data, only: ks_flux_D

  real(dp), dimension(:),   pointer :: dscf => null()

  call init_matel_main_tables()
  call init_matel_thermal_transport()
  call init_matel_orb_XYZ_orb()
  call init_matel_optical_P()

  dscf => ks_flux_D(:,1)        !NOTE: only 1st spin component

  ! This should add NL-part of the zero current
  ! do mu=1,no_u                  ! ?
  !    do nu=1,no_u               ! ?
  !       do I=1,na_u             ! ?
  !          do alp=1,3

  !          enddo
  !       enddo
  !    enddo
  ! enddo

end subroutine compute_Jzero

function SG(i) result (s)
  integer, intent(in) :: i
  character (len=2)   :: s
  s = 'G' // coord_table(i)
end function SG

function RG(i,j) result (s)
  integer, intent(in) :: i, j
  character (len=3)   :: s
  s = 'H' // coord_table(i) // coord_table(j)
end function RG

end module zero_flux_procs
