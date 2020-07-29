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
  ! integer :: mu, nu, I, alp

  private
  public :: compute_Jzero

contains

subroutine compute_Jzero ()
  use ks_flux_data, only: ks_flux_D

  use atm_types,   only: species, species_info
  use basis_types, only: nsp
  use radial,      only: rad_get
  use m_bessph,    only: BESSPH

  type(species_info), pointer :: spp
  real(dp) :: rm, delta_rm, aux_fr, aux_dfdr, vl_int
  integer  :: in, is, npts
  real(dp), allocatable :: rvals(:), aux(:)

  rm = 0.0_dp

  ! find rm and npts for vlocal
  do is = 1, nsp                !
     spp => species(is)
     if (spp%reduced_vlocal%cutoff > rm) then
        rm = spp%reduced_vlocal%cutoff
        npts = spp%reduced_vlocal%n
        delta_rm = spp%reduced_vlocal%delta
     end if
  end do

  allocate(rvals(npts))
  allocate(aux(npts))

  rvals(:) = 0.0_dp
  aux(:) = 0.0_dp

  do in=1,npts                   ! init R values for selected* radial function
     rvals(in) = delta_rm * (in-1)
  enddo

  ! aux (ir) = (rgrid(nt)%r(ir)*upf(nt)%vloc(ir)+2.d0*zv(nt))* besr (ir) * rgrid(nt)%r(ir)
  do is = 1,nsp
     spp => species(is)
     do in=1,npts
        call rad_get(spp%reduced_vlocal, rvals(in), aux_fr, aux_dfdr)
        aux(in) = aux_fr * bessph(0, rvals(in)) * rvals(in) ! bessel_spherical for l=0
     end do
  end do

  call simpson(npts, delta_rm, aux, vl_int)
  print*, "[Jzero] integral of aux", vl_int

  ! cleanup
  deallocate(rvals)
  deallocate(aux)

  ! real(dp), dimension(:),   pointer :: dscf => null()

  ! call init_matel_main_tables()
  ! call init_matel_thermal_transport()
  ! call init_matel_orb_XYZ_orb()
  ! call init_matel_optical_P()

  ! dscf => ks_flux_D(:,1)        !NOTE: only 1st spin component

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

subroutine simpson(mesh, h, func, res)
  !! This Simpson's integration rule is short and dirty.
  !! Constant `h` provided instead of metrics for linear mesh.
  integer,  intent(in) :: mesh
  real(dp), intent(in) :: h, func(mesh)
  real(dp), intent(out):: res
  real(dp) :: f1, f2, f3, r12
  integer  :: i

  res = 0.0d0
  r12 = h * 1.0d0 / 3.0d0
  f3 = func (1) * r12

  do i = 2, mesh - 1, 2
     f1 = f3
     f2 = func (i) * r12
     f3 = func (i + 1) * r12
     res = res + f1 + 4.0d0 * f2 + f3
  end do
end subroutine simpson

end module zero_flux_procs
