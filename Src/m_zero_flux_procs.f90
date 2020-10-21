module zero_flux_procs
  !! Contains procedures for computation of J_zero.
  !! See notes on components of the Kohn-Sham flux for thermal transport.
  use precision, only: dp, grid_p
  ! use atomlist, only: no_u, no_l, indxuo, iaorb, qtot
  use siesta_geom, only: na_u, isa
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
    !! NOTE: Check `tab_local` in QE.
    !! Seems that it should be initialized only once.
    !! Thus, this subroutine -> `compute_tab_local`.
    use m_virtual_step_data
    use zero_flux_data, only: nqxq, pref, dq !, tpiba
    use zero_flux_data, only: tab_local, H_g, u_g
    use zero_flux_data, only: zero_flux_Jzero
    use ks_flux_data,   only: ks_flux_D

    use m_virtual_step_data, only: va_before_move
    use m_virtual_step_data, only: xa_before_move

    use atomlist,    only: qa
    use atm_types,   only: species, species_info
    use basis_types, only: nsp
    use radial,      only: rad_get
    use m_bessph,    only: BESSPH

    use iso_c_binding, only: c_loc, c_f_pointer

    complex(grid_p), dimension (:), pointer :: tc_v => null()
    integer  :: np

    real(dp) :: px, ux, vx, wx, xg
    integer  :: i0, i1, i2, i3
    integer  :: isp, ig, igp, igl, err
    integer  :: a, b, i, ii, it, ia, n_x, n_y, n_z
    type(species_info), pointer :: spp
    real(dp) :: rm, delta_rm, aux_fr, aux_dfdr, vl_int, qi
    real(dp) :: vel(3)          !! temp buffer for (`va_before_move` x 2)
    integer  :: in, is, iq, npts
    real(dp), allocatable :: rvals(:), aux(:)
    real(dp) :: H_g_rad(ngl_vmd,0:1)  !NOTE: tab_local maybe should also moved here
    ! as to an init routine

    !FIXME:
    ! rm = 0.0_dp
    rm = species(1)%reduced_vlocal%cutoff

    ! find rm and npts for vlocal
    do is = 1, nsp                !
       spp => species(is)
       if (spp%reduced_vlocal%cutoff < rm) then !FIXME: check for maximum number of points? UPD: or minimum?
          rm = spp%reduced_vlocal%cutoff
          npts = spp%reduced_vlocal%n
          delta_rm = spp%reduced_vlocal%delta
       end if
    end do

    allocate(rvals(npts))
    allocate(aux(npts))
    allocate(tab_local(nqxq, nsp, 0:1))

    rvals(:) = 0.0_dp
    aux(:) = 0.0_dp

    do in=1,npts                   ! init R values for selected* radial function
       rvals(in) = delta_rm * (in-1)
    enddo

    !NOTE: in QE:
    ! aux (ir) = (rgrid(nt)%r(ir)*upf(nt)%vloc(ir)+2.d0*zv(nt))* besr (ir) * rgrid(nt)%r(ir)
    do is = 1,nsp
       spp => species(is)
       do iq=1,nqxq
          qi = (iq - 1) * dq
          do in=1,npts
             call rad_get(spp%reduced_vlocal, rvals(in), aux_fr, aux_dfdr)
             aux(in) = aux_fr * bessph(0, rvals(in)*qi) * rvals(in) ! bessel_spherical for l=0
          end do
          call simpson(npts, delta_rm, aux, vl_int)
          tab_local(iq, is, 0) = vl_int * pref
          do in=1,npts
             call rad_get(spp%reduced_vlocal, rvals(in), aux_fr, aux_dfdr)
             aux(in) = aux_fr * bessph(1, rvals(in)*qi) * rvals(in) * rvals(in) ! bessel_spherical for l=1
          end do
          call simpson(npts, delta_rm, aux, vl_int)
          tab_local(iq, is, 1) = vl_int * pref
       end do
    end do

    allocate(H_g(ngm_plus_vmd, 3, 3, nsp))
    H_g = 0.0_dp

    !init of radial part
    init_h_g: do it = 1, nsp
       do igl = 1, ngl_vmd
          xg=sqrt(gl_vmd(igl))
          px = xg / dq - int(xg / dq)
          ux = 1.d0 - px
          vx = 2.d0 - px
          wx = 3.d0 - px
          i0 = int(xg / dq) + 1
          i1 = i0 + 1
          i2 = i0 + 2
          i3 = i0 + 3
          ! H_g_rad - radial integral that interpolates tab_local,
          ! tab_local is 1D, independent from the cell.
          H_g_rad(igl,0) = &
               &tab_local(i0, it, 0) * ux * vx * wx / 6.d0 + &
               &tab_local(i1, it, 0) * px * vx * wx / 2.d0 - &
               &tab_local(i2, it, 0) * px * ux * wx / 2.d0 + &
               &tab_local(i3, it, 0) * px * ux * vx / 6.d0
          H_g_rad(igl,1) = &
               &tab_local (i0, it, 1) * ux * vx * wx / 6.d0 + &
               &tab_local (i1, it, 1) * px * vx * wx / 2.d0 - &
               &tab_local (i2, it, 1) * px * ux * wx / 2.d0 + &
               &tab_local (i3, it, 1) * px * ux * vx / 6.d0

       !debug
       ! write(880,*) H_g_rad(igl, 0)
       ! write(881,*) H_g_rad(igl, 1)
       !debug
       end do


       do a=1,3
          do b=1,3
             if (a>=b) then
                do ig = gstart_vmd, ngm_plus_vmd
                   igp = igplus_vmd(ig)     ! get the global index for g-vectors
                   H_g(ig,a,b,it) =  g_vmd(a,igp) * g_vmd(b,igp) / gg_vmd(igp) * &
                        & (sqrt(gg_vmd(igp)) * H_g_rad(igtongl_vmd(ig),1) - &
                        &  2.d0 * qa(it) * pref * 2.d0 / gg_vmd(igp))
                end do
             end if
             if (a==b) then
                do ig = gstart_vmd, ngm_plus_vmd
                   igp = igplus_vmd(ig)     ! get the global index for g-vectors
                   H_g(ig,a,b,it) =  H_g(ig,a,b,it) - &
                        & (H_g_rad(igtongl_vmd(ig),0) - pref*2.d0*qa(it)/gg_vmd(igp))
                end do
                if (gstart_vmd==2) then
                   H_g(1,a,b,it)= - H_g_rad(1,0)
                end if
             end if
          end do
       end do
    end do init_h_g

    !NOTE: A.M.: "questo è necessario?"
    !      V.D.: Looks like yes, since h_a_b should be symmetric
    do a=1,3
       do b=1,3
          if (a>b) then
             do it=1,nsp
                H_g(:,b,a,it)=H_g(:,a,b,it)
             end do
             !NOTE: this seems to be not necessary indeed
             !      (Siesta's ionic flux coinsides with QE already)
             ! I_uno_g(:,b,a)=I_uno_g(:,a,b)
          end if
       end do
    end do

    !debug
    ! write(999,*) "H_g ", "a ", "b ", "ig ", "igtongl "
    ! do a=1,3
    !    do b=1,3
    !       do ig = 1, ngm_plus_vmd
    !          do it=1,nsp
    !             igp = igplus_vmd(ig)     ! get the global index for g-vectors
    !             write(999,*) H_g(ig,a,b,it), a, b, ig, igtongl_vmd(ig)
    !          end do
    !       end do
    !    end do
    ! end do
    !debug

    if (.not.(allocated(u_g))) allocate(u_g(ngm_plus_vmd,3))
    u_g(:,:) = 0.0_dp

    do a = 1,3
       do b = 1,3
          do ig = 1, ngm_plus_vmd
             igp = igplus_vmd(ig)     ! get the global index for g-vectors

             do ia = 1,na_u
                vel(:) = va_before_move(:,ia) * 2.0_dp ! velocity buffer

                !NOTE: The forward FFT in SIESTA is opposite to the one in QE.
                !      That leads to `Rho` in G-space being in opposite phase w/r to QE.
                !      Here I account for it by inverting the direction of FT for `u_g`.
                u_g(ig,a) = u_g(ig,a) - vel(b) * H_g(ig,a,b,isa(ia)) * &
                    & exp(+(0.d0,1.d0) * &  ! <- The plus sign in the exp to correspond with forward FFT in SIESTA
                    &     dot_product(g_vmd(1:3,igp), xa_before_move(1:3,ia)))
             end do
          end do
       end do
    end do

    zero_flux_Jzero(:) = 0.0_dp
    np = product(ntml_vmd)
    call c_f_pointer(c_loc(charge_g), tc_v, [np])

    jzero_local: do a = 1,3
       do ig = gstart_vmd, ngm_plus_vmd
          igp = igplus_vmd(ig)     ! get the global index
          zero_flux_Jzero(a)= zero_flux_Jzero(a) + &
               & 2.d0 * dble(tc_v(igp) * conjg(u_g(ig,a))) ! Mind the factor of 2
       end do                                              ! Summation over G>, Aris M. notes (3.8)
       if (gstart_vmd == 2) then
          zero_flux_Jzero(a) = zero_flux_Jzero(a) + dble(tc_v(1) * conjg(u_g(1,a)))
       end if
    end do jzero_local

    ! cleanup
    nullify(tc_v)
    deallocate(rvals)
    deallocate(aux)
    deallocate(u_g)
    deallocate(H_g)
    deallocate(tab_local)

    !TODO: next - NL-part of the zero current
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
