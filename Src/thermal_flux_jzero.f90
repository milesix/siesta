module thermal_flux_jzero
  !! Contains procedures for computation of J_zero.
  !! See notes on components of the Kohn-Sham flux for thermal transport.
  use precision, only: dp, grid_p
  ! use atomlist, only: no_u, no_l, indxuo, iaorb, qtot
  use siesta_geom, only: na_u, isa
  use atomlist, only: no_u, no_l, iaorb, iphorb
  use sparse_matrices, only: listh, listhptr, numh

  use matel_mod, only: new_matel
  use thermal_flux_data

  implicit none

  real(DP), allocatable, save :: tab_local(:,:,:)
  !! One-dimensional table for radial pseudopotential part.
  real(dp), allocatable, save :: H_g(:,:,:,:)
  !! Table for the values of `h_ab` used in computation
  !! of the local, long-range part of the Zero heat current.
  complex(DP), allocatable, save :: u_g(:,:)

  real(dp), parameter :: dq = 0.01_dp
  !! space between points in the pseudopotential tab.
  !!FIXME: constant?
  real(dp), parameter :: cell_factor = 1.0_dp
  !! maximum expected (linear) cell contraction
  !! during relaxation/MD
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: e2 = 2.0_dp !! electron charge squared

  character(len=1), dimension(3), parameter :: coord_table = [ 'X', 'Y', 'Z' ]

  real(dp)  :: omega
  real(dp)  :: pref   !! `4pi/omega` prefactor
  integer   :: nqxq   !! size of interpolation table
  ! real(dp)  :: tpiba

  private
  public :: init_zero_flux_data, compute_jzero

contains

  subroutine init_zero_flux_data ()
    use fdf, only: fdf_physical
    use utils, only: die
    use cellsubs, only : volcel
    !NOTE: This should be present even with `ion_flux` calculation switched off,
    ! though ofc it should be moved into a proper scope.

    ! real(dp)  :: gcutm  !! potential cut-off

    omega = volcel(cell_vmd)
    pref  = 4.0_dp * pi / omega !! `4pi/omega` prefactor
    ! tpiba = 2.0_dp * pi / gk_setup%alat

    ! gcutm = meshcutoff / tpiba**2
    nqxq  = int((sqrt(gk_setup%meshcutoff)/dq + 4) * cell_factor)

  end subroutine init_zero_flux_data


  subroutine compute_jzero ()
    !! NOTE: Check `tab_local` in QE.
    !! Seems that it should be initialized only once.
    !! Thus, this subroutine -> `compute_tab_local`.
    ! use atomlist,    only: qa
    use atm_types,   only: species, species_info
    use basis_types, only: nsp
    use radial,      only: rad_get
    use m_bessph,    only: BESSPH

    use iso_c_binding, only: c_loc, c_f_pointer

    use siesta_geom
    use atomlist,    only: no_u, no_s, no_l, indxuo, iaorb, iphorb, lasto, qtot, amass, lastkb, iphKB
    use neighbour,   only: jna=>jan, xij, r2ij
    use neighbour,   only: mneighb, reset_neighbour_arrays
    use atmfuncs,    only: rcut, orb_gindex, kbproj_gindex
    use matel_mod,   only: is_orb, is_kb

    ! integer &  ia, iio, ind, io, ioa, is, ix, &  j, ja, jn, jo, joa, js, nnia, norb, ka
    ! real(dp) &  grSij(3), rij, Sij, xinv(3), sum
    ! integer &  ikb, ina, ino, jno, ko, koa, ks, ig, jg, kg, &  nkb, nna, nno, ilm1, ilm2, npoints, &  ir
    ! real(dp) ::  epsk, grSki(3), rki, rmax, rmaxkb, rmaxo
    real(dp) :: rmax, rmaxkb, rmaxo, qa_sp
    integer  :: iio, ind, io, ioa, j, ja, jn, jo, joa, js, jg, nnia, nna, norb, ks, koa, ikb
    integer  :: na, no, nua, nuo, nai, naj
    real(dp) :: Jznl_alt(3)

    complex(grid_p), dimension (:), pointer :: tc_v => null()
    integer  :: np

    real(dp) :: px, ux, vx, wx, xg
    integer  :: i0, i1, i2, i3
    integer  :: isp, ig, igp, igl, err
    integer  :: a, b, ii, it, ia, n_x, n_y, n_z
    integer  :: I_ind, alpha, r_ind, grad_ind
    integer  :: mu, nu, iamu, ianu, inda, iph_nu, iph_mu
    real(dp) :: val, grad(3), J_tmp(3), tmp_nl
    real(dp) :: R_I(3), R_mu(3), R_nu(3), R_Imu(3), R_Inu(3)
    type(species_info), pointer :: spp
    real(dp) :: rm, delta_rm, aux_fr, aux_dfdr, vl_int, qi
    real(dp) :: vel(3)          !! temp buffer for (`va_in`)
    integer  :: in, is, iq, npts
    real(dp), allocatable :: aux(:)
    real(dp) :: H_g_rad(ngl_vmd,0:1)  !NOTE: tab_local maybe should also moved here
    real(dp) :: rbuf            !! temp buffer for R distance values
    ! as to an init routine

    allocate(tab_local(nqxq, nsp, 0:1))
    tab_local(:,:,:) = 0.0_dp

    !NOTE: in QE:
    ! aux (ir) = (rgrid(nt)%r(ir)*upf(nt)%vloc(ir)+2.d0*zv(nt))* besr (ir) * rgrid(nt)%r(ir)
    do is = 1,nsp
       spp => species(is)

       rm = spp%reduced_vlocal%cutoff
       npts = spp%reduced_vlocal%n
       npts =  npts + mod(npts,2)  ! for Simpson needs 2n points
       delta_rm = rm / npts

       allocate(aux(npts))
       aux(:) = 0.0_dp

       do iq=1,nqxq
          qi = (iq - 1) * dq

          do in=1,npts
             rbuf = (delta_rm * (in-1))
             call rad_get(spp%reduced_vlocal, rbuf, aux_fr, aux_dfdr)
             aux(in) = aux_fr * bessph(0, rbuf*qi) * rbuf ! bessel_spherical for l=0
          end do
          call simpson(npts, delta_rm, aux, vl_int)
          tab_local(iq, is, 0) = vl_int * pref

          do in=1,npts
             rbuf = (delta_rm * (in-1))
             call rad_get(spp%reduced_vlocal, rbuf, aux_fr, aux_dfdr)
             aux(in) = aux_fr * bessph(1, rbuf*qi) * rbuf * rbuf ! bessel_spherical for l=1
          end do
          call simpson(npts, delta_rm, aux, vl_int)
          tab_local(iq, is, 1) = vl_int * pref

       end do

       deallocate(aux)

    end do

    allocate(H_g(ngm_plus_vmd, 3, 3, nsp))
    H_g = 0.0_dp

    !init of radial part
    init_h_g: do it = 1, nsp
       do igl = 1, ngl_vmd
          xg = sqrt(gl_vmd(igl))!* tpiba  !(not to) cut below this radial interpolation table
          !NOTE: maybe there is a more accurate way to build it.
          px = xg / dq - int(xg / dq)
          ux = 1.d0 - px
          vx = 2.d0 - px
          wx = 3.d0 - px
          i0 = int(xg / dq) + 1
          i1 = i0 + 1
          i2 = i0 + 2
          i3 = i0 + 3

          if ( i3 <= nqxq ) then
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
          else
             H_g_rad(igl,0:1) = 0.0_dp
          end if

          !debug
          ! write(880,*) H_g_rad(igl, 0)
          ! write(881,*) H_g_rad(igl, 1)
          !debug
       end do

       spp => species(it)
       qa_sp = spp%zval

       do a=1,3
          do b=1,3
             if (a>=b) then
                do ig = gstart_vmd, ngm_plus_vmd
                   igp = igplus_vmd(ig)     ! get the global index for g-vectors
                   H_g(ig,a,b,it) =  g_vmd(a,igp) * g_vmd(b,igp) / gg_vmd(igp) * &
                        & (sqrt(gg_vmd(igp)) * H_g_rad(igtongl_vmd(ig),1) - &
                        &  2.0_dp * e2 * qa_sp * pref / gg_vmd(igp))
                end do
             end if
             if (a==b) then
                do ig = gstart_vmd, ngm_plus_vmd
                   igp = igplus_vmd(ig)     ! get the global index for g-vectors
                   H_g(ig,a,b,it) =  H_g(ig,a,b,it) - &
                        & (H_g_rad(igtongl_vmd(ig),0) - e2*pref*qa_sp/gg_vmd(igp))
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
                ! vel(:) = va_in(:,ia) * 2.0_dp ! velocity buffer
                vel(:) = va_in(:,ia) ! velocity buffer

                !NOTE: The forward FFT in SIESTA is opposite to the one in QE.
                !      That leads to `Rho` in G-space being in opposite phase w/r to QE.
                !      Here I account for it by inverting the direction of FT for `u_g`.
                u_g(ig,a) = u_g(ig,a) - vel(b) * H_g(ig,a,b,isa(ia)) * &
                    & exp(+(0.d0,1.d0) * &  ! <- The plus sign in the exp to correspond with forward FFT in SIESTA
                    &     dot_product(g_vmd(1:3,igp), xa_in(1:3,ia)))
             end do
          end do
       end do
    end do

    np = product(ntml_vmd)
    call c_f_pointer(c_loc(charge_g_base), tc_v, [np])

    ! do ig = gstart_vmd, ngm_plus_vmd
    !    igp = igplus_vmd(ig)     ! get the global index

    !    print*, "[DBG]", ig, igp
    !    print*, "[DBG]", tc_v(igp)
    !    print*, "[DBG]", u_g(ig,:)
    !    print*, "[DBG]", "----------------------------"

    ! end do

    jzero_local: do a = 1,3
       do ig = gstart_vmd, ngm_plus_vmd
          igp = igplus_vmd(ig)     ! get the global index
          gk_results%Jzloc(a)= gk_results%Jzloc(a) + &
               & 2.d0 * dble(tc_v(igp) * conjg(u_g(ig,a))) ! Mind the factor of 2

          write(700+a, *) gk_results%Jzloc(a) !debug
       end do                                              ! Summation over G>, Aris M. notes (3.8)
       if (gstart_vmd == 2) then
          gk_results%Jzloc(a) = gk_results%Jzloc(a) + dble(tc_v(1) * conjg(u_g(1,a)))
       end if
    end do jzero_local

    ! cleanup
    nullify(tc_v)
    deallocate(u_g)
    deallocate(H_g)
    deallocate(tab_local)

    ! This should add NL-part of the zero current
    !FIXME: For an auxiliary supercell this will be erroneous.
    !`na_s != na_u` means aux sc is used.
    jzero_nl: do I_ind=1,na_u
       R_I(1:3) = xa_in(1:3,I_ind)
       vel(1:3) = va_in(1:3,I_ind)

       do mu=1,no_l
          iamu = iaorb(mu)
          is   = isa(iamu)
          ioa  = iphorb(mu)
          iph_mu = orb_gindex(is, ioa)

          R_mu(1:3) = xa_in(1:3,iamu)

          do nu=1,no_l
             ianu = iaorb(nu)
             js   = isa(ianu)
             joa  = iphorb(nu)
             iph_nu = orb_gindex(js, joa)
             ! ianu = iaorb(nu)
             ! iph_nu = iphorb(nu)

             R_nu(1:3) = xa_in(1:3,ianu)

             J_tmp(1:3) = 0.0_dp
             ! is = isa(I_ind)
             ! spp => species(is)

             do ikb = lastkb(I_ind-1)+1,lastkb(I_ind)
                ks = isa(I_ind)
                koa    = iphKB(ikb)
                alpha  = kbproj_gindex(ks,koa)
             ! do inda=1,spp%nprojs
             !    alpha=spp%pj_gindex(inda)

                do grad_ind = 1,3
                   call new_MATEL(coord_table(grad_ind), iph_mu, alpha, (R_I(:)-R_mu(:)), val, grad)
                   tmp_nl = val

                   call new_MATEL(SG(grad_ind), iph_nu, alpha, (R_I(:)-R_nu(:)), val, grad)

                   J_tmp(grad_ind) = J_tmp(grad_ind) - val*tmp_nl*vel(grad_ind)
                end do

                call new_MATEL('S', iph_nu, alpha, (R_I(:)-R_nu(:)), val, grad)
                tmp_nl = val

                do r_ind = 1,3
                   do grad_ind = 1,3
                      call new_MATEL(RG(r_ind,grad_ind), iph_mu, alpha, (R_I(:)-R_mu(:)), val, grad)

                      J_tmp(r_ind) = J_tmp(r_ind) - val*tmp_nl*vel(grad_ind)
                   end do
                end do

             end do

             ! gk_results%Jznl(1:3) = gk_results%Jznl(1:3) + J_tmp(1:3) * 0.5_dp*DM_save(mu,nu,1,1) ! only 1st spin component
             gk_results%Jznl(1:3) = gk_results%Jznl(1:3) + J_tmp(1:3) * 0.5_dp*DM_save(iph_mu,iph_nu,1,1) ! only 1st spin component
          enddo
       enddo
    end do jzero_nl

    gk_results%Jzero(:) = gk_results%Jzloc(:) + gk_results%Jznl(:)

!!!!!
    ! This should really add NL-part of the zero current
    ! accounting for periodically reflected \mu and \nu.

    ! scell   <- siesta_geom
    ! rmax    <- calculated
    ! rmaxo   <- calculated
    ! na      <- siesta_geom
    ! xa      <- siesta_geom | xa_in
    ! nnia    <- output
    ! nna     <- output

    ! Inits
    na =  na_s
    nua = na_u
    no =  no_s
    nuo = no_u

! Find maximum ranges
    rmaxo = 0.0d0
    rmaxkb = 0.0d0
    do ia = 1,na
       is = isa(ia)
       do ikb = lastkb(ia-1)+1,lastkb(ia)
          ioa = iphKB(ikb)
          rmaxkb = max( rmaxkb, rcut(is,ioa) )
       enddo
       do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rmaxo = max( rmaxo, rcut(is,ioa) )
       enddo
    enddo
    rmax = rmaxo + rmaxkb

    ! Initialize neighb subroutine
    call mneighb( scell, 2.0d0*rmaxo, na, xa_in, 0, 0, nnia )

    Jznl_alt(:) = 0.0_dp

    jzero_nl_alt: do I_ind=1,na_u

       ks = isa(I_ind)

       R_I(1:3) = xa_in(1:3,I_ind)
       vel(1:3) = va_in(1:3,I_ind)

       ! Find neighbour atoms
       call mneighb( scell, rmax, na, xa_in, I_ind, 1, nna ) ! < isc = 1 to include self

       do nai = 1, nna
          ia = jna(nai)
          is = isa(ia)

          do naj = 1, nna
             ja = jna(naj)
             js = isa(ja)

                R_Imu(1:3) = xij(1:3, ia)
                R_Inu(1:3) = xij(1:3, ja)

                do mu = lasto(ia-1)+1,lasto(ia)
                   ioa = iphorb(mu)
                   iph_mu = orb_gindex(is,ioa)

                   do nu = lasto(ja-1)+1,lasto(ja)
                      joa = iphorb(nu)
                      iph_nu = orb_gindex(js,joa)

                      J_tmp(1:3) = 0.0_dp

                      do ikb = lastkb(I_ind-1)+1,lastkb(I_ind)
                         koa    = iphKB(ikb)
                         alpha  = kbproj_gindex(ks,koa)

                         ! write(5000, *) iph_mu, iph_nu, alpha
                         ! write(5000, *) is_orb(iph_mu), is_orb(iph_nu), is_kb(alpha)
                         ! write(5000, *)

                         do grad_ind = 1,3
                            call new_MATEL(coord_table(grad_ind), iph_mu, alpha, R_Imu(1:3), val, grad)
                            tmp_nl = val

                            call new_MATEL(SG(grad_ind), iph_nu, alpha, R_Inu(1:3), val, grad)

                            J_tmp(grad_ind) = J_tmp(grad_ind) - val*tmp_nl*vel(grad_ind)
                         end do

                         call new_MATEL('S', iph_nu, alpha, R_Inu(1:3), val, grad)
                         tmp_nl = val

                         do r_ind = 1,3
                            do grad_ind = 1,3
                               call new_MATEL(RG(r_ind,grad_ind), iph_mu, alpha, R_Imu(1:3), val, grad)

                               J_tmp(r_ind) = J_tmp(r_ind) - val*tmp_nl*vel(grad_ind)
                            end do
                         end do
                      end do

                      ! Jznl_alt(1:3) = Jznl_alt(1:3) + J_tmp(1:3) * 0.5_dp*DM_save(mu,nu,1,1) ! only 1st spin component
                      Jznl_alt(1:3) = Jznl_alt(1:3) + J_tmp(1:3) * 0.5_dp*DM_save(iph_mu,iph_nu,1,1) ! only 1st spin component
                   end do
                end do

          end do
       end do

       ! write(5000, *) "------------------------------"

    end do jzero_nl_alt

    print*, "[test_Jznl_alt]", Jznl_alt * 0.0483776900146_dp  !to compare with QE at once
    call reset_neighbour_arrays( )
!!!!!

  end subroutine compute_jzero


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

end module thermal_flux_jzero
