module thermal_flux_jion
  use precision,   only: dp, grid_p
  use siesta_geom, only: na_u
  use atomlist,    only: amass, qa
  use cellsubs,    only: volcel

  use m_erf,       only: erf_func, erfc_func
  use thermal_flux_data

  implicit none

  real(dp), parameter :: amu_Ry = 911.444243104_dp
  !! Coeff. of conversion to atomic Rydberg units
  !! of atomic mass (1822.888486209 * 0.5_dp).
  real(dp), parameter :: coeff_vel_pow2 = 0.00234040089115_dp
  !! Correction coefficient for speed in kinetic energy part (A)
  !! being in Siesta units.

  real(dp), parameter :: e2 = 2.0_dp
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

  private
  public :: init_ion_flux_data, compute_jion

contains

  real(kind=dp) function h_(x)
    real(kind=dp) :: x
    h_=erfc_func(x)/x
  end function h_


  real(kind=dp) function hp_(x)
    real(kind=dp) :: x
    hp_=-(2.d0/sqrt(pi))*(1.d0/x)*exp(-x*x)-1.d0/(x*x)*erfc_func(x)
  end function hp_


  real(kind=dp) function modulus(vector)
    real(kind=dp) ::vector(3)
    modulus=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
  end function modulus


  real(dp) function alpha_0_lr(l_eta, G_sq, flag)
    real(dp)  :: l_eta, G_sq
    integer   :: flag
    if (flag==1) then
       alpha_0_lr = exp(-G_sq/(4.d0*l_eta))/G_sq
    else
       alpha_0_lr = 1.d0/G_sq
    end if
  end function alpha_0_lr


  subroutine pbc(vin,vout)
    !! Apply the minimum image distance in a cubic cell.
    real(dp), intent(in) :: vin(3)
    real(dp), intent(out) :: vout(3)

    integer :: i,n

    vout(:)=0.d0
    do i=1,3

       if (vin(i)>=0) then
          n=int(vin(i)/gk_setup%alat)
       else
          n=int(vin(i)/gk_setup%alat)-1
       end if

       vout(i)=vin(i)-dble(n)*gk_setup%alat
    end do

  end subroutine pbc


  subroutine add_local_first(val,pos,a,b)
    integer,  intent(in) :: a, b
    real(dp), intent(in) :: pos(3)
    real(dp), intent(inout) :: val
    real(dp) :: u(3), u_mod, n(3)
    integer  :: n_x, n_y, n_z
    real(dp)  :: val1

    val1 = 0.d0

    do n_x = -gk_setup%n_max_ewald, gk_setup%n_max_ewald
       do n_y = -gk_setup%n_max_ewald, gk_setup%n_max_ewald
          do n_z = -gk_setup%n_max_ewald, gk_setup%n_max_ewald
             n(1:3) = n_x * cell_vmd(1:3,1) +&
                  & n_y * cell_vmd(1:3,2) +&
                  & n_z * cell_vmd(1:3,3)
             u(1:3) = pos(1:3)-n(1:3)
             u_mod=modulus(u)
             val1=val1+gk_setup%eta_ewald*hp_(sqrt(gk_setup%eta_ewald)*u_mod)*u(a)*u(b)/u_mod
             if (a==b) then
                val1=val1+sqrt(gk_setup%eta_ewald)*h_(sqrt(gk_setup%eta_ewald)*u_mod)
             end if
          end do
       end do
    end do

    val = val + val1

  end subroutine add_local_first


  subroutine add_local_second(val,pos)
    real(dp), intent(in) :: pos(3)
    real(dp), intent(inout) :: val
    real(dp) :: modul, n(3), erf_val
    integer  :: n_x, n_y, n_z
    real(DP) :: val1

    val1 = 0.d0

    do n_x = -gk_setup%n_max_ewald, gk_setup%n_max_ewald
       do n_y = -gk_setup%n_max_ewald, gk_setup%n_max_ewald
          do n_z = -gk_setup%n_max_ewald, gk_setup%n_max_ewald
             n(1:3) = n_x * cell_vmd(1:3,1) +&
                  & n_y * cell_vmd(1:3,2) +&
                  & n_z * cell_vmd(1:3,3)
             modul = modulus(pos(1:3)-n(1:3))
             erf_val = erfc_func(sqrt(gk_setup%eta_ewald)*modul)
             val1 = val1 + erf_val/modul
          end do
       end do
    end do

    val = val + val1

  end subroutine add_local_second


  subroutine init_I_prime ()
    integer  :: ig, igp, n_x, n_y, n_z
    real(dp) :: n(3), modul, erf_value, omega

    omega = volcel(cell_vmd)
    I_prime = 0.d0

    do n_x= -gk_setup%n_max_ewald, gk_setup%n_max_ewald
       do n_y= -gk_setup%n_max_ewald, gk_setup%n_max_ewald
          do n_z= -gk_setup%n_max_ewald, gk_setup%n_max_ewald
             if ((n_x.ne.0).or.(n_y.ne.0).or.(n_z.ne.0)) then
                n(1:3) = n_x * cell_vmd(1:3,1) +&
                     & n_y * cell_vmd(1:3,2) +&
                     & n_z * cell_vmd(1:3,3)
                modul=modulus(n(:))
                erf_value=erfc_func(sqrt(gk_setup%eta_ewald)*modul)
                I_prime = I_prime + erf_value / modul
             end if
          end do
       end do
    end do

    I_prime=I_prime-2.d0*sqrt(gk_setup%eta_ewald/pi)

    I_prime_rec=0.d0
    do ig=gstart_vmd, ngm_plus_vmd ! start from vector whose module is >0
       igp = igplus_vmd(ig)     ! get the global index
       I_prime_rec=I_prime_rec+2.d0*(4.d0 * pi)/omega *&
            & alpha_0_lr(gk_setup%eta_ewald,gg_vmd(igp),1) ! QE: scaled, gg(igm)*tpiba2
    end do

    I_prime_rec=I_prime_rec-pi/(gk_setup%eta_ewald*omega)
    I_prime=I_prime+I_prime_rec

  end subroutine init_I_prime


  subroutine init_I_first_g ()
    integer  :: ig, igp, a, b
    real(dp) :: omega

    omega = volcel(cell_vmd)
    allocate(I_first_g(ngm_plus_vmd,3,3))

    do a=1,3
       do b=1,3
          if ( a >= b ) then
             do ig=gstart_vmd, ngm_plus_vmd
                igp = igplus_vmd(ig)     ! get the global index for g-vectors
                I_first_g(ig,a,b)=(4.d0*pi)/omega *&
                     & exp(-gg_vmd(igp)/(4.d0*gk_setup%eta_ewald))/(gg_vmd(igp)) *&
                     & (2 + gg_vmd(igp)/(2.d0*gk_setup%eta_ewald)) * &
                     & g_vmd(a,igp) * g_vmd(b,igp) / gg_vmd(igp)
             end do
             if (gstart_vmd==2) I_first_g(1,a,b) = 0.d0
          end if
       end do
    end do

  end subroutine init_I_first_g


  subroutine init_I_second_g ()
    integer  :: ig, igp
    real(dp) :: omega

    omega = volcel(cell_vmd)
    allocate(I_second_g(ngm_plus_vmd))

    do ig=gstart_vmd, ngm_plus_vmd
       igp = igplus_vmd(ig)     ! get the global index for g-vectors
       I_second_g(ig)=(4.d0 * pi)/omega * alpha_0_lr(gk_setup%eta_ewald,gg_vmd(igp),1)
    end do

    if (gstart_vmd==2) I_second_g(1) = -pi/(gk_setup%eta_ewald*omega)

  end subroutine init_I_second_g


  subroutine I_first_value(y,x)
    real(DP), intent(in)  :: x(3)
    real(DP), intent(out) :: y
    integer  :: ig, igp
    real(DP) :: scalar

    y = 0.d0
    do ig=gstart_vmd, ngm_plus_vmd
       igp = igplus_vmd(ig)     ! get the global index for g-vectors
       if (gg_vmd(igp)/(4.d0*gk_setup%eta_ewald) > 20.d0) exit
       y = y+2.d0*( I_second_g(ig)*cos(dot_product(g_vmd(1:3,igp), x(1:3))) )
    end do

    if (gstart_vmd==2) y = y + I_second_g(1)
    call add_local_second(y,x)

  end subroutine I_first_value


  subroutine I_second_value(y,x,a,b)
    real(dp), intent(out) :: y
    real(dp),intent(in)   :: x(3)
    integer,intent(in)    :: a,b
    integer  :: ig, igp
    real(dp) :: comp_iso

    y = 0.d0
    do ig=gstart_vmd, ngm_plus_vmd
       igp = igplus_vmd(ig)     ! get the global index for g-vectors
       if (gg_vmd(igp)/(4.d0*gk_setup%eta_ewald) > 20.d0) exit
       y = y + 2.d0*(I_first_g(ig,a,b)*cos(dot_product(g_vmd(1:3,igp), x(1:3))))
    end do

    if (gstart_vmd==2) y = y + I_first_g(1,a,b)
    call add_local_first(y,x,a,b)

    ! if (a==b) then
    !    call I_dos_value(comp_iso,x)
    !    y = y - comp_iso
    ! end if

  end subroutine I_second_value


  subroutine init_ion_flux_data ()
    call init_I_prime()
    call init_I_first_g()
    call init_I_second_g()
  end subroutine init_ion_flux_data


  subroutine compute_jion_a()
    real(dp) :: vel(3)
    !! Temp. storage for a velocity of an ion
    !! in Rydberg units (2x Hartree).
    integer :: ia

    do ia=1,na_u
       ! vel(:) = va_before_move(:,ia) * 2.0_dp
       ! vel(:) = va_before_move(:,ia)
       vel(:) = va_in(:,ia)

       gk_results%jion_a(:) = gk_results%jion_a(:) +&
            &  vel(:) * (0.5_dp * amu_Ry * amass(ia)) * (vel(1)**2 + vel(2)**2 + vel(3)**2) &
            &  * coeff_vel_pow2
    end do
  end subroutine compute_jion_a


  subroutine compute_jion_b()
    real(dp) :: vel(3)
    integer  :: ia

    do ia=1,na_u
       ! vel(:) = va_before_move(:,ia) * 2.0_dp
       vel(:) = va_in(:,ia)

       gk_results%jion_b(:) = gk_results%jion_b(:) + &
            & 2./3. * e2 * qa(ia)**2 * vel(:) * I_prime
    end do

  end subroutine compute_jion_b


  subroutine compute_jion_cde()
    real(dp) :: u(3), u_pbc(3), buff
    real(dp) :: vel_i(3), vel_j(3)
    integer  :: ia, ja, a, b

    do ia=1,na_u
       do ja=1,na_u
          if ( ia > ja ) then
             ! vel_i(:) = va_before_move(:,ia) * 2.0_dp
             ! vel_j(:) = va_before_move(:,ja) * 2.0_dp
             vel_i(:) = va_in(:,ia)
             vel_j(:) = va_in(:,ja)
             u(:) = xa_in(:,ia) - xa_in(:,ja)

             call pbc(u(1:3), u_pbc(1:3))
             call I_first_value(buff, u_pbc)

             gk_results%jion_c(:) = gk_results%jion_c(:) + &
                  & 1.0_dp * e2 * qa(ia) * qa(ja) * &
                  & (vel_i(:) + vel_j(:)) * buff

             do a=1,3
                do b=1,3
                   if (a > b) then
                      call I_second_value(buff, u_pbc, a, b)

                      gk_results%jion_e(a) = gk_results%jion_e(a) - &
                           & 1./2. * e2 * qa(ia) * qa(ja) * &
                           & (vel_i(b) + vel_j(b)) * buff

                      gk_results%jion_e(b) = gk_results%jion_e(b) - &
                           & 1./2. * e2 * qa(ia) * qa(ja) * &
                           & (vel_i(a) + vel_j(a)) * buff
                   end if

                   if (a == b) then
                      call I_second_value(buff, u_pbc, a, b)

                      gk_results%jion_d(a) = gk_results%jion_d(a) - &
                           & 1./2. * e2 * qa(ia) * qa(ja) * &
                           & (vel_i(b) + vel_j(b)) * buff
                   end if
                end do
             end do
          end if
       end do
    end do

  end subroutine compute_jion_cde


  subroutine compute_jion
    call compute_jion_a()
    gk_results%jion(:) = gk_results%jion(:) + gk_results%jion_a(:)

    !>NOTE:
    !> The following contribution behaves like a non-diffusive mass flux.
    !> According to Davide Tisi, It does not contribute to the resulting
    !> thermal conductivity and is therefore excluded from final flux
    !> calculation by default.
    !>TODO:
    !> Add optional contribution of B-flux to the final J.
    call compute_jion_b()
    gk_results%jion(:) = gk_results%jion(:) + gk_results%jion_b(:)

    call compute_jion_cde()
    gk_results%jion(:) = gk_results%jion(:) + gk_results%jion_c(:)
    gk_results%jion(:) = gk_results%jion(:) + gk_results%jion_d(:)
    gk_results%jion(:) = gk_results%jion(:) + gk_results%jion_e(:)
  end subroutine compute_jion

end module thermal_flux_jion
