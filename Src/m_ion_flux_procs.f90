module ion_flux_procs
  !! Contains procedures for computation of components of J_ion part of heat flux.
  use precision, only: dp
  use siesta_geom, only: na_u
  use atomlist, only: amass, qa

  use m_virtual_step_data, only: va_before_move
  use m_virtual_step_data, only: xa_before_move

  implicit none

  integer :: ia, ja, a, b, iu, iv, ix
  real(dp), parameter :: amu_Ry = 911.444243104_dp
  !! Coeff. of conversion to atomic Rydberg units
  !! of atomic mass (1822.888486209 * 0.5_dp).
  real(dp), parameter :: coeff_vel_pow2 = 0.00234040089115_dp
  !! Correction coefficient for speed in kinetic energy part (A)
  !! being in Siesta units.

  real(dp), parameter :: e2 = 2.0_dp

  private
  public :: compute_Jion_a, compute_Jion_b
  public :: compute_Jion_cde

contains

  subroutine compute_Jion_a()
    use ion_flux_data, only: ion_flux_a, ion_flux_Jion

    real(dp) :: vel(3)
    !! Temp. storage for a velocity of an ion
    !! in Rydberg units (2x Hartree).

    do ia=1,na_u
       ! vel(:) = va_before_move(:,ia) * 2.0_dp
       vel(:) = va_before_move(:,ia)

       ion_flux_a(:) = ion_flux_a(:) +&
            &  vel(:) * (0.5_dp * amu_Ry * amass(ia)) * (vel(1)**2 + vel(2)**2 + vel(3)**2) &
            &  * coeff_vel_pow2
    end do
  end subroutine compute_Jion_a


  subroutine compute_Jion_b()
    use ion_flux_data, only: ion_flux_b, I_prime

    real(dp) :: vel(3)

    do ia=1,na_u
       ! vel(:) = va_before_move(:,ia) * 2.0_dp
       vel(:) = va_before_move(:,ia)

       ion_flux_b(:) = ion_flux_b(:) + &
            & 2./3. * e2 * qa(ia)**2 * vel(:) * I_prime
    end do

  end subroutine compute_Jion_b


  subroutine compute_Jion_cde()
    use ion_flux_data, only: ion_flux_c, ion_flux_d
    use ion_flux_data, only: ion_flux_e, ion_flux_Jion
    use ion_flux_data, only: pbc, I_uno_value, I_dos_value

    real(dp) :: u(3), u_pbc(3), buff
    real(dp) :: vel_i(3), vel_j(3)

    do ia=1,na_u
       do ja=1,na_u
          if ( ia > ja ) then
             ! vel_i(:) = va_before_move(:,ia) * 2.0_dp
             ! vel_j(:) = va_before_move(:,ja) * 2.0_dp
             vel_i(:) = va_before_move(:,ia)
             vel_j(:) = va_before_move(:,ja)
             u(:) = xa_before_move(:,ia) - xa_before_move(:,ja)

             call pbc(u(1:3), u_pbc(1:3))
             call I_dos_value(buff, u_pbc)

             ion_flux_c(:) = ion_flux_c(:) + &
                  & 1.0_dp * e2 * qa(ia) * qa(ja) * &
                  & (vel_i(:) + vel_j(:)) * buff

             do a=1,3
                do b=1,3
                   if (a > b) then
                      call I_uno_value(buff, u_pbc, a, b)

                      ion_flux_e(a) = ion_flux_e(a) - &
                           & 1./2. * e2 * qa(ia) * qa(ja) * &
                           & (vel_i(b) + vel_j(b)) * buff

                      ion_flux_e(b) = ion_flux_e(b) - &
                           & 1./2. * e2 * qa(ia) * qa(ja) * &
                           & (vel_i(a) + vel_j(a)) * buff
                   end if

                   if (a == b) then
                      call I_uno_value(buff, u_pbc, a, b)

                      ion_flux_d(a) = ion_flux_d(a) - &
                           & 1./2. * e2 * qa(ia) * qa(ja) * &
                           & (vel_i(b) + vel_j(b)) * buff
                   end if
                end do
             end do
          end if
       end do
    end do

  end subroutine compute_Jion_cde

end module ion_flux_procs
