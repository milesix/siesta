module hartree_flux_procs
  use precision, only: dp, grid_p
  use alloc,     only : re_alloc, de_alloc
  use cellsubs,  only : volcel  ! Finds unit cell volume

  use hartree_flux_data
  use m_virtual_step_data

  implicit none

contains

  subroutine init_Jhart_BASE()
    !! Invoked on the `BASE` step.
    !! Stores Hartree potential value from the `charge_g`
    !! for calculation of Jhart.
    real(dp) :: PI, PI8, volume
    real(dp), allocatable :: fac(:) !! scale factor for Hartree potential
    integer  :: ig, np

    PI  = 4.0_dp * atan(1.0_dp)
    PI8 = PI * 8._dp
    volume = volcel(cell_vmd)
    np = product(ntml_vmd)

    allocate(fac(np)); fac(:) = 0.0_dp

    ! Assign non-zero factor value for g-vectors in G>.
    ! Here we can use our G_plus-to-global indexes from `igplus_vmd`.
    factors_assignment: do ig=1,ngm_plus_vmd
       fac(igplus_vmd(ig)) = PI8 / gg_vmd(igplus_vmd(ig)) / volume
    end do factors_assignment

    ! Allocate and store scaled Hartree potential:
    call re_alloc(Vhart_deriv,1,2,1,np,'Vhart_deriv','init_Jhart')
    Vhart_deriv(1,:) = charge_g(1,:) * fac(:)
    Vhart_deriv(2,:) = charge_g(2,:) * fac(:)
    !NOTE: Do not rewrite `charge_g` - we need it more than once.

    deallocate(fac)             ! cleanup
  end subroutine init_Jhart_BASE


  subroutine compute_Jhart_VIRTUAL()
    !! Invoked on the `VIRTUAL` step.
    !! This means that V_hart(G)=0 if G is in the "negative" hemisphere
    !! I assume that this is done to do further operations in the same way
    !! as QE. In particular using a half-space avoids having to find -G
    !! as G's are generated.
    !! Just DO NOT use this as the "real" V_H(G) later on
    !! Operations logic here:
    !! 1. get the derivative of Hartree potential;
    !! 2. grab C-to-Fortran pointers to the potential and its derivative
    !! (as complex numbers);
    !! 3. second loop over indeces of the "positive" hemisphere
    !! to compute the Hartree flux component
    !! according to the formula 1.17 in Aris Marcolongo's manual.
    !!@note
    !! The sign of the terms for the sum in [1.17] is opposite to ours:
    !! that is due to opposite phases of FFT used in SIESTA and QE.
    !!@endnote
    use siesta_options, only: virtual_dt
    use iso_c_binding, only: c_loc, c_f_pointer

    complex(grid_p), dimension (:), pointer :: tc_v => null()
    complex(grid_p), dimension (:), pointer :: tc_dot => null()

    real(dp) :: PI, PI8, volume
    real(dp), allocatable :: fac(:) !! scale factor for Hartree potential
    integer  :: ig, igp, np

    PI  = 4.0_dp * atan(1.0_dp)
    PI8 = PI * 8._dp
    volume = volcel(cell_vmd)
    np = product(ntml_vmd)

    allocate(fac(np)); fac(:) = 0.0_dp

    ! Assign non-zero factor value for g-vectors in G>.
    ! Here we can use our G_plus-to-global indexes from `igplus_vmd`.
    factors_assignment: do ig=1,ngm_plus_vmd
       igp = igplus_vmd(ig)     ! get the global index
       fac(igp) = PI8 / gg_vmd(igp) / volume
    end do factors_assignment

    ! Obtain V_Hartree derivative:
    Vhart_deriv(1,:) = (charge_g(1,:)*fac(:) - Vhart_deriv(1,:))/virtual_dt
    Vhart_deriv(2,:) = (charge_g(2,:)*fac(:) - Vhart_deriv(2,:))/virtual_dt

    h_flux_Jhart(:) = 0.d0      ! init flux accumulator

    ! In A. Marcolongo's notes, they take the semi-sum of base+virtual
    ! for tc_v...
    call c_f_pointer(c_loc(charge_g), tc_v, [np])
    call c_f_pointer(c_loc(Vhart_deriv), tc_dot, [np])

    reduce_hartree_flux: do ig = gstart_vmd, ngm_plus_vmd
       igp = igplus_vmd(ig)     ! get the global index
       h_flux_Jhart(1:3) = h_flux_Jhart(1:3) + &
            & dimag(conjg(-tc_dot(igp))*tc_v(igp)*fac(igp)) * &
            & g_vmd(1:3,igp) * volume / 4.d0 / PI
    end do reduce_hartree_flux

    nullify(tc_v)               ! cleanup
    nullify(tc_dot)
    deallocate(fac)
  end subroutine compute_Jhart_VIRTUAL

end module hartree_flux_procs
