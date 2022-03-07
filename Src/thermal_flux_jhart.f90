module thermal_flux_jhart
  use precision, only: dp, grid_p
  use cellsubs,  only : volcel  ! Finds unit cell volume

  use thermal_flux_data

  implicit none

contains

  subroutine store_vhart_substep()
    !! Invoked on the each substep of the Thermal Flux run flow.
    !! Stores Hartree potential value from the `charge_g'
    !! for calculation of its derivative for Jhart flux component.

    real(dp) :: PI, PI8, volume
    real(dp), allocatable :: fact(:) !! scale factor for Hartree potential
    integer  :: ig, np

    PI  = 4.0_dp * atan(1.0_dp)
    PI8 = PI * 8._dp
    volume = volcel(cell_vmd)
    np = product(ntml_vmd)

    allocate(fact(np)); fact(:) = 0.0_dp

    ! Assign non-zero factor value for g-vectors in G>.
    ! Here we can use our G_plus-to-global indexes from `igplus_vmd`.
    factors_assignment: do ig=1,ngm_plus_vmd
       fact(igplus_vmd(ig)) = PI8 / gg_vmd(igplus_vmd(ig)) / volume
    end do factors_assignment

    ! Allocate (if not yet) and store scaled Hartree potential:
    if(.not.allocated(Vhart_save)) allocate(Vhart_save(2,np,gk_setup%dpoints))

    Vhart_save(1,:,substep) = charge_g(1,:) * fact(:)
    Vhart_save(2,:,substep) = charge_g(2,:) * fact(:)
    !NOTE: Do not rewrite `charge_g` - we need it more than once.

    deallocate(fact)             ! cleanup
  end subroutine store_vhart_substep


  subroutine compute_jhart()
    !! Assuming that V_hart(G)=0 if G is in the "negative" hemisphere
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

    use iso_c_binding, only: c_loc, c_f_pointer

    complex(grid_p), dimension (:), pointer :: tc_v => null()
    complex(grid_p), dimension (:), pointer :: tc_dot => null()

    real(dp) :: PI, PI8, volume
    real(dp), allocatable :: fact(:) !! scale factor for Hartree potential
    integer  :: ig, igp, np

    PI  = 4.0_dp * atan(1.0_dp)
    PI8 = PI * 8._dp
    volume = volcel(cell_vmd)
    np = product(ntml_vmd)

    allocate(fact(np)); fact(:) = 0.0_dp

    ! Assign non-zero factor value for g-vectors in G>.
    ! Here we can use our G_plus-to-global indexes from `igplus_vmd`.
    factors_assignment: do ig=1,ngm_plus_vmd
       igp = igplus_vmd(ig)     ! get the global index
       fact(igp) = PI8 / gg_vmd(igp) / volume
    end do factors_assignment

    ! In A. Marcolongo's notes, they take the semi-sum of base+virtual
    ! for tc_v...
    call c_f_pointer(c_loc(charge_g_base), tc_v, [np])
    call c_f_pointer(c_loc(Vhart_save(:,:,2)), tc_dot, [np])

    !NOTE: The forward FFT in SIESTA is opposite to the one in QE.
    !      That leads to `Rho' in G-space being in opposite phase w/r to QE.
    !      Here I account for it by inverting the sign of each product contribution
    !      ( as `dimag(conjg(tc_dot(igp))*tc_v(igp))` changes sign
    !        with the change of phase of multipliers ).
    reduce_hartree_flux: do ig = gstart_vmd, ngm_plus_vmd
       igp = igplus_vmd(ig)     ! get the global index
       gk_results%Jhart(1:3) = gk_results%Jhart(1:3) + &         ! <- The plus (inverted) sign is due to
            & dimag(conjg(tc_dot(igp))*tc_v(igp)*fact(igp)) * &  !    the difference in phase of `charge_g'
            & g_vmd(1:3,igp) * volume / 4.d0 / PI                !    w/r to QE.
    end do reduce_hartree_flux

    nullify(tc_v)
    nullify(tc_dot)
    deallocate(fact)
  end subroutine compute_jhart

end module thermal_flux_jhart
