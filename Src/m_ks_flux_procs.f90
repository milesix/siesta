module ks_flux_procs
  !! Contains procedures for computation of J_ks.
  !! See notes on components of the Kohn-Sham flux for thermal transport.
  use precision, only: dp
  use atomlist, only: no_u, no_l, indxuo, iaorb, qtot
  use densematrix, only: psi ! To get at the c^i_\mu coefficients of the wavefunctions
  use sparse_matrices, only: listh, listhptr, numh
  use m_virtual_step_data, only: Dfull, Dderiv
  !!TODO: are the dimensions of Dfull & Dderiv correct?

  implicit none

  private
  public :: compute_Jks

contains

subroutine compute_psi_hat_c (psi_hat_c, n_wfs)
  !! We need to multiply three matrices:  coeffs, rmatrix, dscf_hat
  !!
  !! psi_hat_c = matmul(coeffs,rmatrix)
  !! psi_hat_c = matmul(psi_hat_c,dscf_hat)
  !!
  !! coeffs is distributed in parallel, with no_u rows and no_l columns.
  !!
  !! rmatrix is distributed in parallel, with no_l rows and no_u columns.
  !! dscf_hat is distributed in parallel, with no_l rows and no_u columns.
  !! rmatrix is formally sparse, but dscf_hat IS NOT, since we need the
  !! "full density matrix" and not the "sparse section" that is enough for
  !! the computation of the charge density.
  !!
  !! We will do first a serial version. This means that no_l = no_u.
  !! BUT, we might want to deal with only a subset of the wavefunctions,
  !! so each processor will keep only n_wfs of them (n_wfs <= no_l)
  !!
  !! Then:
  !!
  !! coeffs has no_u rows and n_wfs columns.
  !!
  !! rmatrix has no_l rows and no_u columns.
  !! dscf_hat has no_l rows and no_u columns.
  !! these two are also sparse.
  !!
  !! psi_hat_c will have no_u rows and n_wfs columnns

  ! To get the density and <phi_1|r|phi_2> matrix
  use ks_flux_data, only: ks_flux_D, ks_flux_Rmat

  ! Output: coefficients of psi_hat_c, stored in module ks_flux_data
  ! use ks_flux_data, only: psi_hat_c

  use,intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

  real(dp), allocatable, intent(inout) :: psi_hat_c(:,:,:)
  integer, intent(in) :: n_wfs
  !! Number of wavefunctions to process

  real(dp), dimension(:), pointer :: rmatrix => null()
  real(dp), dimension(:), pointer :: dscf_hat => null()

  !> Coordinate index regarding the selected dimension of Rmat.
  !> It will define X/Y/Z component of resulting Jks.
  integer :: idx
  integer :: nu, mu, alpha, iw, ij, k, col

  !> First dimension: coeff index,
  !> Second dimension: wavefunction index (only no_l of them)
  real(dp), pointer :: coeffs(:,:)
  real(dp), allocatable :: tmp_g(:)

  call c_f_pointer(c_loc(psi(1)), coeffs, [no_u, no_l]) ! for new compilers
  ! explicit pointing to psi(1) -> psi; avoid gcc<=4.8

  allocate (tmp_g(no_l))

  dscf_hat => ks_flux_D(:,1)    !NOTE: Just first spin component
  ! computing and storage of every of the 3 components
  ! should be available to order in the .fdf

  psi_hat_c(:,:,:) = 0.0        ! Init result to zeros outside main loop

  DO idx=1,3
     rmatrix => ks_flux_Rmat(:,idx)


  do iw = 1,n_wfs
     tmp_g(:) = 0.0
     do nu = 1,no_l
        do mu = 1,no_u
          do ij = 1, numh(nu)
             k = listhptr(nu) + ij
             col = listh(k)
             if (col == mu) then
               tmp_g(nu) = tmp_g(nu) + coeffs(mu,iw) * rmatrix(k)
               exit
             endif
          end do
        end do
     end do

     ! Use dense, complete, DM, but note that, for the spinless case, DM = 2 \sum_occ { |psi><psi| }
     ! We need then to use a factor of 1/2

     do alpha = 1,no_u
        do nu = 1,no_l

           if (alpha == nu) then
              psi_hat_c(alpha,iw,idx) = psi_hat_c(alpha,iw,idx) + &
                   tmp_g(nu) * (1.0 - 0.5_dp * Dfull(nu,alpha,1)) ! <- only 1st spin component
           else
              psi_hat_c(alpha,iw,idx) = psi_hat_c(alpha,iw,idx) + &
                   tmp_g(nu) * (- 0.5_dp * Dfull(nu,alpha,1)) ! <- only 1st spin component
           end if

        end do
     end do

  end do ! iw
  END DO ! idx

  deallocate(tmp_g)

end subroutine compute_psi_hat_c

!
!  This routine is very complex. It has to be re-checked
!
subroutine compute_psi_dot_c (psi_dot_c, n_wfs)
  use siesta_options, only : virtual_dt
  ! Take proper velocities from `BASE` step:
  use m_virtual_step_data, only: va_before_move

  ! Entities stored during the BASE MD-step:
  use ks_flux_data, only: ks_flux_D, ks_flux_S, ks_flux_gS
  ! Derivative of the density (computed earlier during VIRTUAL step)
  use ks_flux_data, only: Dscf_deriv

  ! Output: coefficients of psi_dot_c, stored in module ks_flux_data
  ! use ks_flux_data, only: psi_dot_c

  use,intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

  implicit none

  real(dp), allocatable, intent(inout) :: psi_dot_c(:,:)
  integer, intent(in) :: n_wfs
  !! Number of wavefunctions to process

  real(dp), dimension(:),   pointer :: dscf_hat => null()
  real(dp), dimension(:,:), pointer :: gradS => null()
  real(dp), dimension(:),   pointer :: S => null()

  integer  :: iw
  integer  :: alpha, beta, lambda, nu, mu
  integer  :: ind, col, sec_ind, sec_col
  real(dp) :: acc_left

  !> First dimension: coeff index,
  !> Second dimension: wavefunction index (only no_l of them)
  real(dp), pointer :: coeffs(:,:)

  real(dp), allocatable :: tmp_left(:)
  real(dp), allocatable :: tmp_right(:)

  call c_f_pointer(c_loc(psi(1)), coeffs, [no_u, no_l]) ! for new compilers
  ! explicit pointing to psi(1) -> psi
  ! avoid gcc<=4.8

  allocate (tmp_left(no_u))
  allocate (tmp_right(no_u))

  ! Finite-difference time derivative of Dfull:
  Dderiv(:,:,1) = (Dfull(:,:,1) - Dderiv(:,:,1))/virtual_dt !NOTE: Just first spin component

  dscf_hat => ks_flux_D(:,1)         !NOTE: Just first spin component
  S => ks_flux_S(:)
  gradS => ks_flux_gS(:,:)
  ! computing and storage of every of the 3 components
  ! should be available to order in the .fdf

  psi_dot_c(:,:) = 0.0          ! Init result to zeros outside main loop

  ! call timer("explicit_matmul_Jks", 1)

  ! 1st row of multiplication schematic:
  do iw = 1,n_wfs
     ! Store matrix multiplication with coeffs in the
     ! "right" temporary buffer array:
     tmp_right(:) = 0.0
     do nu=1,no_u
        do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
           alpha = listh(ind)
           tmp_right(nu) = tmp_right(nu) + &
                & S(ind) * coeffs(alpha,iw)
        end do
     end do

     ! Use "left" buffer array to compute center-right multiplication,
     ! then store the resulting column in the "right" buffer array:
     tmp_left(:) = 0.0
     do mu=1,no_u
        ! do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
        !    nu = listh(ind)
        !    tmp_left(mu) = tmp_left(mu)+Dscf_deriv(ind,1)*tmp_right(nu)
        !    !NOTE: Only first spin component here now--^
        ! end do
        do nu=1,no_l
           tmp_left(mu) = tmp_left(mu)+ 0.5_dp * Dderiv(mu,nu,1)*tmp_right(nu)
         !NOTE: Only first spin component here now--^
        end do
     end do
     tmp_right(:) = tmp_left(:)

     ! Compute the "left" matrix multiplication,
     ! then add to the resulting array:
     do lambda=1,no_u
        do mu=1,no_u
           acc_left = 0.0
           ! do ind = (listhptr(lambda)+1),&
           !      & listhptr(lambda) + numh(lambda)
           !    beta = listh(ind)
           do beta=1,no_l       ! Dfull is dense now
              do sec_ind = (listhptr(beta)+1),&
                   & listhptr(beta) + numh(beta)
                 sec_col = listh(sec_ind)
                 if ( sec_col .eq. mu ) then
                    if (beta.eq.lambda) then
                       acc_left = acc_left + &
                            &(1.0_dp - 0.5_dp * Dfull(lambda,beta,1)) * S(sec_ind)
                            ! &(1.0 - dscf_hat(ind)) * S(sec_ind)
                    else
                       acc_left = acc_left + &
                            &(0.0_dp - 0.5_dp * Dfull(lambda,beta,1)) * S(sec_ind)
                            ! &(0.0 - dscf_hat(ind)) * S(sec_ind)
                    end if
                    exit
                 end if
              end do
           end do

           ! add dot-product element to the resulting array:
           psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw) +&
                & acc_left * tmp_right(mu)
        end do
     end do
  end do                        ! <-- end of 1st matmul row

  ! 2nd row of multiplication schematic:
  do iw = 1,n_wfs
     ! Store matrix multiplication with coeffs in the
     ! "right" temporary buffer array:
     tmp_right(:) = 0.0
     do nu=1,no_u
        do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
           alpha = listh(ind)
           tmp_right(nu) = tmp_right(nu) + &
                & S(ind) * coeffs(alpha,iw)
        end do
     end do

     ! Use "left" buffer array to compute center-right multiplication,
     ! then store the resulting column in the "right" buffer array:
     tmp_left(:) = 0.0
     do mu=1,no_u
        ! do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
        !    nu = listh(ind)
        !    tmp_left(mu) = tmp_left(mu) + dscf_hat(ind) * tmp_right(nu)
        ! end do
        do nu=1,no_l
           tmp_left(mu) = tmp_left(mu) + 0.5_dp * Dfull(mu,nu,1) * tmp_right(nu)
        end do
     end do
     tmp_right(:) = tmp_left(:)

     ! Compute the "left" matrix multiplication,
     ! then add to the resulting array:
     do lambda=1,no_u
        do mu=1,no_u
           acc_left = 0.0
           do beta=1,no_l     ! Dfull is dense now
           ! do ind = (listhptr(lambda)+1),&
           !      & listhptr(lambda) + numh(lambda)
           !    beta = listh(ind)
              do sec_ind = (listhptr(beta)+1),&
                   & listhptr(beta) + numh(beta)
                 sec_col = listh(sec_ind)
                 if ( sec_col .eq. mu ) then
                    if (beta.eq.lambda) then
                       acc_left = acc_left + &
                            ! &(1.0 - dscf_hat(ind)) * &
                            &(1.0_dp - 0.5_dp * Dfull(lambda,beta,1)) * &
                            &(sum(gradS(sec_ind,:) * &
                            & va_before_move(:,iaorb(mu))))
                    else
                       acc_left = acc_left + &
                            ! &(0.0 - dscf_hat(ind)) * &
                            &(0.0 - 0.5_dp * Dfull(lambda,beta,1)) * &
                            &(sum(gradS(sec_ind,:) * &
                            & va_before_move(:,iaorb(mu))))
                    end if
                    exit
                 end if
              end do
           end do

           ! add dot-product element to the resulting array:
           psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw) +&
                & acc_left * tmp_right(mu)
        end do
     end do
  end do                        ! <-- end of 2nd matmul row

  ! 3rd row of multiplication schematic:
  do iw = 1,n_wfs
     ! Store matrix multiplication with coeffs in the
     ! "right" temporary buffer array:
     tmp_right(:) = 0.0
     do nu=1,no_u
        do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
           alpha = listh(ind)
           tmp_right(nu) = tmp_right(nu) +  &
                & (0.0 - sum(gradS(ind,:) * &
                & va_before_move(:,iaorb(nu)))) * &
                & coeffs(alpha,iw)
        end do
     end do

     ! Use "left" buffer array to compute center-right multiplication,
     ! then store the resulting column in the "right" buffer array:
     tmp_left(:) = 0.0
     do mu=1,no_u
        ! do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
        !    nu = listh(ind)
        !    tmp_left(mu) = tmp_left(mu) + dscf_hat(ind) * tmp_right(nu)
        ! end do
        do nu=1,no_l
           tmp_left(mu) = tmp_left(mu) + 0.5_dp * Dfull(mu,nu,1) * tmp_right(nu)
        end do
     end do
     tmp_right(:) = tmp_left(:)

     ! Compute the "left" matrix multiplication,
     ! then add to the resulting array:
     do lambda=1,no_u
        do mu=1,no_u
           acc_left = 0.0
           do beta=1,no_l     ! Dfull is dense now
           ! do ind = (listhptr(lambda)+1),&
           !      & listhptr(lambda) + numh(lambda)
           !    beta = listh(ind)
              do sec_ind = (listhptr(beta)+1),&
                   & listhptr(beta) + numh(beta)
                 sec_col = listh(sec_ind)
                 if ( sec_col .eq. mu ) then
                    if (beta.eq.lambda) then
                       acc_left = acc_left + &
                            ! &(1.0 - dscf_hat(ind)) * S(sec_ind)
                            &(1.0_dp  - 0.5_dp * Dfull(lambda,beta,1)) * S(sec_ind)
                    else
                       acc_left = acc_left + &
                            ! &(0.0 - dscf_hat(ind)) * S(sec_ind)
                            &(0.0_dp - 0.5_dp * Dfull(lambda,beta,1)) * S(sec_ind)
                    end if
                    exit
                 end if
              end do
           end do

           ! add dot-product element to the resulting array:
           psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw) +&
                & acc_left * tmp_right(mu)
        end do
     end do
  end do                        ! <-- end of 3rd matmul row
  ! call timer("explicit_matmul_Jks", 2)
  ! call timer("explicit_matmul_Jks", 3)

  deallocate(tmp_left)
  deallocate(tmp_right)
end subroutine compute_psi_dot_c


subroutine compute_Jks ()

  use sparse_matrices, only: H, S
  use ks_flux_data, only: psi_hat_c, psi_dot_c, ks_flux_Jks

  use m_eo, only: eo

  use,intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

  integer :: idx
  integer :: no_occ_wfs
  integer :: iw, mu, nu, j, ind, col

  !> Selected element of matrix H
  !> (only a factor, not the whole sandwitch!)

  real(dp):: H_el_S

  allocate(psi_hat_c(no_u,no_l,3))
  allocate(psi_dot_c(no_u,no_l))

  no_occ_wfs = ceiling(qtot/2)

  call compute_psi_hat_c (psi_hat_c, no_occ_wfs)
  call compute_psi_dot_c (psi_dot_c, no_occ_wfs)

  ks_flux_Jks(:) = 0.0          ! Init result to zeros outside main loop

  ! Parallel operation note: we need psi_hat_c and psi_dot_c coefficients for all
  ! occupied wavefunctions in all processors.
  !
  do iw = 1,no_occ_wfs

     do mu = 1,no_l
        do j = 1,numh(mu)
           ind = j + listhptr(mu)
           col = listh(ind)
           nu = indxuo(col)

           H_el_S = H(ind,1) + eo(iw,1,1) * S(ind)   ! We assume Gamma point with no
                                                    ! aux supercell
                                                    ! Otherwise, fold to H(k=0), S(k=0)
           do idx = 1,3

              ks_flux_Jks(idx) = ks_flux_Jks(idx) + &
                    psi_hat_c(mu,iw,idx) * H_el_s * psi_dot_c(nu,iw)

           end do
        end do
     end do
  end do

  ks_flux_Jks(:) = 2.0_dp * ks_flux_Jks(:)   ! For spin
  !
  ! For parallel operation, now reduce ks_flux_Jks over all processors
  !
  Print *, "[Jks] ", ks_flux_Jks

  deallocate(psi_hat_c)
  deallocate(psi_dot_c)

end subroutine compute_Jks


end module ks_flux_procs
