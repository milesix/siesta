module ks_flux_procs
  !! Contains procedures for computation of J_ks.
  !! See notes on components of the Kohn-Sham flux for thermal transport.
  use precision, only: dp
  use atomlist, only: no_u, no_l, indxuo, iaorb, qtot
  use densematrix, only: psi ! To get at the c^i_\mu coefficients of the wavefunctions
  use sparse_matrices, only: listh, listhptr, numh

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
  !! these two are also sparse.
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

  ! allocate(psi_hat_c(no_u,no_l,3))

  allocate (tmp_g(no_u))

  dscf_hat => ks_flux_D(:,1)    !NOTE: Just first spin component
  ! computing and storage of every of the 3 components
  ! should be available to order in the .fdf

  psi_hat_c(:,:,:) = 0.0        ! Init result to zeros outside main loop

  DO idx=1,3
  rmatrix => ks_flux_Rmat(:,idx)

  do iw = 1,n_wfs
     tmp_g(:) = 0.0
     do nu = 1,no_u
        do mu = 1,no_l
          do ij = 1, numh(mu)
             k = listhptr(mu) + ij
             col = listh(k)
             if (col == nu) then
               tmp_g(nu) = tmp_g(nu) + coeffs(mu,iw) * rmatrix(k)
               exit
             endif
          end do
        end do
     end do

     do alpha = 1,no_u
        do nu = 1,no_l
          do ij = 1, numh(nu)
             k = listhptr(nu) + ij
             col = listh(k)
             if (col == alpha) then
               if (alpha == nu) then
                 psi_hat_c(alpha,iw,idx) = psi_hat_c(alpha,iw,idx) + &
                      tmp_g(nu) * (1.0 - dscf_hat(k))
               else
                 psi_hat_c(alpha,iw,idx) = psi_hat_c(alpha,iw,idx) + &
                      tmp_g(nu) * (- dscf_hat(k))
               end if
               exit
             endif
          end do
        end do
     end do
  end do
  END DO

  deallocate(tmp_g)

end subroutine compute_psi_hat_c


subroutine compute_psi_dot_c (psi_dot_c, n_wfs)
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
        do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
           nu = listh(ind)
           tmp_left(mu) = tmp_left(mu)+Dscf_deriv(ind,1)*tmp_right(nu)
           !NOTE: Only first spin component here now--^
        end do
     end do
     tmp_right(:) = tmp_left(:)

     ! Compute the "left" matrix multiplication,
     ! then add to the resulting array:
     do lambda=1,no_u
        do mu=1,no_u
           acc_left = 0.0
           do ind = (listhptr(lambda)+1),&
                & listhptr(lambda) + numh(lambda)
              beta = listh(ind)
              do sec_ind = (listhptr(beta)+1),&
                   & listhptr(beta) + numh(beta)
                 sec_col = listh(sec_ind)
                 if ( sec_col .eq. mu ) then
                    if (beta.eq.lambda) then
                       acc_left = acc_left + &
                            &(1.0 - dscf_hat(ind)) * S(sec_ind)
                    else
                       acc_left = acc_left + &
                            &(0.0 - dscf_hat(ind)) * S(sec_ind)
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
        do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
           nu = listh(ind)
           tmp_left(mu) = tmp_left(mu) + dscf_hat(ind) * tmp_right(nu)
        end do
     end do
     tmp_right(:) = tmp_left(:)

     ! Compute the "left" matrix multiplication,
     ! then add to the resulting array:
     do lambda=1,no_u
        do mu=1,no_u
           acc_left = 0.0
           do ind = (listhptr(lambda)+1),&
                & listhptr(lambda) + numh(lambda)
              beta = listh(ind)
              do sec_ind = (listhptr(beta)+1),&
                   & listhptr(beta) + numh(beta)
                 sec_col = listh(sec_ind)
                 if ( sec_col .eq. mu ) then
                    if (beta.eq.lambda) then
                       acc_left = acc_left + &
                            &(1.0 - dscf_hat(ind)) * &
                            &(sum(gradS(sec_ind,:) * &
                            & va_before_move(:,iaorb(mu))))
                    else
                       acc_left = acc_left + &
                            &(0.0 - dscf_hat(ind)) * &
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
        do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
           nu = listh(ind)
           tmp_left(mu) = tmp_left(mu) + dscf_hat(ind) * tmp_right(nu)
        end do
     end do
     tmp_right(:) = tmp_left(:)

     ! Compute the "left" matrix multiplication,
     ! then add to the resulting array:
     do lambda=1,no_u
        do mu=1,no_u
           acc_left = 0.0
           do ind = (listhptr(lambda)+1),&
                & listhptr(lambda) + numh(lambda)
              beta = listh(ind)
              do sec_ind = (listhptr(beta)+1),&
                   & listhptr(beta) + numh(beta)
                 sec_col = listh(sec_ind)
                 if ( sec_col .eq. mu ) then
                    if (beta.eq.lambda) then
                       acc_left = acc_left + &
                            &(1.0 - dscf_hat(ind)) * S(sec_ind)
                    else
                       acc_left = acc_left + &
                            &(0.0 - dscf_hat(ind)) * S(sec_ind)
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

  use sparse_matrices, only: H
  use ks_flux_data, only: psi_hat_c, psi_dot_c, ks_flux_Jks

  use m_eo, only: eo

  use,intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

  integer :: idx
  integer :: no_occ_wfs
  integer :: iw, alpha, beta, j, ind, col

  !> Selected element of matrix H
  !> (only a factor, not the whole sandwitch!)
  real(dp):: H_el

  !> First dimension: coeff index,
  !> Second dimension: wavefunction index (only no_l of them)
  real(dp), pointer :: coeffs(:,:)
  real(dp), allocatable :: res_array(:)

  call c_f_pointer(c_loc(psi(1)), coeffs, [no_u, no_l])

  allocate(psi_hat_c(no_u,no_l,3))
  allocate(psi_dot_c(no_u,no_l))

  no_occ_wfs = ceiling(qtot/2)

  call compute_psi_hat_c (psi_hat_c, no_occ_wfs)
  call compute_psi_dot_c (psi_dot_c, no_occ_wfs)

  ks_flux_Jks(:) = 0.0          ! Init result to zeros outside main loop

  do iw = 1,no_occ_wfs
     do alpha = 1,no_u
        do j = 1,numh(alpha)
           ind = j + listhptr(alpha)
           col = listh(ind)
           beta = indxuo(col)

           if (alpha.eq.beta) then
              H_el = H(ind,1) + eo(iw,1,1) ! <- correct?
           else
              H_el = H(ind,1) ! <- NOTE: one spin component
           endif

           do idx = 1,3

              ks_flux_Jks(idx) = ks_flux_Jks(idx) + &
                   & psi_hat_c(alpha,iw,idx) * &
                   & H_el * &
                   & psi_dot_c(beta,iw)
           end do
        end do
     end do
  end do

  Print *, "[Jks] ", ks_flux_Jks

  deallocate(psi_hat_c)
  deallocate(psi_dot_c)

end subroutine compute_Jks


end module ks_flux_procs