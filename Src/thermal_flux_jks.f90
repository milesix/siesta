module thermal_flux_jks
  !! Contains procedures for computation of J_ks. See notes on components
  !! of the Kohn-Sham flux for thermal transport. Brief reference for
  !! initial implementation given below (discussed objects could be renamed):
  !!
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
  !!

  use precision,       only: dp, grid_p
  use atomlist,        only: no_u, no_l, indxuo, iaorb, qtot, amass
  use siesta_geom,     only: na_u
  use densematrix,     only: psi ! To get at the c^i_\mu coefficients of the wavefunctions
  use sparse_matrices, only: listh, listhptr, numh
  ! use sparse_matrices, only: H, S, Dscf
  use sparse_matrices, only: gradS, Rmat

  use thermal_flux_data

  implicit none

  private
  public :: compute_jks

  real(dp), allocatable       :: psi_hat_c(:,:,:)
  real(dp), allocatable       :: psi_dot_c(:,:)

contains

  subroutine compute_psi_hat_c (n_wfs)
    integer, intent(in) :: n_wfs
    !! Number of wavefunctions to process

    !> Coordinate index regarding the selected dimension of Rmat.
    !> It will define X/Y/Z component of resulting Jks.
    integer :: idx
    integer :: nu, mu, alpha, iw, ij, k, col, ind

    real(dp), allocatable :: tmp_g(:)

    allocate (tmp_g(no_l))

    ! computing and storage of every of the 3 components
    ! should be available to order in the .fdf

    psi_hat_c(:,:,:) = 0.0        ! Init result to zeros outside main loop

    do idx=1,3
       do iw = 1,n_wfs

          tmp_g(:) = 0.0

          do nu = 1,no_l      ! g(nu) distributed...
             do ij = 1, numh(nu)
                k = listhptr(nu) + ij
                col = listh(k)

                tmp_g(nu) = tmp_g(nu) + psi_base(col,iw) * Rmat(k,idx)
             end do
          end do

          ! Use dense, complete, DM, but note that, for the spinless case, DM = 2 \sum_occ { |psi><psi| }
          ! We need then to use a factor of 1/2

          do alpha = 1,no_u
             do nu = 1,no_l   ! *** WARNING: local vs global?
                !! nu_g = LocalToGlobalOrb(nu...)
                !!           if (alpha == nu) then    ! alpha == nu_g   ; rest depending on
                ! whether Dfull is distributed or not
                psi_hat_c(alpha,iw,idx) = psi_hat_c(alpha,iw,idx) + &
                     tmp_g(nu) * (Sinv(nu,alpha) - 0.5_dp * DM_save(nu,alpha,1,1)) ! <- only 1st spin component
                !!           else
                !!              psi_hat_c(alpha,iw,idx) = psi_hat_c(alpha,iw,idx) + &
                !!                   tmp_g(nu) * (- 0.5_dp * Dfull(nu,alpha,1)) ! <- only 1st spin component
                !!           end if

             end do
          end do
          ! MPI: reduce (sum) psi_hat_c() over processes

       end do ! iw
    end do ! idx

    deallocate(tmp_g)
  end subroutine compute_psi_hat_c



  subroutine compute_psi_dot_c (n_wfs)

    integer, intent(in) :: n_wfs
    !! Number of wavefunctions to process
    integer  :: iw, i, ix, ij
    integer  :: alpha, beta, lambda, nu, mu, ind
    real(dp) :: acc_left, norm, tr0, tr_dot, tr_sdot

    real(dp), allocatable :: tmp(:), tmp2(:), tmp3(:), tmp4(:,:)

    allocate (tmp(no_u))
    allocate (tmp2(no_u))
    allocate (tmp3(no_u))
    allocate (tmp4(3,no_u))

    psi_dot_c(:,:) = 0.0_dp          ! Init result to zeros outside main loop

    ! multiplication scheme A
    do iw = 1,n_wfs

       tmp(:) = 0.0_dp

       do nu = 1,no_u
          do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
             alpha = listh(ind)
             tmp(nu) = tmp(nu) + S_base(nu,alpha) * psi_base(alpha,iw)
          end do
       end do

       tmp2(:) = 0.0_dp

       do mu = 1,no_u
          do nu = 1,no_u
             tmp2(mu) = tmp2(mu) + 0.5_dp * DM_save(mu,nu,1,2) * tmp(nu) ! using derivative of DM
          end do
       end do

       tmp3(:) = 0.0_dp

       do beta = 1,no_u
          do ind = (listhptr(beta)+1), listhptr(beta) + numh(beta)
             mu = listh(ind)
             tmp3(beta) = tmp3(beta) + S_base(beta,mu) * tmp2(mu)
          end do
       end do

       do lambda=1,no_u
          do beta = 1,no_u
             psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw)&
                  & + (Sinv(lambda,beta) - 0.5_dp * DM_save(lambda,beta,1,1)) * tmp3(beta) ! using full DM
          end do
       end do
    end do
    ! end of part A


    ! multiplication scheme B
    do iw = 1,n_wfs

       tmp(:) = 0.0_dp

       do nu = 1,no_u
          do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
             alpha = listh(ind)
             tmp(nu) = tmp(nu) + S_base(nu,alpha) * psi_base(alpha,iw)
          end do
       end do

       tmp2(:) = 0.0_dp

       do mu = 1,no_u
          do nu = 1,no_u
             tmp2(mu) = tmp2(mu) + 0.5_dp * DM_save(mu,nu,1,1) * tmp(nu) ! using full DM
          end do
       end do

       tmp3(:) = 0.0_dp

       do beta = 1,no_u
          do ind = (listhptr(beta)+1), listhptr(beta) + numh(beta)
             mu = listh(ind)
             tmp3(beta) = tmp3(beta) + sum(gradS(:,ind)*va_in(:,iaorb(mu))) * tmp2(mu)
          end do
       end do

       do lambda=1,no_u
          do beta = 1,no_u
             psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw)&
                  & + (Sinv(lambda,beta) - 0.5_dp * DM_save(lambda,beta,1,1)) * tmp3(beta) ! using full DM
          end do
       end do
    end do
    ! end of part B


    ! multiplication scheme C
    do iw = 1,n_wfs

       tmp4(:,:) = 0.0_dp

       do nu = 1,no_u
          do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
             alpha = listh(ind)
             tmp4(:,nu) = tmp4(:,nu) + gradS(:,ind) * psi_base(alpha,iw)
          end do
       end do

       tmp2(:) = 0.0_dp

       do mu = 1,no_u
          do nu = 1,no_u
             tmp2(mu) = tmp2(mu) - 0.5_dp * DM_save(mu,nu,1,1) & ! using full DM
                  &* sum(va_in(:,iaorb(nu))*tmp4(:,nu))
          end do
       end do

       tmp3(:) = 0.0_dp

       do beta = 1,no_u
          do ind = (listhptr(beta)+1), listhptr(beta) + numh(beta)
             mu = listh(ind)
             tmp3(beta) = tmp3(beta) + S_base(beta,mu) * tmp2(mu)
          end do
       end do

       do lambda=1,no_u
          do beta = 1,no_u
             psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw)&
                  & + (Sinv(lambda,beta) - 0.5_dp * DM_save(lambda,beta,1,1)) * tmp3(beta) ! using full DM
          end do
       end do


    end do
    ! end of part C

    deallocate(tmp)
    deallocate(tmp2)
    deallocate(tmp3)
    deallocate(tmp4)

  end subroutine compute_psi_dot_c


  subroutine compute_jks ()
    integer :: no_occ_wfs

    allocate(psi_hat_c(no_u,no_l,3))
    allocate(psi_dot_c(no_u,no_l))

    no_occ_wfs = nint(qtot/2)

    call compute_psi_hat_c (no_occ_wfs)
    call compute_psi_dot_c (no_occ_wfs)

    deallocate(psi_hat_c)
    deallocate(psi_dot_c)

  end subroutine compute_jks

end module thermal_flux_jks
