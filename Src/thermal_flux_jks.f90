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
  !! use sparse_matrices, only: H, S, Dscf
  !! use sparse_matrices, only: gradS, Rmat

  use thermal_flux_data     ! H_base, eo_base, gradS_base, Sinv, etc

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
    real(dp) :: tr_hat, tr_shat

    real(dp), allocatable :: tmp_g(:)
    real(dp), allocatable :: work(:), amat(:,:)
    integer, allocatable :: ipiv(:)
    real(dp), parameter :: alpha_reg = 0.001_dp ! To regularize (H-e)^-1
    integer :: i, j, lwork, info

    allocate (tmp_g(no_l))
    allocate(work(no_u), ipiv(no_u), amat(no_u,no_u))

    ! computing and storage of every of the 3 components
    ! should be available to order in the .fdf

    psi_hat_c(:,:,:) = 0.0_dp        ! Init result to zeros outside main loop

    ! Compute first P_c [H,r] |Psi_iw>
    do idx=1,3
       do iw = 1,n_wfs

          tmp_g(:) = 0.0

          do nu = 1,no_l      ! g(nu) distributed...
             do ij = 1, numh(nu)
                k = listhptr(nu) + ij
                col = listh(k)

                !! tmp_g(nu) = tmp_g(nu) + psi_base(col,iw) * Rmat_base(k,idx)
                ! New equations: Use [H,r] trick and linear equations.
                ! So far without non-local potential term, so Gamma_mu_nu in notes
                ! is just gradS_mu_nu
                tmp_g(nu) = tmp_g(nu) + psi_base(col,iw) * gradS_base(idx,k)
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

          !
          ! *** NOTE: This should be moved after the real computation of psi_hat_c
          !
          if (gk_setup%verbose_output) then
             !TEST-hat
             tr_shat = 0.0_dp
             tr_hat = 0.0_dp
             do mu = 1,no_u
                ! do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
                !    nu = listh(ind)
                do nu = 1,no_u
                   tr_shat = tr_shat + S_base(mu,nu) * psi_base(mu,iw) * psi_hat_c(nu,iw,idx)
                   tr_hat = tr_hat + S_base(mu,nu) * psi_hat_c(mu,iw,idx) * psi_hat_c(nu,iw,idx)
                end do
             end do

             print *, "[TEST-hat]", iw, idx, tr_hat, tr_shat
             !TEST END
          end if

       end do ! iw
    end do ! idx

    deallocate(tmp_g)

    do iw = 1,n_wfs

    ! Build linear equations
    ! (H - eps_i + alpha*P_v) (Psi_hat) = psi_hat_c above
    ! Use dsysv for AX=B, with the solution X overwriting B (just what we need)
    !
    ! A = H_base - eo_base(i) 1 + alpha*(0.5_dp Dfull)

       do i = 1, no_u
          do j = 1, no_u
             amat(i,j) = H_base(i,j) + alpha_reg*DM_save(i,j,1,1)
          enddo
          amat(i,i) = amat(i,i) - eo_base(iw)
       enddo

!!$subroutine dsysv(	character 	UPLO,
!!$integer 	N,
!!$integer 	NRHS,
!!$double precision, dimension( lda, * ) 	A,
!!$integer 	LDA,
!!$integer, dimension( * ) 	IPIV,
!!$double precision, dimension( ldb, * ) 	B,
!!$integer 	LDB,
!!$double precision, dimension( * ) 	WORK,
!!$integer 	LWORK,
!!$integer 	INFO
!!$)
       ! This will do idx=1..3 in one go.
       ! Note the computation for the leading dimension of B. The size of
       ! the second dimension of psi_hat_c is n_wfs.
       lwork = 4*no_u
       call dsysv('U', no_u, 3, amat, no_u, ipiv, psi_hat_c(1,iw,1), no_u*n_wfs, &
            work, lwork, info)
       if (info < 0) then
          write(0,*) "Argument had an illegal value: ", -info
          call die("Error in dsysv call")
       else if (info > 0) then
          call die("A is singular in dsysv call")
       endif

    enddo ! iw

    deallocate(ipiv, work, amat)

  end subroutine compute_psi_hat_c



  subroutine compute_psi_dot_c (n_wfs)

    integer, intent(in) :: n_wfs
    !! Number of wavefunctions to process
    integer  :: iw, i, ix, ij
    integer  :: alpha, beta, lambda, nu, mu, ind
    real(dp) :: tr0, tr_dot, tr_sdot

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
          ! do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
          !    alpha = listh(ind)
          do alpha = 1,no_u
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
          ! do ind = (listhptr(beta)+1), listhptr(beta) + numh(beta)
          !    mu = listh(ind)
          do mu = 1,no_u
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
          ! do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
          !    alpha = listh(ind)
          do alpha = 1,no_u
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
          ! do mu = 1,no_u
             tmp3(beta) = tmp3(beta) + sum(gradS_base(:,ind)*va_in(:,iaorb(mu))) * tmp2(mu)
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
          ! do alpha = 1,no_u
             tmp4(:,nu) = tmp4(:,nu) + gradS_base(:,ind) * psi_base(alpha,iw)
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
          ! do ind = (listhptr(beta)+1), listhptr(beta) + numh(beta)
          !    mu = listh(ind)
          do mu = 1,no_u
             tmp3(beta) = tmp3(beta) + S_base(beta,mu) * tmp2(mu)
          end do
       end do

       do lambda=1,no_u
          do beta = 1,no_u
             psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw)&
                  & + (Sinv(lambda,beta) - 0.5_dp * DM_save(lambda,beta,1,1)) * tmp3(beta) ! using full DM
          end do
       end do

       if (gk_setup%verbose_output) then
          !TEST-dot
          tr_sdot = 0.0_dp
          tr_dot = 0.0_dp
          tr0 = 0.0_dp
          do mu = 1,no_u
             ! do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
             !    nu = listh(ind)
             do nu = 1,no_u
                tr_sdot = tr_sdot + S_base(mu,nu) * psi_base(mu,iw) * psi_dot_c(nu,iw)
                tr0 = tr0 + S_base(mu,nu) * psi_base(mu,iw) * psi_base(nu,iw)
                tr_dot = tr_dot + S_base(mu,nu) * psi_dot_c(mu,iw) * psi_dot_c(nu,iw)
             end do
          end do

          print *, "[TEST-dot]", iw, tr0, tr_dot, tr_sdot
       end if

    end do
    ! end of part C

    deallocate(tmp)
    deallocate(tmp2)
    deallocate(tmp3)
    deallocate(tmp4)

  end subroutine compute_psi_dot_c


  subroutine compute_jks ()
    integer :: no_occ_wfs
    integer :: idx, iw, mu, nu, j, ind, col

    !> Selected element of matrix H
    !> (only a factor, not the whole sandwitch!)
    real(dp):: H_el_S, S_el_S

    ! real(dp):: ks_flux_Jks(3)
    real(dp):: ks_flux_Jks_A(3), ks_flux_Jks_B(3), ks_flux_Jele(3)

    allocate(psi_hat_c(no_u,no_l,3))
    allocate(psi_dot_c(no_u,no_l))

    no_occ_wfs = nint(qtot/2)

    call compute_psi_hat_c (no_occ_wfs)
    call compute_psi_dot_c (no_occ_wfs)

    ! ks_flux_Jks(:) = 0.0          ! Init result to zeros outside main loop
    ks_flux_Jks_A(:) = 0.0          ! Init result to zeros outside main loop
    ks_flux_Jks_B(:) = 0.0          ! Init result to zeros outside main loop
    ks_flux_Jele(:) = 0.0          ! Init result to zeros outside main loop

    ! Parallel operation note: we need psi_hat_c and psi_dot_c coefficients for all
    ! occupied wavefunctions in all processors.
    !
    do iw = 1,no_occ_wfs

       do mu = 1,no_l
          do nu = 1,no_l

             ! H_el_S = H_base(mu,nu) + eo_base(iw) * S_base(mu,nu)   ! These are (k=0) matrices
             H_el_S = H_base(mu,nu)
             S_el_S = eo_base(iw) * S_base(mu,nu)   ! These are (k=0) matrices

             do idx = 1,3
                ! ks_flux_Jks(idx) = ks_flux_Jks(idx) + &
                !      psi_hat_c(mu,iw,idx) * H_el_s * psi_dot_c(nu,iw)

                ks_flux_Jks_A(idx) = ks_flux_Jks_A(idx) + &
                     psi_hat_c(mu,iw,idx) * H_el_S * psi_dot_c(nu,iw)

                ks_flux_Jks_B(idx) = ks_flux_Jks_B(idx) + &
                     psi_hat_c(mu,iw,idx) * S_el_S * psi_dot_c(nu,iw)

                ks_flux_Jele(idx) = ks_flux_Jele(idx) + &
                     psi_hat_c(mu,iw,idx) * S_base(mu,nu) * psi_dot_c(nu,iw)
             end do

          end do
       end do

    end do

    ! For parallel operation, reduce ks_flux_J* over all processors...

    ! Save results for the Kohn-Sham flux component:
    ! ks_flux_Jks(:) = 2.0_dp * ks_flux_Jks(:)   ! x2 For spin
    gk_results%Jks_A(:) = 2.0_dp * ks_flux_Jks_A(:)   ! x2 For spin
    gk_results%Jks_B(:) = 2.0_dp * ks_flux_Jks_B(:)   ! x2 For spin
    gk_results%Jks(:)   = 2.0_dp * (ks_flux_Jks_A(:) + ks_flux_Jks_B(:))
    gk_results%Jele(:)  = 2.0_dp * 2.0_dp * ks_flux_Jele(:)   ! Why the multiplication by 2*2??
    !
    deallocate(psi_hat_c)
    deallocate(psi_dot_c)

  end subroutine compute_jks

end module thermal_flux_jks
