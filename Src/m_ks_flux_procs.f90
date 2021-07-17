module ks_flux_procs
  !! Contains procedures for computation of J_ks.
  !! See notes on components of the Kohn-Sham flux for thermal transport.
  use precision, only: sp, dp
  use atomlist, only: no_u, no_l, indxuo, iaorb, qtot, amass
  use densematrix, only: psi ! To get at the c^i_\mu coefficients of the wavefunctions
  use sparse_matrices, only: listh, listhptr, numh
  use m_virtual_step_data, only: Dfull, Dderiv, istep_vmd
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
  use ks_flux_data, only: ks_flux_S
  use m_virtual_step_data, only: Sinv
  
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
  integer :: nu, mu, alpha, iw, ij, k, col, ind
  real(dp) :: norm, tr_hat, tr_shat
  
  !> First dimension: coeff index,
  !> Second dimension: wavefunction index (only no_l of them)
  real(dp), pointer :: coeffs(:,:)
  real(dp), pointer :: overlap(:)
  real(dp), allocatable :: tmp_g(:)

  call c_f_pointer(c_loc(psi(1)), coeffs, [no_u, no_l]) ! for new compilers
  ! explicit pointing to psi(1) -> psi; avoid gcc<=4.8

  allocate (tmp_g(no_l))

  dscf_hat => ks_flux_D(:,1)    !NOTE: Just first spin component
  ! computing and storage of every of the 3 components
  ! should be available to order in the .fdf

  psi_hat_c(:,:,:) = 0.0        ! Init result to zeros outside main loop

  overlap => ks_flux_S
  DO idx=1,3
     rmatrix => ks_flux_Rmat(:,idx)


  do iw = 1,n_wfs
     tmp_g(:) = 0.0
     do nu = 1,no_l      ! g(nu) distributed...
        ! This can be simplified
        ! Remove this loop...
        do mu = 1,no_u
          !
          do ij = 1, numh(nu)
             k = listhptr(nu) + ij
             col = listh(k)
             ! ... and simply do
             ! tmp_g(nu) = tmp_g(nu) + coeffs(col,iw) * rmatrix(k)
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
        do nu = 1,no_l   ! *** WARNING: local vs global?
           !! nu_g = LocalToGlobalOrb(nu...)
    !!           if (alpha == nu) then    ! alpha == nu_g   ; rest depending on
                                    ! whether Dfull is distributed or not
              psi_hat_c(alpha,iw,idx) = psi_hat_c(alpha,iw,idx) + &
                   tmp_g(nu) * (sinv(nu,alpha) - 0.5_dp * Dfull(nu,alpha,1)) ! <- only 1st spin component
    !!           else
    !!              psi_hat_c(alpha,iw,idx) = psi_hat_c(alpha,iw,idx) + &
    !!                   tmp_g(nu) * (- 0.5_dp * Dfull(nu,alpha,1)) ! <- only 1st spin component
    !!           end if

        end do
     end do
     ! MPI: reduce (sum) psi_hat_c() over processes
     
!!$     ! Normalize here
!!$     ! norm = c*S*c
!!$     tmp_g(:) = 0.0
!!$     do nu = 1,no_l
!!$        do ij = 1, numh(nu)
!!$             k = listhptr(nu) + ij
!!$             col = listh(k)
!!$             tmp_g(nu) = tmp_g(nu) + overlap(k) * psi_hat_c(col,iw,idx)
!!$        end do
!!$     end do
!!$     norm = 0.0_dp
!!$     do mu = 1, no_l
!!$        norm = norm + psi_hat_c(mu,iw,idx) * tmp_g(mu)
!!$     enddo
!!$!!!     psi_hat_c(:,iw,idx) = psi_hat_c(:,iw,idx) / (sqrt(norm) + 1.0e-8_dp)

     !TEST
     tr_shat = 0.0_dp
     tr_hat = 0.0_dp
     do mu = 1,no_u
        do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
           nu = listh(ind)
           tr_shat = tr_shat + overlap(ind) * coeffs(mu,iw) * psi_hat_c(nu,iw,idx)
           tr_hat = tr_hat + overlap(ind) * psi_hat_c(mu,iw,idx) * psi_hat_c(nu,iw,idx)
        end do
     end do

     print *, "[TEST-hat]", iw, idx, tr_hat, tr_shat

!TEST END
  end do ! iw
  END DO ! idx

  deallocate(tmp_g)

end subroutine compute_psi_hat_c

!
!  This routine is very complex. It has to be re-checked
!
subroutine compute_psi_dot_c (psi_dot_c, n_wfs)
  use siesta_geom,  only : na_u
  use siesta_options, only : virtual_dt
  ! Take proper velocities from `BASE` step:
  use m_virtual_step_data, only: va_before_move

  ! Entities stored during the BASE MD-step:
  use ks_flux_data, only: ks_flux_D, ks_flux_S, ks_flux_gS

  ! Output: coefficients of psi_dot_c, stored in module ks_flux_data
  ! use ks_flux_data, only: psi_dot_c

  use siesta_options, only: virtual_md_Jks_A
  use siesta_options, only: virtual_md_Jks_B
  use siesta_options, only: virtual_md_Jks_C

  use,intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

  implicit none

  real(dp), allocatable, intent(inout) :: psi_dot_c(:,:)
  integer, intent(in) :: n_wfs
  !! Number of wavefunctions to process

  real(dp), dimension(:),   pointer :: dscf_hat => null()
  real(dp), dimension(:,:), pointer :: gradS => null()
  real(dp), dimension(:),   pointer :: S => null()

  integer  :: iw, i, ix, ij
  integer  :: alpha, beta, lambda, nu, mu
  ! integer  :: ind, col, sec_ind, sec_col
  real(dp) :: acc_left, norm, tr0, tr_dot, tr_sdot

  !> First dimension: coeff index,
  !> Second dimension: wavefunction index (only no_l of them)
  real(dp), pointer :: coeffs(:,:)

  integer :: ind, ind2
  integer :: k, col
  real(dp), allocatable :: tmp(:), tmp2(:), tmp3(:), tmp4(:,:)
  real(dp) :: pcm(3), massi, mtot
  ! real(dp), allocatable :: tmp_left(:)
  ! real(dp), allocatable :: tmp_right(:)

  real(dp), allocatable :: tmp_g(:)

  allocate (tmp_g(no_l))

  call c_f_pointer(c_loc(psi(1)), coeffs, [no_u, no_l]) ! for new compilers
  ! explicit pointing to psi(1) -> psi
  ! avoid gcc<=4.8

  allocate (tmp(no_u))
  allocate (tmp2(no_u))
  allocate (tmp3(no_u))
  allocate (tmp4(3,no_u))
  ! allocate (tmp_left(no_u))
  ! allocate (tmp_right(no_u))

  ! dscf_hat => ks_flux_D(:,1)         !NOTE: Just first spin component
  S => ks_flux_S(:)
  gradS => ks_flux_gS(:,:)
  ! computing and storage of every of the 3 components
  ! should be available to order in the .fdf

  psi_dot_c(:,:) = 0.0_dp          ! Init result to zeros outside main loop

  ! call timer("explicit_matmul_Jks", 1)


  ! multiplication scheme A
  do iw = 1,n_wfs

     tmp(:) = 0.0_dp

     do nu = 1,no_u
        do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
           alpha = listh(ind)
           tmp(nu) = tmp(nu) + S(ind) * coeffs(alpha,iw)
        end do
     end do

     tmp2(:) = 0.0_dp
! Dscf_deriv(ind,1)
     do mu = 1,no_u
        do nu = 1,no_u
        ! do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
        !    nu = listh(ind)
           tmp2(mu) = tmp2(mu) + 0.5_dp * Dderiv(mu,nu,1) * tmp(nu)
        end do
     end do

     !! WARNING: Temporary change!!!
     tmp2(:) = 0.0_dp
     
     tmp3(:) = 0.0_dp

     do beta = 1,no_u
        do ind = (listhptr(beta)+1), listhptr(beta) + numh(beta)
           mu = listh(ind)
           tmp3(beta) = tmp3(beta) + S(ind) * tmp2(mu)
        end do
     end do

     do lambda=1,no_u
        do beta = 1,no_u
           if (beta.eq.lambda) then
              psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw)&
                   & + (1.0_dp - 0.5_dp * Dfull(lambda,beta,1)) * tmp3(beta)&
                   & * virtual_md_Jks_A
           else
              psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw)&
                   & + (0.0_dp - 0.5_dp * Dfull(lambda,beta,1)) * tmp3(beta)&
                   & * virtual_md_Jks_A
           end if
        end do
     end do
  end do                        ! end of part A


  ! multiplication scheme B
  do iw = 1,n_wfs

     tmp(:) = 0.0_dp

     do nu = 1,no_u
        do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
           alpha = listh(ind)
           tmp(nu) = tmp(nu) + S(ind) * coeffs(alpha,iw)
        end do
     end do

     tmp2(:) = 0.0_dp

     do mu = 1,no_u
        do nu = 1,no_u
           tmp2(mu) = tmp2(mu) + 0.5_dp * Dfull(mu,nu,1) * tmp(nu)
        end do
     end do

     tmp3(:) = 0.0_dp

     do beta = 1,no_u
        do ind = (listhptr(beta)+1), listhptr(beta) + numh(beta)
           mu = listh(ind)
           tmp3(beta) = tmp3(beta) + sum(gradS(:,ind)*va_before_move(:,iaorb(mu))) * tmp2(mu)
        end do
     end do

     do lambda=1,no_u
        do beta = 1,no_u
           if (beta.eq.lambda) then
              psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw)&
                   & + (1.0_dp - 0.5_dp * Dfull(lambda,beta,1)) * tmp3(beta)&
                   & * virtual_md_Jks_B
           else
              psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw)&
                   & + (0.0_dp - 0.5_dp * Dfull(lambda,beta,1)) * tmp3(beta)&
                   & * virtual_md_Jks_B
           end if
        end do
     end do
  end do                        ! end of part B


  ! multiplication scheme C
  do iw = 1,n_wfs

     tmp4(:,:) = 0.0_dp

     do nu = 1,no_u
        do ind = (listhptr(nu)+1), listhptr(nu) + numh(nu)
           alpha = listh(ind)
           tmp4(:,nu) = tmp4(:,nu) + gradS(:,ind) * coeffs(alpha,iw)
        end do
     end do
! va_before_move(:,iaorb(alpha))
     tmp2(:) = 0.0_dp

     do mu = 1,no_u
        do nu = 1,no_u
           tmp2(mu) = tmp2(mu) - 0.5_dp * Dfull(mu,nu,1) &
                &* sum(va_before_move(:,iaorb(nu))*tmp4(:,nu))
        end do
     end do

     tmp3(:) = 0.0_dp

     do beta = 1,no_u
        do ind = (listhptr(beta)+1), listhptr(beta) + numh(beta)
           mu = listh(ind)
           tmp3(beta) = tmp3(beta) + S(ind) * tmp2(mu)
        end do
     end do

     do lambda=1,no_u
        do beta = 1,no_u
           if (beta.eq.lambda) then
              psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw)&
                   & + (1.0_dp - 0.5_dp * Dfull(lambda,beta,1)) * tmp3(beta)&
                   & * virtual_md_Jks_C
           else
              psi_dot_c(lambda,iw) = psi_dot_c(lambda,iw)&
                   & + (0.0_dp - 0.5_dp * Dfull(lambda,beta,1)) * tmp3(beta)&
                   & * virtual_md_Jks_C
           end if
        end do
     end do

!!$     ! Normalize here
!!$     ! norm = c*S*c
!!$     tmp_g(:) = 0.0
!!$     do nu = 1,no_l
!!$        do ij = 1, numh(nu)
!!$             k = listhptr(nu) + ij
!!$             col = listh(k)
!!$             tmp_g(nu) = tmp_g(nu) + S(k) * psi_dot_c(col,iw)
!!$        end do
!!$     end do
!!$     norm = 0.0_dp
!!$     do mu = 1, no_l
!!$        norm = norm + psi_dot_c(mu,iw) * tmp_g(mu)
!!$     enddo
!!$!!     psi_dot_c(:,iw) = psi_dot_c(:,iw) / (sqrt(norm) + 1.0e-8_dp)

     !TEST
     tr_sdot = 0.0_dp
     tr_dot = 0.0_dp
     tr0 = 0.0_dp
     do mu = 1,no_u
        do ind = (listhptr(mu)+1), listhptr(mu) + numh(mu)
           nu = listh(ind)
           tr_sdot = tr_sdot + S(ind) * coeffs(mu,iw) * psi_dot_c(nu,iw)
           tr0 = tr0 + S(ind) * coeffs(mu,iw) * coeffs(nu,iw)
           tr_dot = tr_dot + S(ind) * psi_dot_c(mu,iw) * psi_dot_c(nu,iw)
        end do
     end do

     print *, "[TEST-dot]", iw, tr0, tr_dot, tr_sdot

!TEST END

  end do                        ! end of part C

  ! call timer("explicit_matmul_Jks", 2)
  ! call timer("explicit_matmul_Jks", 3)

  deallocate(tmp)
  deallocate(tmp2)
  deallocate(tmp3)
  deallocate(tmp4)
  deallocate(tmp_g)

  ! deallocate(tmp_left)
  ! deallocate(tmp_right)
end subroutine compute_psi_dot_c


subroutine compute_Jks ()

  ! use sparse_matrices, only: H
  ! use sparse_matrices, only: S

  use ks_flux_data, only: ks_flux_S, ks_flux_H
  use ks_flux_data, only: psi_hat_c, psi_dot_c, ks_flux_Jks, ks_flux_Jele

  use siesta_options, only: virtual_md_verbose
  use m_eo, only: eo

  use,intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

  integer :: idx
  integer :: no_occ_wfs
  integer :: iw, mu, nu, j, ind, col

  integer :: iu
  external io_assign, io_close

  real(dp), dimension(:),   pointer :: S => null()
  real(dp), dimension(:,:), pointer :: H => null()

  !> Selected element of matrix H
  !> (only a factor, not the whole sandwitch!)

  real(dp):: H_el_S

  character(len=1024) :: fname

  allocate(psi_hat_c(no_u,no_l,3))
  allocate(psi_dot_c(no_u,no_l))

  no_occ_wfs = nint(qtot/2)

  call compute_psi_hat_c (psi_hat_c, no_occ_wfs)
  call compute_psi_dot_c (psi_dot_c, no_occ_wfs)

  S => ks_flux_S(:)
  H => ks_flux_H(:,:)

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


  ks_flux_Jele(:) = 0.0          ! Init result to zeros outside main loop

  do iw = 1,no_occ_wfs

     do mu = 1,no_l
        do j = 1,numh(mu)
           ind = j + listhptr(mu)
           col = listh(ind)
           nu = indxuo(col)

           do idx = 1,3

              ks_flux_Jele(idx) = ks_flux_Jele(idx) + &
                    psi_hat_c(mu,iw,idx) * S(ind) * psi_dot_c(nu,iw)

           end do
        end do
     end do
  end do

  ks_flux_Jele(:) = 2.0_dp * 2.0_dp * ks_flux_Jele(:)

  Print *, "[Jele] ", ks_flux_Jele

  if (virtual_md_verbose) then
  ! Output of `pseudo-wfs-objects`:

     write(fname, "(A8,I3.3,A5)") "Psi_dot_", istep_vmd, ".WFSX"
     call write_ks_psi(psi_dot_c, fname)

     do idx=1,3
        write(fname, "(A8,I1,A1,I3.3,A5)")&
             & "Psi_hat-", idx, "_", istep_vmd, ".WFSX"
        call write_ks_psi(psi_hat_c(:,:,idx), fname)
     end do
  end if

  deallocate(psi_hat_c)
  deallocate(psi_dot_c)

end subroutine compute_Jks


subroutine write_ks_psi(ks_psi, fname)
  !! `WFSX`-file format follows `writewave.F`

  use kpoint_scf_m, only : kpoint_scf
  use atmfuncs,     only : symfio, cnfigfio, labelfis, nofis
  use atomlist,     only : iaorb, iphorb
  use siesta_geom,  only : isa
  use m_eo,         only : eo
  use units,        only : eV

  real(SP), dimension(:), allocatable :: aux   !! NOTE SP

  real(dp), intent(in) :: ks_psi(:,:)
  character(len=1024), intent(in)   :: fname

  integer :: nuotot
  integer :: iu, j
  integer :: iw, idx, no_occ_wfs
  external io_assign, io_close

  nuotot = no_u
  no_occ_wfs = nint(qtot/2)

  allocate(aux(nuotot))         ! 1-d since 1 spin component now

  call io_assign(iu)
  open(iu,file=trim(fname),form="unformatted",&
       & action='write',status='unknown')
  ! Writing header
  write (iu) 1, 1               ! nk, gamma
  write (iu) 1                  ! spin%H
  write (iu) nuotot
  write(iu) (iaorb(j),labelfis(isa(iaorb(j))),&
       & iphorb(j), cnfigfio(isa(iaorb(j)),iphorb(j)),&
       & symfio(isa(iaorb(j)),iphorb(j)), j=1,nuotot)

  write(iu) 1, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp ! ik, k(1:3), kpoint_weight - at Gamma
  write(iu) 1              ! ispin
  write(iu) no_occ_wfs     ! nwflist(ik)

  do iw = 1,no_occ_wfs
     ! coerce input array to sp:
     do j = 1,nuotot
        aux(j) = real(ks_psi(j, iw), kind=sp)
     enddo

     write(iu) iw               ! indwf
     write(iu) eo(iw,1,1)/eV    ! eo here left for compatibility
     write(iu) aux(:)           ! (aux(1:,j), j=1,nuotot)
  end do

  call io_close(iu)
  deallocate(aux)

end subroutine write_ks_psi


end module ks_flux_procs
