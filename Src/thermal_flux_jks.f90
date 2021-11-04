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
    use sparse_matrices, only: maxnh

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

    real(dp), allocatable :: gradS(:,:)

    allocate (tmp_g(no_l))
    allocate(work(no_u), ipiv(no_u), amat(no_u,no_u))
    allocate(gradS(maxnh, 3))

    ! computing and storage of every of the 3 components
    ! should be available to order in the .fdf

    psi_hat_c(:,:,:) = 0.0_dp        ! Init result to zeros outside main loop

    call momentum_operat(gradS)

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
                ! tmp_g(nu) = tmp_g(nu) + psi_base(col,iw) * gradS_base(idx,k)

                write(2001,*) k, gradS_base(:,k)
                write(2002,*) k, gradS(k,:)

                tmp_g(nu) = tmp_g(nu) - psi_base(col,iw) * gradS(k,idx)

                ! if ((col.le.n_wfs).and.(col.ne.nu)) then
                !    tmp_g(nu) = tmp_g(nu) + psi_base(col,iw) * gradS_base(idx,k)
                        ! & / (eo_base(col) - eo_base(nu))
                ! end if
             end do
          end do

          ! do alpha = 1,no_u
          !    psi_hat_c(alpha,iw,idx) = tmp_g(alpha)
          ! end do

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

  subroutine momentum_operat(gradS)
    use atomlist,        only: no_u, no_s, no_l, indxuo, iaorb, iphorb, lasto, qtot, amass, lastkb, iphKB
    use siesta_geom
    use densematrix,     only: psi ! To get at the c^i_\mu coefficients of the wavefunctions
    use sparse_matrices, only: listh, listhptr, numh

    use precision,    only : dp
    use parallel,     only : Node, Nodes
    use atmparams,    only : lmx2, nzetmx, nsemx
    use atmfuncs,     only : epskb, lofio, mofio, rcut, rphiatm
    use atmfuncs,     only : orb_gindex, kbproj_gindex
    use atm_types,    only : nspecies
    use parallelsubs, only : GlobalToLocalOrb, LocalToGlobalOrb
    use alloc,        only : re_alloc, de_alloc
    use sys,          only : die
    use neighbour,    only : jna=>jan, xij, r2ij
    use neighbour,    only : mneighb, reset_neighbour_arrays
    use matel_mod,    only : new_matel

    integer :: ia, iio, ind, io, ioa, is, ix
    integer :: j, ja, jn, jo, joa, js, nnia, norb, ka

    real(dp) :: grSij(3), rij, Sij, xinv(3), sum

    integer :: ikb, ina, ino, jno, ko, koa, ks, ig, jg, kg
    integer ::  nkb, nna, nno, ilm1, ilm2, npoints, ir

    real(dp) ::  epsk, grSki(3), rki, rmax, rmaxkb, rmaxo
    real(dp) ::  Sik, Sjk, Sikr2, Sikr3, Sikr4, Sjkr2, Sjkr3, Sjkr4
    real(dp) ::  dintg2(3), dint1, dint2, dintgmod2
    real(dp) ::  dintg1(3), dintgmod1
    real(dp) ::  phi1, phi2, dphi1dr, dphi2dr, Sir0(3), r

    integer,  pointer, save :: iano(:)
    integer,  pointer, save :: iono(:)
    integer,           save :: maxkba = 25
    integer,           save :: maxno = 1000
!N      logical,  pointer, save :: calculated(:,:,:)
    logical,  pointer, save :: listed(:)
    logical,  pointer, save :: needed(:)
    logical                 :: within
    real(dp),          save :: dx = 0.01d0
!N      real(dp), pointer, save :: Pij(:,:,:)
    real(dp), pointer, save :: Si(:,:)
    real(dp), pointer, save :: Ski(:,:,:)
    real(dp),          save :: tiny = 1.0d-9
    real(dp), pointer, save :: Vi(:,:)

    integer   ::  na, no, nua, nuo

    real(dp), allocatable, intent(inout) :: gradS(:,:)

    gradS(:,:) = 0.0_dp

    na =  na_s
    nua = na_u
    no =  no_s
    nuo = no_u

! Allocate arrays
      norb = lmx2*nzetmx*nsemx
      call re_alloc( listed, 1, no, 'listed', 'temp_jks' )
      call re_alloc( needed, 1, no, 'needed', 'temp_jks' )
      call re_alloc( Si,     1, no, 1, 3, 'Si',     'temp_jks' )
      call re_alloc( Vi,     1, no, 1, 3, 'Vi',     'temp_jks' )
      call re_alloc( iano, 1, maxno, 'iano',  'temp_jks' )
      call re_alloc( iono, 1, maxno, 'iono',  'temp_jks' )
      call re_alloc( Ski, 1, 4, 1, maxkba, 1, maxno, 'Ski', 'temp_jks' )

      !------

! Initialize arrayd Vi only once
        no = lasto(na)
        do jo = 1,no
          Vi(jo,:) = 0.0d0
          listed(jo) = .false.
          needed(jo) = .false.
        enddo

! Find out which orbitals are needed on this node
        do iio = 1,nuo
          call LocalToGlobalOrb(iio,Node,Nodes,io)
          needed(io) = .true.
          do j = 1,numh(iio)
            ind = listhptr(iio) + j
            jo = listh(ind)
            needed(jo) = .true.
          enddo
        enddo

! Loop on atoms with KB projectors
        do ka = 1,na
          ks = isa(ka)
          nkb = lastkb(ka) - lastkb(ka-1)
          if (nkb.gt.maxkba) then
            maxkba = nkb + 10
            call re_alloc( Ski, 1, 4, 1, maxkba, 1, maxno, copy=.true., &
     &                     name='Ski', routine='temp_jks' )
          endif

! Find neighbour atoms
          call mneighb( scell, rmax, na, xa, ka, 0, nna )

! Find neighbour orbitals
          nno = 0
          do ina = 1,nna
            ia = jna(ina)

            is = isa(ia)
            rki = sqrt(r2ij(ina))
            do io = lasto(ia-1)+1,lasto(ia)
              if (needed(io)) then
                call GlobalToLocalOrb(io,Node,Nodes,iio)
                ioa = iphorb(io)
                ig = orb_gindex(is,ioa)

! Find if orbital is within range
                within = .false.
                do ko = lastkb(ka-1)+1,lastkb(ka)
                  koa = iphKB(ko)
                  if ( rki .lt. rcut(is,ioa)+rcut(ks,koa) ) &
     &              within = .true.
                enddo

! Find overlap between neighbour orbitals and KB projectors
                if (within) then
                  nno = nno + 1
                  if (nno.gt.maxno) then
                    maxno = nno + 500
                    call re_alloc( iano, 1, maxno, copy=.true., &
     &                             name='iano', routine='temp_jks' )
                    call re_alloc( iono, 1, maxno, copy=.true., &
     &                             name='iono', routine='temp_jks' )
                    call re_alloc( Ski, 1, 2, 1, maxkba, 1, maxno, &
     &                             copy=.true., name='Ski', &
     &                             routine='temp_jks' )
                  endif
                  iono(nno) = io
                  iano(nno) = ia
                  ikb = 0
                  do ko = lastkb(ka-1)+1,lastkb(ka)
                    ikb = ikb + 1
                    koa = iphKB(ko)
                    kg = kbproj_gindex(ks,koa)
                    do ix = 1,3
                     xinv(ix) = - xij(ix,ina)
                    enddo
                    call new_MATEL('S', ig, kg, xinv, Ski(1,ikb,nno), grSki)

                    call new_MATEL('X', ig, kg, xinv, Ski(2,ikb,nno), grSki)
                    call new_MATEL('Y', ig, kg, xinv, Ski(3,ikb,nno) , grSki)
                    call new_MATEL('Z', ig, kg, xinv, Ski(4,ikb,nno), grSki)
                  enddo
                endif
              endif
            enddo

          enddo

! Loop on neighbour orbitals
          do ino = 1,nno
            io = iono(ino)
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio .gt. 0) then
              ia = iano(ino)
              if (ia .le. nua) then

! Scatter filter of desired matrix elements
                do j = 1,numh(iio)
                  ind = listhptr(iio) + j
                  jo = listh(ind)
                  listed(jo) = .true.
                enddo

! Find matrix elements with other neighbour orbitals
                do jno = 1,nno
                  jo = iono(jno)
                  if (listed(jo)) then

! Loop on KB projectors
                    ikb = 0
                    do ko = lastkb(ka-1)+1,lastkb(ka)
                      ikb = ikb + 1
                      koa = iphKB(ko)
                      epsk = epskb(ks,koa)
                      Sik = Ski(1,ikb,ino)
                      ! Sikr= Ski(2,ikb,ino)
                      Sikr2= Ski(2,ikb,ino)
                      Sikr3= Ski(3,ikb,ino)
                      Sikr4= Ski(4,ikb,ino)
                      Sjk = Ski(1,ikb,jno)
                      ! Sjkr= Ski(2,ikb,jno)
                      Sjkr2= Ski(2,ikb,jno)
                      Sjkr3= Ski(3,ikb,jno)
                      Sjkr4= Ski(4,ikb,jno)
                      ! Vi(jo) = Vi(jo) + epsk * (Sik * Sjkr - Sikr * Sjk)
                      Vi(jo,1) = Vi(jo,1) + epsk * (Sik * Sjkr2 - Sikr2 * Sjk)
                      Vi(jo,2) = Vi(jo,2) + epsk * (Sik * Sjkr3 - Sikr3 * Sjk)
                      Vi(jo,3) = Vi(jo,3) + epsk * (Sik * Sjkr4 - Sikr4 * Sjk)
                    enddo

                  endif
                enddo

! Pick up contributions to H and restore Di and Vi
                do j = 1,numh(iio)
                  ind = listhptr(iio) + j
                  jo = listh(ind)
                  ! Sp(ind) = Sp(ind) + Vi(jo)
                  gradS(ind,1) = gradS(ind,1) + Vi(jo,1)
                  gradS(ind,2) = gradS(ind,2) + Vi(jo,2)
                  gradS(ind,3) = gradS(ind,3) + Vi(jo,3)
                  Vi(jo,:) = 0.0d0
                  listed(jo) = .false.
                enddo

              endif
            endif
          enddo
        enddo

      !------

    ! Initialize neighb subroutine
    call mneighb( scell, 2.0d0*rmaxo, na, xa, 0, 0, nnia )

      do jo = 1,no
         Si(jo,:) = 0.0d0
      enddo

      do ia = 1,nua

        is = isa(ia)
        call mneighb( scell, 2.0d0*rmaxo, na, xa, ia, 0, nnia )

        do io = lasto(ia-1)+1,lasto(ia)
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio .gt. 0) then
            ioa = iphorb(io)
            ig = orb_gindex(is,ioa)
            do jn = 1,nnia
              do ix = 1,3
                xinv(ix) = - xij(ix,jn)
              enddo
              ja = jna(jn)
              rij = sqrt( r2ij(jn) )
              do jo = lasto(ja-1)+1,lasto(ja)
                joa = iphorb(jo)
                js = isa(ja)
                jg = orb_gindex(js,joa)

                if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then

                    if (rij.lt.tiny) then
! Perform the direct computation of the matrix element of the momentum
! within the same atom
!N                     if ( .not.calculated(joa,ioa,is) ) then
                       ilm1 = lofio(is,ioa)**2 + lofio(is,ioa) + &
     &                      mofio(is,ioa) + 1
                       ilm2 = lofio(is,joa)**2 + lofio(is,joa) + &
     &                      mofio(is,joa) + 1
                       call intgry(ilm1,ilm2,dintg2)
                       call intyyr(ilm1,ilm2,dintg1)
                       dintgmod1 = dintg1(1)**2 + dintg1(2)**2 + &
     &                      dintg1(3)**2
                       dintgmod2 = dintg2(1)**2 + dintg2(2)**2 + &
     &                      dintg2(3)**2
                       Sir0(:) = 0.0d0
                       if ((dintgmod2.gt.tiny).or.(dintgmod1.gt.tiny)) &
     &                      then
                          dint1 = 0.0d0
                          dint2 = 0.0d0
                          npoints = int(max(rcut(is,ioa),rcut(is,joa)) &
     &                         /dx) + 2
                          do ir = 1,npoints
                             r = dx*(ir-1)
                             call rphiatm(is,ioa,r,phi1,dphi1dr)
                             call rphiatm(is,joa,r,phi2,dphi2dr)
                             dint1 = dint1 + dx*phi1*dphi2dr*r**2
                             dint2 = dint2 + dx*phi1*phi2*r
                          enddo
!     The factor of two because we use Ry for the Hamiltonian
                          Sir0(1) = dint1*dintg1(1)+dint2*dintg2(1)
     ! &                  -2.0d0*(dk(1)*(dint1*dintg1(1)+dint2*dintg2(1))+ &
     ! &                      dk(2)*(dint1*dintg1(2)+dint2*dintg2(2))+ &
     ! &                      dk(3)*(dint1*dintg1(3)+dint2*dintg2(3))) &
                          Sir0(2) = dint1*dintg1(2)+dint2*dintg2(2)
                          Sir0(3) = dint1*dintg1(3)+dint2*dintg2(3)
                       endif
!N     Pij(ioa,joa,is) = - Sir0
!N     Pij(joa,ioa,is) =   Sir0
                       Si(jo,:) = Sir0(:)
!N                       calculated(ioa,joa,is) = .true.
!N                       calculated(joa,ioa,is) = .true.
!N                    endif

                    else
! Matrix elements between different atoms are taken from the
! gradient of the overlap
                      call new_MATEL('S', ig, jg, xij(1:3,jn), &
     &                           Sij, grSij )
! The factor of two because we use Ry for the Hamiltonian
                      Si(jo,1) = grSij(1)
                      Si(jo,2) = grSij(2)
                      Si(jo,3) = grSij(3)
     ! &                  2.0d0*(grSij(1)*dk(1) &
     ! &             +           grSij(2)*dk(2) &
     ! &             +           grSij(3)*dk(3))
                  endif
                endif
              enddo
            enddo
            do j = 1,numh(iio)
              ind = listhptr(iio) + j
              jo = listh(ind)
              gradS(ind,1) = gradS(ind,1) + Si(jo,1)
              gradS(ind,2) = gradS(ind,2) + Si(jo,2)
              gradS(ind,3) = gradS(ind,3) + Si(jo,3)
              Si(jo,:) = 0.0d0
            enddo
          endif
        enddo
      enddo

      call reset_neighbour_arrays( )
  end subroutine momentum_operat

end module thermal_flux_jks
