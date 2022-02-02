! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
SUBROUTINE PDOS3K( NUO, NO, MAXUO, MAXNH, &
    MAXO, NUMH, LISTHPTR, LISTH, H, S, &
    E1, E2, NHIST, SIGMA, &
    XIJ, INDXUO, NK, KPOINT, WK, EO, &
    HAUX, SAUX, PSI, DTOT, DPR, NUOTOT )

  ! **********************************************************************
  ! Find the density of states projected onto the atomic orbitals
  !     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
  ! where n run over all the bands between two given energies
  ! Written by J. Junquera and E. Artacho. Nov' 99
  ! Spin-orbit coupling version by J. Ferrer, October 2007
  ! ****  INPUT  *********************************************************
  ! INTEGER NUO               : Number of atomic orbitals in the unit cell
  ! INTEGER NO                : Number of atomic orbitals in the supercell
  ! INTEGER MAXUO             : Maximum number of atomic orbitals in the unit cell
  ! INTEGER MAXNH             : Maximum number of orbitals interacting
  !                             with any orbital
  ! INTEGER MAXO              : First dimension of eo
  ! INTEGER NUMH(NUO)         : Number of nonzero elements of each row
  !                             of hamiltonian matrix
  ! INTEGER LISTHPTR(NUO)     : Pointer to each row (-1) of the
  !                             hamiltonian matrix
  ! INTEGER LISTH(MAXNH)      : Nonzero hamiltonian-matrix element
  !                             column indexes for each matrix row
  ! REAL*8  H(MAXNH,4)    : Hamiltonian in sparse format
  ! REAL*8  S(MAXNH)          : Overlap in sparse format
  ! REAL*8  E1, E2            : Energy range for density-matrix states
  !                             (to find local density of states)
  !                             Not used if e1 > e2
  ! INTEGER NHIST             : Number of the subdivisions of the histogram
  ! REAL*8  SIGMA             : Width of the gaussian to expand the eigenvectors
  ! REAL*8  XIJ(3,MAXNH)      : Vectors between orbital centers (sparse)
  !                             (not used if only gamma point)
  ! INTEGER INDXUO(NO)        : Index of equivalent orbital in unit cell
  ! INTEGER NK                : Number of k points
  ! REAL*8  KPOINT(3,NK)      : k point vectors
  ! REAL*8  WK(NK)            : Weights for k points
  ! REAL*8  EO(MAXO,2,NK): Eigenvalues
  ! INTEGER NUOTOT            : Total number of orbitals per unit cell
  ! ****  AUXILIARY  *****************************************************
  ! REAL*8  HAUX(2,NUO,NUO)   : Auxiliary space for the hamiltonian matrix
  ! REAL*8  SAUX(2,NUO,NUO)   : Auxiliary space for the overlap matrix
  ! REAL*8  PSI(2,NUO,NUO)    : Auxiliary space for the eigenvectors
  ! ****  OUTPUT  ********************************************************
  ! REAL*8  DTOT(NHIST,4)   : Total density of states
  ! REAL*8  DPR(NHIST,NUO,4): Proyected density of states
  ! **********************************************************************

  use precision
  use parallel,     only : Node, Nodes, BlockSize
  use parallelsubs, only : GetNodeOrbs, LocalToGlobalOrb
  use units,        only : pi
  use alloc,        only : re_alloc, de_alloc
#ifdef MPI
  use mpi_siesta
#endif
  use sys, only : die

  IMPLICIT NONE

  INTEGER :: NUO, NO, MAXUO, MAXNH, NK, MAXO, NHIST, NUOTOT

  INTEGER :: NUMH(NUO), LISTHPTR(NUO), LISTH(MAXNH), INDXUO(NO)

  real(dp) :: H(MAXNH,8), S(MAXNH), E1, E2, SIGMA, &
      XIJ(3,MAXNH), KPOINT(3,NK), EO(MAXO*2,NK), &
      DTOT(NHIST,4), DPR(NHIST,NUOTOT,4), WK(NK)
  complex(dp), target :: psi(2,nuotot,2*nuo)

  complex(dp) :: Haux(2,nuotot,2,nuo), Saux(2,nuotot,2,nuo)
  complex(dp), pointer :: caux(:,:)

  complex(dp) kphs

  EXTERNAL :: CDIAG


  ! Internal variables ---------------------------------------------------
  INTEGER :: IK, ISPIN, IUO, JUO, IO, J, JO, IHIST, IBAND, IND, IERROR

  real(dp) :: KXIJ, CKXIJ, SKXIJ, DELTA, ENER, DIFF, GAUSS, NORM, WKSUM

  complex(dp) :: D11, D12, D22

  complex(dp), pointer :: Spr(:,:) => null()

#ifdef MPI
  integer :: BNode, Bnuo, ibandg, maxnuo, MPIerror
  complex(dp), pointer :: Sloc(:,:)
  real(dp) :: tmp(nhist,nuotot,4)
#endif

  ! START -----------------------------------------------------------------

  ! Initialize some variables
  delta = (e2 - e1)/(nhist-1)

  call re_alloc(Spr, 1, nuotot, 1, nuo, name='Spr', routine='pdos3k')

  !     Find eigenvalues at every k point
  do ik = 1,nk

    Saux = cmplx(0.0_dp,0.0_dp,dp)
    Haux = cmplx(0.0_dp,0.0_dp,dp)

    do iuo = 1,nuo
      do j = 1,numh(iuo)
        ind = listhptr(iuo) + j
        jo = listh(ind)
        juo = indxuo(jo)
        kxij = kpoint(1,ik) * xij(1,ind) + &
            kpoint(2,ik) * xij(2,ind) + &
            kpoint(3,ik) * xij(3,ind)
        kphs = exp(cmplx(0.0_dp, -kxij, dp))

        Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)   * kphs
        Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)   * kphs
        Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + cmplx(H(ind,1), H(ind,5),dp) * kphs
        Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + cmplx(H(ind,2), H(ind,6),dp) * kphs
        Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo) + cmplx(H(ind,3), - H(ind,4),dp) * kphs
        Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo) + cmplx(H(ind,7), + H(ind,8),dp) * kphs

      enddo
    enddo

    ! Diagonalize for each k point. Note duplication of problem size
    call cdiag( Haux, Saux, 2*nuotot, 2*nuo, 2*nuotot, &
        eo(:,ik), psi, 2*nuotot, 1, ierror, 2*BlockSize )
    if (ierror.gt.0) then
      call die('Terminating due to failed diagonalisation')
    elseif (ierror.lt.0) then
      ! Repeat diagonalisation with increased memory to handle clustering
      Saux = cmplx(0.0_dp,0.0_dp,dp)
      Haux = cmplx(0.0_dp,0.0_dp,dp)

      do iuo = 1,nuo
        do j = 1,numh(iuo)
          ind = listhptr(iuo) + j
          jo = listh(ind)
          juo = indxuo(jo)
          kxij = kpoint(1,ik) * xij(1,ind) + &
              kpoint(2,ik) * xij(2,ind) + &
              kpoint(3,ik) * xij(3,ind)
          kphs = exp(cmplx(0.0_dp, -kxij, dp))

          Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)   * kphs
          Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)   * kphs
          Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + cmplx(H(ind,1), H(ind,5),dp) * kphs
          Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + cmplx(H(ind,2), H(ind,6),dp) * kphs
          Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo) + cmplx(H(ind,3), - H(ind,4),dp) * kphs
          Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo) + cmplx(H(ind,7), + H(ind,8),dp) * kphs

        enddo
      enddo

      call cdiag( Haux, Saux, 2*nuotot, 2*nuo, 2*nuotot, &
          eo(:,ik), psi, 2*nuotot, 1, ierror, 2*BlockSize )
    endif

    ! Recalculate again the overlap matrix in k-space

    Spr = cmplx(0.0_dp,0.0_dp,dp)
    do iuo = 1,nuo
      do j = 1,numh(iuo)
        ind = listhptr(iuo) + j
        jo = listh(ind)
        juo = indxuo(jo)
        kxij = kpoint(1,ik) * xij(1,ind) + &
            kpoint(2,ik) * xij(2,ind) + &
            kpoint(3,ik) * xij(3,ind)
        ! Since we are doing element wise multiplications (and not dot-products)
        ! we might as well setup the transpose S(k)^T == S(-k) because this will
        ! mean that we can do a simpler multiplication further down
        kphs = exp(cmplx(0.0_dp, kxij, dp))
        Spr(juo,iuo) = Spr(juo,iuo) + S(ind) * kphs
      enddo
    enddo

#ifdef MPI
    ! Find maximum number of orbitals per node
    call MPI_AllReduce(nuo,maxnuo,1,MPI_integer,MPI_max, MPI_Comm_World,MPIerror)

    ! Allocate workspace array for broadcast overlap matrix
    nullify( Sloc )
    call re_alloc( Sloc, 1, nuotot, 1, maxnuo, name='Sloc', routine='pdos3k' )

    ! Loop over nodes broadcasting overlap matrix
    do BNode = 0,Nodes-1

      ! Find out how many orbitals there are on the broadcast node
      call GetNodeOrbs(nuotot,BNode,Nodes,Bnuo)

      ! Transfer data
      if (Node.eq.BNode) then
        Sloc(1:nuotot,1:Bnuo) = Spr(1:nuotot,1:Bnuo)
      endif
      call MPI_Bcast(Sloc(1,1),nuotot*Bnuo, MPI_double_complex,BNode,MPI_Comm_World,MPIerror)

      ! Loop over all the energy range

      do ihist = 1, nhist
        ener = E1 + (ihist - 1) * delta
        do iband = 1, nuo*2
          call LocalToGlobalOrb((iband+1)/2,Node,Nodes,ibandg)
          ibandg = ibandg * 2 - mod(iband, 2) 
          diff = (ener - eo(ibandg,ik))**2 / (sigma ** 2)
          if (diff .gt. 15.0d0) cycle
          gauss = exp(-diff) * wk(ik)
          caux => psi(:,:,iband) ! c_{up,j}, c_{down,j}
          do jo = 1, Bnuo
            call LocalToGlobalOrb(jo,BNode,Nodes,juo)
            do io = 1, nuotot
              D11 = caux(1,io) * conjg(caux(1,juo)) * Sloc(io,jo)
              D22 = caux(2,io) * conjg(caux(2,juo)) * Sloc(io,jo)
              D12 = caux(1,io) * conjg(caux(2,juo)) * Sloc(io,jo)
              !D21 = caux(2,io) * conjg(caux(1,juo)) * Sloc(io,jo)

              D11 = gauss*D11
              D22 = gauss*D22
              D12 = gauss*D12

              dpr(ihist,juo,1) = dpr(ihist,juo,1) + real(D11,dp)
              dpr(ihist,juo,2) = dpr(ihist,juo,2) + real(D22,dp)
              dpr(ihist,juo,3) = dpr(ihist,juo,3) + real(D12,dp)
              dpr(ihist,juo,4) = dpr(ihist,juo,4) - aimag(D12)
            enddo
          enddo
        enddo
      enddo

    enddo !BNode

    ! Free workspace array for overlap
    call de_alloc(Sloc, 'Sloc', 'pdos3k')

#else
    ! Loop over all the energy range

    do ihist = 1, nhist
      ener = E1 + (ihist - 1) * delta
      do iband = 1, nuo*2
        diff = (ener - eo(iband,ik))**2 / (sigma ** 2)
        if (diff .gt. 15.0d0) cycle
        gauss = exp(-diff) * wk(ik)
        caux => psi(:,:,iband) ! c_{up,j}, c_{down,j}
        do io = 1, nuotot
          do jo = 1, nuotot
            D11 = caux(1,io) * conjg(caux(1,jo)) * Spr(io,jo)
            D22 = caux(2,io) * conjg(caux(2,jo)) * Spr(io,jo)
            D12 = caux(1,io) * conjg(caux(2,jo)) * Spr(io,jo)
            !D21 = caux(2,io) * conjg(caux(1,jo)) * Spr(io,jo)

            D11 = gauss*D11
            D22 = gauss*D22
            D12 = gauss*D12

            dpr(ihist,jo,1) = dpr(ihist,jo,1) + real(D11,dp)
            dpr(ihist,jo,2) = dpr(ihist,jo,2) + real(D22,dp)
            dpr(ihist,jo,3) = dpr(ihist,jo,3) + real(D12,dp)
            dpr(ihist,jo,4) = dpr(ihist,jo,4) - aimag(D12)
          enddo
        enddo
      enddo
    enddo

#endif

  enddo  ! nk

  call de_alloc(Spr, 'Spr', 'pdos3k')

#ifdef MPI

  ! Global reduction of DPR matrix
  tmp = 0.0d0
  call MPI_AllReduce(dpr(1,1,1),tmp(1,1,1),nhist*nuotot*4, &
      MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
  dpr = tmp

#endif

  wksum = 0.0d0
  do ik = 1, nk
    wksum = wksum + wk(ik)
  enddo
  norm = sigma * sqrt(pi) * wksum
  dpr = dpr /norm

  do ihist = 1, nhist
    dtot(ihist,1) = sum(dpr(ihist,:,1))
    dtot(ihist,2) = sum(dpr(ihist,:,2))
    dtot(ihist,3) = sum(dpr(ihist,:,3))
    dtot(ihist,4) = sum(dpr(ihist,:,4))
  enddo

end SUBROUTINE PDOS3K
