! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
subroutine pdosg( nspin, nuo, no, maxspn, maxnh, &
    maxo, numh, listhptr, listh, H, S, &
    E1, E2, nhist, sigma, indxuo, eo, &
    haux, saux, psi, dtot, dpr, nuotot )
  ! **********************************************************************
  ! Find the density of states projected onto the atomic orbitals
  !     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
  ! where n run over all the bands between two given energies
  ! Written by J. Junquera and E. Artacho. Nov' 99
  ! Gamma point version adapted from PDOSK by Julian Gale. Feb' 03
  ! ****  INPUT  *********************************************************
  ! integer nspin             : Number of spin components (1 or 2)
  ! integer nuo               : Number of atomic orbitals in the unit cell
  ! integer no                : Number of atomic orbitals in the supercell
  ! integer maxspn            : Second dimension of eo and qo 
  !                             (maximum number of differents spin polarizations)
  ! integer maxnh             : Maximum number of orbitals interacting
  !                             with any orbital
  ! integer maxo              : First dimension of eo
  ! integer numh(nuo)         : Number of nonzero elements of each row
  !                             of hamiltonian matrix
  ! integer listhptr(nuo)     : Pointer to each row (-1) of the
  !                             hamiltonian matrix
  ! integer listh(maxnh)      : Nonzero hamiltonian-matrix element
  !                             column indexes for each matrix row
  ! real*8  H(maxnh,nspin)    : Hamiltonian in sparse format
  ! real*8  S(maxnh)          : Overlap in sparse format
  ! real*8  E1, E2            : Energy range for density-matrix states
  !                             (to find local density of states)
  !                             Not used if e1 > e2
  ! integer nhist             : Number of the subdivisions of the histogram
  ! real*8  sigma             : Width of the gaussian to expand the eigenvectors
  ! integer indxuo(no)        : Index of equivalent orbital in unit cell
  ! real*8  eo(maxo,maxspn)   : Eigenvalues
  ! integer nuotot            : Total number of orbitals per unit cell
  ! ****  AUXILIARY  *****************************************************
  ! real*8  haux(nuo,nuo)     : Auxiliary space for the hamiltonian matrix
  ! real*8  saux(nuo,nuo)     : Auxiliary space for the overlap matrix
  ! real*8  psi(nuo,nuo)      : Auxiliary space for the eigenvectors
  ! ****  OUTPUT  ********************************************************
  ! real*8  dtot(nhist,2)   : Total density of states
  ! real*8  dpr(nhist,nuo,2): Proyected density of states
  ! **********************************************************************

  use precision
  use parallel,     only : Node, Nodes, BlockSize
  use parallelsubs, only : GetNodeOrbs, LocalToGlobalOrb
  use units,        only : pi
  use alloc,        only : re_alloc, de_alloc
#ifdef MPI
  use mpi_siesta
#endif
  use sys,          only : die

  implicit none

  integer :: nspin, nuo, no, maxspn, maxnh, maxo, nhist, nuotot

  integer :: numh(nuo), listhptr(nuo), listh(maxnh), indxuo(no)

  real(dp) :: H(maxnh,nspin), S(maxnh), E1, E2, sigma, eo(maxo,maxspn), &
      haux(nuotot,nuo), saux(nuotot,nuo), psi(nuotot,nuo), &
      dtot(nhist,2), dpr(nhist,nuotot,2) 

  ! Internal variables ---------------------------------------------------
  integer :: ispin, iuo, juo, j, jo, ihist, iband, ind, ierror

  real(dp) :: delta, ener, diff, pipj1, gauss, norm

#ifdef MPI
  integer :: BNode, Bnuo, ibandg, maxnuo, MPIerror
  real(dp), dimension(:,:), pointer ::  Sloc
  real(dp), dimension(:,:,:), pointer :: tmp
#endif

  external :: rdiag

  ! Initialize some variables
  delta = (E2 - E1)/(nhist-1)

  ! Solve eigenvalue problem for each k-point
  do ispin = 1, nspin

    call setup_Gamma()

    ! Diagonalize at the Gamma point
    call rdiag( haux, saux, nuotot, nuo, nuotot, eo(1,ispin), &
        psi, nuotot, 1, ierror, BlockSize)

    ! Check error flag and take appropriate action
    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation')
    elseif ( ierror < 0 ) then
      call setup_Gamma()
      call rdiag( haux, saux, nuotot, nuo, nuotot, eo(1,ispin), &
          psi, nuotot, 1, ierror, BlockSize)
    endif

    call setup_S()

#ifdef MPI
    ! Find maximum number of orbitals per node
    call MPI_AllReduce(nuo,maxnuo,1,MPI_integer,MPI_max, &
        MPI_Comm_World,MPIerror)

    ! Allocate workspace array for broadcast overlap matrix
    nullify( Sloc )
    call re_alloc( Sloc, 1, nuotot, 1, maxnuo, name='Sloc', routine='pdosg' )

    ! Loop over nodes broadcasting overlap matrix
    do BNode = 0,Nodes-1

      ! Find out how many orbitals there are on the broadcast node
      call GetNodeOrbs(nuotot,BNode,Nodes,Bnuo)

      ! Transfer data
      if (Node.eq.BNode) then
        Sloc(1:nuotot,1:Bnuo) = Saux(1:nuotot,1:Bnuo)
      endif
      call MPI_Bcast(Sloc(1,1),nuotot*Bnuo, &
          MPI_double_precision,BNode,MPI_Comm_World,MPIerror)

      ! Loop over all the energy range
      do ihist = 1, nhist
        ener = E1 + (ihist - 1) * delta
        do iband = 1, nuo
          call LocalToGlobalOrb(iband,Node,Nodes,ibandg)
          diff = (ener - eo(ibandG,ispin))**2 / (sigma ** 2)
          if (diff < 15.0d0) then
            gauss = exp(-diff)
            if (Node.eq.BNode) then
              ! Only add once to dtot - not everytime loop over processors is executed
              dtot(ihist,ispin) = dtot(ihist,ispin) + gauss
            endif
            do jo = 1, Bnuo
              call LocalToGlobalOrb(jo,BNode,Nodes,juo)
              do iuo = 1, nuotot
                ! Solo para los Juo que satisfagan el criterio del record...
                pipj1 = psi(iuo,iband) * psi(juo,iband)
                dpr(ihist,juo,ispin) = dpr(ihist,juo,ispin) + pipj1*gauss*Sloc(iuo,jo)
              enddo
            enddo
          endif
        enddo

      enddo

      ! End loop over broadcast nodes
    enddo

    ! Free workspace array for overlap
    call de_alloc( Sloc, name='Sloc' )

#else
    ! Loop over all the energy range
    do ihist = 1, nhist
      ener = E1 + (ihist - 1) * delta
      do iband = 1, nuo
        diff = (ener - eo(iband,ispin))**2 / (sigma ** 2)
        if (diff < 15.0d0) then
          gauss = exp(-diff)
          dtot(ihist,ispin) = dtot(ihist,ispin) + gauss
          do iuo = 1, nuotot
            do juo = 1, nuotot
              pipj1 = psi(iuo,iband) * psi(juo,iband)
              dpr(ihist,juo,ispin) = dpr(ihist,juo,ispin) + pipj1*gauss*saux(iuo,juo)
            enddo
          enddo
        endif
      enddo

    enddo
#endif

  enddo

#ifdef MPI
  ! Allocate workspace array for global reduction
  nullify( tmp )
  call re_alloc( tmp, 1, nhist, 1, max(nuotot,nspin), 
  &               1, nspin, name='tmp', routine='pdosg' )

  ! Global reduction of dpr matrix
  tmp(1:nhist,1:nuotot,1:nspin) = 0.0d0
  call MPI_AllReduce(dpr(1,1,1),tmp(1,1,1),nhist*nuotot*nspin, &
      MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
  dpr(1:nhist,1:nuotot,1:nspin) = tmp(1:nhist,1:nuotot,1:nspin)

  ! Global reduction of dtot matrix
  tmp(1:nhist,1:nspin,1) = 0.0d0
  call MPI_AllReduce(dtot(1,1),tmp(1,1,1),nhist*nspin, &
      MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
  dtot(1:nhist,1:nspin) = tmp(1:nhist,1:nspin,1)

  ! Free workspace array for global reduction
  call de_alloc( tmp, name='tmp' )
#endif

  norm = sigma * sqrt(pi)

  do ihist = 1, nhist
    do ispin = 1, nspin
      dtot(ihist,ispin) = dtot(ihist,ispin) / norm
      do iuo = 1, nuotot
        dpr(ihist,iuo,ispin) = dpr(ihist,iuo,ispin) /norm
      enddo
    enddo
  enddo

contains

  subroutine setup_Gamma()

    ! Initialize auxiliary variables
    do iuo = 1,nuo
      do juo = 1,nuotot
        Saux(juo,iuo) = 0.0d0
        Haux(juo,iuo) = 0.0d0
      enddo
    enddo

    do iuo = 1, nuo
      do j = 1, numh(iuo)
        ind = listhptr(iuo) + j
        jo  = listh(ind)
        juo = indxuo(jo)
        ! Calculate the Hamiltonian and the overlap
        Saux(juo,iuo) = Saux(1,juo,iuo) + S(ind)
        Haux(juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin)
      enddo
    enddo

  end subroutine setup_Gamma

  subroutine setup_S()

    ! Initialize auxiliary variables 
    do iuo = 1,nuo
      do juo = 1,nuotot
        Saux(juo,iuo) = 0.0d0
      end do
    end do

    do iuo = 1, nuo
      do j = 1, numh(iuo)
        ind = listhptr(iuo) + j
        jo  = listh(ind)
        juo = indxuo(jo)
        Saux(juo,iuo) = Saux(juo,iuo) + S(ind)
      enddo
    enddo

  end subroutine setup_S

end subroutine pdosg
