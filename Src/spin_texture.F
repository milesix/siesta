!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine spin_texture( nuo, no, maxo, maxnh, numh, 
     .                   listhptr, listh, H, S,  
     .                   xij, indxuo, nbk, maxbk, kpt, ef, 
     .                   nuotot, iscf, neigwanted)

!
      use precision
      use sys
      use parallel,      only : Node, Nodes, BlockSize
      use parallelsubs,  only : LocalToGlobalOrb
      use m_spin, only: spin

#ifdef MPI
      use mpi_siesta
#endif
      implicit           none
C *********************************************************************
C Calculates the spin texture for selected k-points.  
C Adapted from diag3k by R.Robles, January 2018
C **************************** INPUT **********************************
C integer nuo                 : Number of basis orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer maxuo               : Maximum number of basis orbitals
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxo                : First dimension of eo and qo
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxnh,4)          : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer nbk                 : Number of k points
C real*8  kpt(3,nbk)          : k point vectors
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C integer iscf                : SCF cycle number
C integer neigwanted          : Number of eigenvalues wanted
C *************************** AUXILIARY *******************************
C complex*16 Haux(2,nuotot,2,nuo): Auxiliary space for the hamiltonian matrix
C complex*16 Saux(2,nuotot,2,nuo): Auxiliary space for the overlap matrix
C complex*16 psi(2,nuotot,2*nuo) : Auxiliary space for the eigenvectors
C complex*16 Dke(2,nuotot,nuo,2*nuotot) : Aux. space that may be the same as Haux
C complex*16 caux(2,nuotot)      : Extra auxiliary space
C *************************** UNITS ***********************************
C xij and kpt must be in reciprocal coordinates of each other.
C *************************** PARALLEL ********************************
C The auxiliary arrays are now no longer symmetry and so the order
C of referencing has been changed in several places to reflect this.
C *********************************************************************
!
!     INPUT / OUTPUT
      
      integer maxuo, maxnh, maxo, maxbk
      integer nbk, no, nuo, no_u, nuotot, iscf
      
      integer indxuo(no), numh(nuo)
      integer listh(maxnh)
      integer listhptr(nuo)
      
      real(dp) H(maxnh,8), S(maxnh)
      real(dp) st(4), st_f(4)
      real(dp) kpt(3,nbk), wk(nbk), xij(3,maxnh)
      real(dp) ef, eo(maxo*2,nbk)
      integer :: neigwanted
      
      complex(dp), dimension(2,nuotot,2,nuo) :: Haux, Saux
      complex(dp), dimension(2,nuotot,2,nuo,maxo*2) :: Dke
      complex(dp), dimension(2,nuotot,2*nuo) :: psi
      complex(dp), dimension(2,nuotot) :: caux
      
!     TEMPOS, INTERNAL VARIABLES etc.
      
      integer  BNode, BTest, ie, ierror, iie, ik, ind, io, iio
      integer  iuo, j, jo, juo, neigneeded, Node0, iounit
      real(dp) eV

      parameter ( eV = 1.d0 / 13.60580d0 )
      
!     Haux(js,juo,is,iuo) = <js,juo|H|is,iuo>
!     Indices is and js are for spin components
!     Indices iuo and juo are for orbital components
      complex(dp)                                 :: cicj
      complex(dp)                                 :: D11, D22, D12, D21
      complex(dp)                                 :: kphs
      real(dp)                                    :: kxij
#ifdef MPI
      integer            :: MPIerror
#endif

      external              cdiag
!***********************************************************************
!     B E G I N
!***********************************************************************

C Start time counter
      call timer( 'spin_texture', 1 )

C Check parameter maxbk
      if (nbk .gt. maxbk) then
        if (Node.eq.0) then
          write(6,'(/,a,/,a)')
     .       'spin_texture: WARNING: parameter maxbk too small',
     .       'spin_texture: No spin_texture calculation performed'
        endif
        goto 999
      endif

!     Find spin texture at every k point

      Node0 = Node
      iounit = 700 + Node

! Spin texture is written in file spin_texture.dat

      open(unit=iounit,file='spin_texture.dat',status='unknown')
      if(Node0.eq.0)
     .write(iounit,'(A,i4,A,i4,A,f12.6)') "k-points = ", nbk, 
     .     ", eigenvalues = ", 2*neigwanted, ", E-fermi (eV) = ", ef/eV
      do ik = 1,nbk
      
        Saux = dcmplx(0.0_dp,0.0_dp)
        Haux = dcmplx(0.0_dp,0.0_dp)
      
        do iuo = 1,nuo
          do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            kxij = kpt(1,ik) * xij(1,ind) +
     .             kpt(2,ik) * xij(2,ind) +
     .             kpt(3,ik) * xij(3,ind)

            kphs = cdexp(dcmplx(0.0_dp,-1.0_dp)*kxij)
            Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)*kphs
            Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)*kphs
            Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) +
     .                          dcmplx(H(ind,1), H(ind,5))*kphs
            Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) +
     .                          dcmplx(H(ind,2), H(ind,6))*kphs
            Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo) +
     .                          dcmplx(H(ind,3),-H(ind,4))*kphs
            Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo) +
     .                          dcmplx(H(ind,7), H(ind,8))*kphs
          enddo
        enddo
      
!       Find eigenvalues
!       Note duplication since the spins are mixed

        call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,eo(1,ik),psi,
     .             -2*neigwanted,iscf,ierror, 2*BlockSize)

        if (ierror.ne.0) then
          call die('Terminating due to failed diagonalisation')
        endif
      enddo
      
      do ik = 1,nbk
        if(Node0.eq.0) then
          write(iounit,'(/A,i4,A,3f12.6,A)') 
     .      "k-point ", ik, " (",(kpt(j,ik),j=1,3)," ) (Bohr^-1)"
          write(iounit,'(A)') 
     .      "     ie     e(eV)      Sx      Sy      Sz  "
        endif

        neigneeded = 2*neigwanted
      
!       Find eigenvectors 
        Saux = dcmplx(0.0_dp,0.0_dp)
        Haux = dcmplx(0.0_dp,0.0_dp)
        do iuo = 1,nuo
          do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            kxij = kpt(1,ik) * xij(1,ind) +
     .             kpt(2,ik) * xij(2,ind) +
     .             kpt(3,ik) * xij(3,ind)

            kphs = cdexp(dcmplx(0.0_dp, -1.0_dp)*kxij)
            Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)   * kphs
            Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)   * kphs
            Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + 
     .                            dcmplx(H(ind,1), H(ind,5)) * kphs
            Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + 
     .                            dcmplx(H(ind,2), H(ind,6)) * kphs
            Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo)
     .                        + dcmplx(H(ind,3), - H(ind,4)) * kphs
            Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo)
     .                        + dcmplx(H(ind,7), + H(ind,8)) * kphs
          enddo
        enddo

        caux(:,:) = dcmplx(0.0_dp,0.0_dp)
        call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,caux,psi,
     .             neigneeded,iscf,ierror, 2*BlockSize)
      
!       Check error flag and take appropriate action
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.lt.0) then
!         Repeat diagonalisation with increased memory to handle clustering
          Saux = dcmplx(0.0_dp,0.0_dp)
          Haux = dcmplx(0.0_dp,0.0_dp)
          do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              kxij = kpt(1,ik) * xij(1,ind) +
     .               kpt(2,ik) * xij(2,ind) +
     .               kpt(3,ik) * xij(3,ind)

              kphs = cdexp(dcmplx(0.0_dp, -1.0_dp)*kxij)
              Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind)   * kphs
              Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind)   * kphs
              Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + 
     .                              dcmplx(H(ind,1), H(ind,5)) * kphs
              Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + 
     .                              dcmplx(H(ind,2), H(ind,6)) * kphs
              Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo)
     .                          + dcmplx(H(ind,3), - H(ind,4)) * kphs
              Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo)
     .                          + dcmplx(H(ind,7), + H(ind,8)) * kphs
            enddo
          enddo
          call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,caux,psi,
     .               neigneeded,iscf,ierror, 2*BlockSize)
        endif
      
! Loop over energies to calculate Dke

        Dke  = dcmplx(0.0_dp,0.0_dp)
!
        BNode = 0
! Node = 0?
        Node = Node0
        iie = 0
        do ie = 1,2*nuotot
          if (Node.eq.BNode) then
            iie = iie + 1
          endif
      
          caux(:,:) = dcmplx(0.0_dp,0.0_dp)
            if (Node.eq.BNode) then
              do j = 1,nuotot
                caux(1,j)= psi(1,j,iie) ! c_{i,up}
                caux(2,j)= psi(2,j,iie) ! c_{i,dn}
              enddo
            endif
#ifdef MPI
            call MPI_Bcast(caux(1,1),2*nuotot,MPI_double_complex,BNode,
     .           MPI_Comm_World,MPIerror)
#endif
            do iuo = 1,nuo
              call LocalToGlobalOrb(iuo,Node,Nodes,iio)
              do juo = 1,nuotot
      
!------- 1,1 -----------------------------------------------------------
               cicj = caux(1,iio) * dconjg(caux(1,juo))
               Dke(1,juo,1,iuo,ie) = Dke(1,juo,1,iuo,ie) + cicj
!------- 2,2 -----------------------------------------------------------
               cicj = caux(2,iio) * dconjg(caux(2,juo))
               Dke(2,juo,2,iuo,ie) = Dke(2,juo,2,iuo,ie) + cicj
!------- 1,2 -----------------------------------------------------------
               cicj = caux(1,iio) * dconjg(caux(2,juo))
               Dke(1,juo,2,iuo,ie) = Dke(1,juo,2,iuo,ie) + cicj
!------- 2,1 -----------------------------------------------------------
               cicj = caux(2,iio) * dconjg(caux(1,juo))
               Dke(2,juo,1,iuo,ie) = Dke(2,juo,1,iuo,ie) + cicj
              enddo
            enddo
!        
          BTest = ie/(2*BlockSize)
          if (BTest*2*BlockSize.eq.ie) then
            BNode = BNode + 1
            if (BNode .gt. Nodes-1) BNode = 0
          endif
        enddo

! Loop in energy to calculate spin texture (st).

       do ie = 1,2*nuotot
        st(:) = 0.0_dp
        do iuo = 1,nuo
        do j = 1,numh(iuo)
          ind = listhptr(iuo) + j
          jo = listh(ind)
          juo = indxuo(jo)
          kxij = kpt(1,ik) * xij(1,ind) +
     .           kpt(2,ik) * xij(2,ind) +
     .           kpt(3,ik) * xij(3,ind)
          kphs = cdexp(dcmplx(0.0_dp,-1.0_dp)*kxij)
      
          D11 = Dke(1,juo,1,iuo,ie) * kphs
          D22 = Dke(2,juo,2,iuo,ie) * kphs
          D12 = Dke(1,juo,2,iuo,ie) * kphs
          D21 = Dke(2,juo,1,iuo,ie) * kphs
      
          st(1) = st(1) + S(ind) * (dreal(D11) + dreal(D22))
          st(2) = st(2) + S(ind) * 2.0d0 * dreal(D12)
          st(3) = st(3) - S(ind) * 2.0d0 * dimag(D12)
          st(4) = st(4) + S(ind) * (dreal(D11) - dreal(D22))
        enddo
        enddo
#ifdef MPI
        call MPI_Barrier(MPI_Comm_World, MPIError) ! Wait for everybody
        call MPI_Reduce(st(1),st_f(1),4,
     .    MPI_double_precision,MPI_sum,0,MPI_Comm_World,MPIerror)
#endif
        if(Node0.eq.0) then
          write(iounit,'(i7,f12.5,3f8.4)') 
     .      ie,eo(ie,ik)/eV,(st_f(j),j=2,4)
        endif
       enddo
      enddo

C This is the only exit point
  999 continue
      call timer( 'spin_texture', 2 )
      
!***********************************************************************      
      return
      end subroutine spin_texture
!***********************************************************************      
      
