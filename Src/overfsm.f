! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_overfsm

      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut, orb_gindex
      use neighbour,     only : jna=>jan, r2ij, xij, mneighb,
     &                          reset_neighbour_arrays
      use alloc,         only : re_alloc, de_alloc
      use m_new_matel,   only : new_matel

      implicit none

      public :: overfsm
      private

      CONTAINS

      subroutine overfsm(nua, na, no, scell, xa, indxua, rmaxo,
     .                   maxnh, maxnd, lasto, iphorb, isa, 
     .                   numd, listdptr, listd, numh, listhptr, listh, 
     .                   nspin, Escf, fa, stress, S, DS, D2S)
C *********************************************************************
C Overlap matrix and orthonormalization contribution to forces and stress.
C Energies in Ry. Lengths in Bohr.
C Writen by J.Soler and P.Ordejon, August-October'96, June'98
C Added linres features by L.Riches and S.Illera (2016) 
C **************************** INPUT **********************************
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer no               : Number of orbitals in supercell
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
c integer indxua(na)       : Index of equivalent atom in unit cell
C real*8  rmaxo            : Maximum cutoff for atomic orbitals
C integer maxnh            : First dimension of S and listh
C integer maxnd            : First dimension of Escf and listd
C integer lasto(0:na)      : Last orbital index of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom
C integer isa(na)          : Species index of each atom
C integer numd(nuotot)     : Number of nonzero elements of each row
C                            of Escf
C integer listdptr(nuotot) : Pointer to start of rows (-1) of Escf
C integer listd(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of Escf
C integer numh(nuotot)     : Number of nonzero elements of each row
C                            of the overlap matrix
C integer listhptr(nuotot) : Pointer to start of rows (-1) of overlap
C                            matrix
C integer listh(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of the overlap matrix
C integer nspin            : Number of spin components of Escf
C integer Escf(maxnd,nspin): Energy-Density matrix
C **************************** OUTPUT *********************************
C real*8  S(maxnh)         : Sparse overlap matrix
C ********************** INPUT and OUTPUT *****************************
C real*8  fa(3,nua)        : Atomic forces (Orthog. part added to input)
C real*8  stress(3,3)      : Stress tensor (Orthog. part added to input)
C Linres variables-----------------------------------------------------
C real*8 dS(maxnh,3)       : First derivative of the overlap matrix 
C real*8 d2S(maxnh,3,3)    : Second derivative of the overlap matrix
C *********************************************************************

      integer, intent(in) ::
     . maxnd, maxnh, na, no, nspin, nua

      integer, intent(in) ::
     . indxua(na), iphorb(no), isa(na), lasto(0:na), 
     . listd(*), numd(*), listh(maxnh), numh(*), listdptr(*),
     . listhptr(*)

      real(dp) , intent(inout)  :: fa(3,nua), stress(3,3)
      real(dp) , intent(in)     :: scell(3,3), Escf(maxnd,nspin), 
     .                 rmaxo, xa(3,na)
      real(dp), intent(out)     :: S(maxnh)
      real(dp), optional, intent(out) :: DS(:,:), D2S(:,:,:) !Linres line

C Internal variables ......................................................
  
      integer
     .  ia, ind, io, ioa, is, ispin, ix, iio, ig, jg,
     .  j, ja, jn, jo, joa, js, jua, jx, nnia

      real(dp)
     .  fij(3), grSij(3) , rij, Sij, volcel, volume

      real(dp), dimension(:), pointer  ::  Di, Si

C LINRES ------------------------------------------------------------
       real(dp) :: gr2Sij(3,3)
       real(dp), dimension(:,:), pointer ::  dSi
       real(dp), dimension(:,:,:), pointer :: d2Si
C -------------------------------------------------------------------

      external  timer
C ......................

C Start timer
      call timer( 'overfsm', 1 )

      volume = nua * volcel(scell) / na

C Initialize neighb subroutine 
      call mneighb( scell, 2.d0*rmaxo, na, xa, 0, 0, nnia)

C Allocate local memory
      nullify( Di )
      call re_alloc( Di, 1, no, 'Di', 'overfsm' )
      nullify( Si )
      call re_alloc( Si, 1, no, 'Si', 'overfsm' )

      Si(1:no) = 0.0d0  ! AG: Superfluous
      Di(1:no) = 0.0d0

C LINRES ------------------------------------------------------------
       nullify( dSi )
       call re_alloc( dSi, 1, no, 1, 3, 'dSi', 'overfsm' )
       nullify( d2Si )
       call re_alloc( d2Si, 1, no, 1, 3, 1, 3, 'd2Si', 'overfsm' )
C -------------------------------------------------------------------

      do ia = 1,nua
        is = isa(ia)
        call mneighb( scell, 2.d0*rmaxo, na, xa, ia, 0, nnia )
        do io = lasto(ia-1)+1,lasto(ia)

C Is this orbital on this Node?
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio.gt.0) then

C Valid orbital 
            ioa = iphorb(io)
            ig = orb_gindex(is,ioa)
            do j = 1,numd(iio)
              ind = listdptr(iio)+j
              jo = listd(ind)
              do ispin = 1,nspin  
                Di(jo) = Di(jo) + Escf(ind,ispin)
              enddo
            enddo
            do jn = 1,nnia
              ja = jna(jn)
              jua = indxua(ja)
              rij = sqrt( r2ij(jn) )
              do jo = lasto(ja-1)+1,lasto(ja)
                joa = iphorb(jo)
                js = isa(ja)
                if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then
                  jg = orb_gindex(js,joa)
                  if (present(d2S)) then
                     call new_MATEL( 'S', ig, jg, xij(1:3,jn),
     .                            Sij, grSij , gr2Sij)
                  else
                     call new_MATEL( 'S', ig, jg, xij(1:3,jn),
     .                            Sij, grSij)
                  endif
                  Si(jo) = Si(jo) + Sij
                  do ix = 1,3
                    if ((.not. present (dS)) .or. (.not.
     &                             present (d2S))) then !Lin line
                      fij(ix) = (- Di(jo)) * grSij(ix)
                      fa(ix,ia)  = fa(ix,ia)  + fij(ix)
                      fa(ix,jua) = fa(ix,jua) - fij(ix)
                    endif
C LINRES ------------------------------------------------------------
                    if (present(DS)) then
                     dSi(jo,ix) = dSi(jo,ix) + grSij(ix)
                    endif
C -------------------------------------------------------------------
                    do jx = 1,3
                      if ((.not. present (dS)) .or. (.not.
     &                             present (d2S))) then !Lin line
                        stress(jx,ix) = stress(jx,ix) +
     .                              xij(jx,jn) * fij(ix) / volume
                      endif
C LINRES ------------------------------------------------------------
                      if (present(D2S)) then
                         d2Si(jo,ix,jx) = d2Si(jo,ix,jx) + gr2Sij(ix,jx)
                      endif
C ------------------------------------------------------------------
                    enddo
                  enddo
                endif
              enddo
            enddo
            do j = 1,numd(iio)
              jo = listd(listdptr(iio)+j)
              Di(jo) = 0.0d0
            enddo
            do j = 1,numh(iio)
              ind = listhptr(iio)+j
              jo = listh(ind)
              S(ind) = Si(jo)
              Si(jo) = 0.0d0
C LINRES ------------------------------------------------------------
              if (present(DS)) then
                do ix = 1,3
                  DS(ind,ix) = dSi(jo,ix)
                    dSi(jo,ix) = 0.0d0

                   do jx = 1,3
                    D2S(ind,ix,jx) =  d2Si(jo,ix,jx)
                    d2Si(jo,ix,jx) = 0.0d0
                   enddo
                enddo
               !if(ia.eq.2) print *, iio,';',jo,';',D2S(ind,1,1)
              endif
C -------------------------------------------------------------------
            enddo
          endif
        enddo
      enddo

C Deallocate local memory
!      call new_MATEL( 'S', 0, 0, 0, 0, xij, Sij, grSij )
      call reset_neighbour_arrays( )
      call de_alloc( Si, 'Si', 'overfsm' )
      call de_alloc( Di, 'Di', 'overfsm' )
C LINRES ------------------------------------------------------------
      if (present(dS)) then 
        call de_alloc( dSi, 'dSi', 'overfsm' )
        call de_alloc( d2Si, 'd2Si', 'overfsm' )
      endif
C -------------------------------------------------------------------



C Finish timer
      call timer( 'overfsm', 2 )

      end subroutine overfsm
      end module m_overfsm




