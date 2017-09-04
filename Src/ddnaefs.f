! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_ddnaefs
      public :: ddnaefs
      CONTAINS
      subroutine ddnaefs( nua, na, scell, xa, indxua, rmaxv,
     .                   isa, dynmat)
C *********************************************************************
C Correction of Neutral Atom energies, forces and stress due to the
C overlap between ionic (bare pseudopotential) charges.
C Energies in Ry. Lengths in Bohr.
C Written by J.Soler and P.Ordejon, August-October'96, June'98 
C Modified by DSP Aug., 1998
C Adapted to compute dynmat terms (7) for LINRES, LR 2015
C **************************** INPUT **********************************
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
c integer indxua(na)       : Index of equivalent atom in unit cell
C real*8  rmaxv            : Maximum cutoff for NA potential
C integer isa(na)          : Species index of each atom
C **************************** OUTPUT *********************************
C ********************** INPUT and OUTPUT *****************************
C dynmat(3,nua,3,nua) 	   : Dynamical matrix
C *********************************************************************
C The following function must exits:
C
C    INTEGER  FUNCTION  IZOFIS(IS) : Returns the atomic number
C Input:
C      INTEGER IS : Specie index
C
C *********************************************************************

      use precision
      use atmfuncs,  only: izofis, psover, vna_gindex
      use neighbour, only: mneighb, jna=>jan, xij, r2ij,
     &                     reset_neighbour_arrays
      use m_new_matel,   only : new_matel

      implicit none

      integer  na, nua

      integer  indxua(na), isa(na)

      real(dp)
     .  scell(3,3), rmaxv, xa(3,na), 
     .  dynmat(3,nua,3,nua)


C Internal variables ......................................................
      integer
     .  ia, is, ix, ja, jn, js, jx, jua, nnia, i, j, ig, jg

      real(dp)
     .  dvdr, fij(3), rij, r2min, vij, volcel, volume, d2vdr,
     .  gr2ij(3,3), pi
       
      parameter ( r2min = 1.d-15 )
C ......................
      call mneighb( scell, 2.0d0*rmaxv, na, xa, 0, 0, nnia )

      volume = nua * volcel(scell) / na

      pi = 4.0_dp * atan(1.0_dp)

 
      do ia = 1,nua

C Find neighbour atoms
        call mneighb( scell, 2.0d0*rmaxv, na, xa, ia, 0, nnia )
        do jn = 1,nnia
          ja = jna(jn)
          jua = indxua(ja)
          is = isa(ia)
          js = isa(ja) 

          if (r2ij(jn).gt.r2min .and.
     .        izofis(is).gt.0 .and. izofis(js).gt.0) then

            ig=vna_gindex(is)
            jg=vna_gindex(js)

            call new_matel( 'T', ig, jg, xij(:,jn), vij, fij, gr2ij ) 

            rij = sqrt( r2ij(jn) )
            call psover( is, js, rij, vij, dvdr, d2vdr )
            do ix = 1,3
             do jx = 1,3
              gr2ij(ix,jx) = gr2ij(ix,jx) / (16.0_dp*pi)
              if (r2ij(jn).gt.r2min) then
                dynmat(jx,jua,ix,ia) = dynmat(jx,jua,ix,ia)
     .          + d2vdr*xij(ix,jn)*xij(jx,jn) / (rij*rij) / 2.0_dp
     .          - dvdr*xij(ix,jn)*xij(jx,jn) / (rij**3) /2.0_dp
                dynmat(jx,ia,ix,jua) = dynmat(jx,ia,ix,jua)
     .          + d2vdr*xij(ix,jn)*xij(jx,jn) / (rij*rij) / 2.0_dp
     .          - dvdr*xij(ix,jn)*xij(jx,jn) / (rij**3) /2.0_dp
                dynmat(jx,ia,ix,ia) = dynmat(jx,ia,ix,ia)-
     .            d2vdr*xij(ix,jn)*xij(jx,jn) / (rij*rij) / 2.0_dp
     .          + dvdr*xij(ix,jn)*xij(jx,jn) / (rij**3) /2.0_dp
                dynmat(jx,jua,ix,jua) = dynmat(jx,jua,ix,jua)-
     .            d2vdr*xij(ix,jn)*xij(jx,jn) / (rij*rij) / 2.0_dp
     .          + dvdr*xij(ix,jn)*xij(jx,jn) / (rij**3) /2.0_dp
                if (ix.eq.jx) then
                  dynmat(jx,jua,ix,ia) = dynmat(jx,jua,ix,ia) +
     .               dvdr / rij /2.0_dp
                  dynmat(jx,ia,ix,jua) = dynmat(jx,ia,ix,jua) +
     .               dvdr / rij /2.0_dp
                  dynmat(jx,ia,ix,ia) = dynmat(jx,ia,ix,ia) -
     .               dvdr / rij /2.0_dp
                  dynmat(jx,jua,ix,jua) = dynmat(jx,jua,ix,jua) -
     .               dvdr / rij /2.0_dp
                endif
              endif
              dynmat(jx,jua,ix,ia)=dynmat(jx,jua,ix,ia)+gr2ij(ix,jx)
              dynmat(jx,ia,ix,jua)=dynmat(jx,ia,ix,jua)+gr2ij(ix,jx)
              dynmat(jx,ia,ix,ia)=dynmat(jx,ia,ix,ia)-gr2ij(ix,jx)
              dynmat(jx,jua,ix,jua)=dynmat(jx,jua,ix,jua)-gr2ij(ix,jx)
             enddo
            enddo
          endif
        enddo
      enddo
      call reset_neighbour_arrays( )
      end subroutine ddnaefs
      end module m_ddnaefs
