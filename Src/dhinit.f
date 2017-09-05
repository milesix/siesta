        subroutine dhinit(ialr, na, maxnh, maxnd, nspin,lasto, 
     &             listh, listhptr, numh, dS, dSmat, dHmat, dT)

C *************************   LINEAR RESPONSE   ************************
C Initializes first order changed overlap matrix elements and
C constructs first order change in non-self-consistent hamiltonian 
C coming from kinetic part.
C Junquera, Ordejon and Pruneda 
C SIESTA implementation L. Riches, March '15
C Cheked and commented S. Illera, April '16
***************************INPUT*************************************
C INTEGER IALR		      : index of displaced atom
C INTEGER NA		      : number of atoms in the supercell
C INTEGER MAXNH		      : Maximum number of orbitals interacting
C INTEGER MAXND		      : Maximum number of nonzero elements of 
C                               each row of density matrix
C INTEGER NSPIN		      : number of different spin polarizations
C INTEGER LASTO		      : position of last orbital of each atom
C INTEGER LISTH(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C INTEGER LISTHPTR(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C INTEGER NUMH		      : Number of nonzero elements of each row 
C                               of density matrix
C REAL*8 DS		      : change in the overlap matrix elements as 
C				given by SIESTA
C REAL*8 DT		      : change in the kinetic elements 
****************************OUTPUT**********************************
C REAL*8 DSMAT		      : overlap changed matrix elements for moving IALR
C REAL*8 DHMAT		      : perturbed hamiltonian 
*******************************************************************

C Modules-----------------------------------------------------------------
!        use m_iodmlr,      only: read_dmlr
        use precision,     only: dp
        use atomlist,      only: no_l, indxuo, indxua, iaorb
        use alloc,         only: re_alloc, de_alloc
        use parallel,      only: Node, Nodes
        use parallelsubs,  only: GlobalToLocalOrb
        use listsc_module, only: listsc
        use alloc,         only: re_alloc
        use mesh,          only: iatfold, nmsc, cmesh
        use siesta_geom,   only: xa
C -------------------------------------------------------------------------

        implicit none

C Internal variable types and dimensions
        integer      :: nspin_read, maxnd_read, j, jo, 
     &                  io, iio, ialr ,na, lasto(0:na),  ind, ix,
     &                  listh(maxnh), listhptr(*), numh(*), juo,
     &                   ind2, k, ko, ispin
        logical      :: found

        integer      :: maxnh, maxnd, nspin, i, ia, iua, iu, ju
        real(dp)     :: DS(maxnh,3), DSMAT(maxnh,3), 
     &                  DHMAT(maxnh,3, nspin), DT(maxnh,3) 
        real(dp)     :: displaat(3), dist(3), qxij, qpoint(3), cqxij, pi
        
      call timer( 'dHinit', 1 )  
C Read change in the density matrix from file ----------------------- 
!     call read_dmlr ( maxnh, no_l, nspin, numh,
!    &                     listhptr, listh, dDscf, found )

      qpoint=(/0.0,0.0,0.0/)
      pi = 4*atan(1.0_dp)

      dSmat = 0.0_dp
C Loop on orbitals of atom IALR ----------------------------------------

      do io = lasto(IALR-1)+1, lasto(IALR) !orbitals from moving atom
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        if (iio .gt. 0) then !local orbital 
          do j = 1,numh(iio) !interacting orbitals
            ind = listhptr(iio)+j !position to io-j interaction 
            jo = listh(ind)
            juo = indxuo(jo)
            do ix = 1,3 ! calculate mu,nu elements (mu is in the unit cell)
              DSMAT(ind,ix) = DSMAT(ind,ix) - DS(ind,ix)
              do ispin = 1,nspin
                DHMAT(ind,ix,ispin) = DHMAT(ind,ix,ispin)
     .                               - DT(ind,ix)
              enddo
            enddo
C         Idetify neighbour index of orbital io relative to orbital jo -
            do k = 1, numh(juo)
              ind2 = listhptr(juo)+ k
              ko = listsc(jo,juo,listh(ind2)) 
              if ( ko .EQ. iio) then !add contributions nu,mu elements
                ia  = iaorb(listh(ind2))
                iua = indxua(ia)
                do ix = 1, 3
                  displaat(ix) = 
     &                 (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+ 
     &                 (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+ 
     &                 (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
                end do
                dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
                qxij    = qpoint(1) *  dist(1)  + 
     &                    qpoint(2) *  dist(2)  + 
     &                    qpoint(3) *  dist(3)
                cqxij = cos(qxij)

                do ix = 1,3
                  DSMAT(ind2,ix) = DSMAT(ind2,ix) -
     &                           DS(ind,ix)
                  DSMAT(ind2,ix) = cqxij * DSMAT(ind2,ix)
                  do ispin = 1,nspin
                    DHMAT(ind2,ix,ispin) = DHMAT(ind2,ix,ispin)
     &                                 - DT(ind,ix)
                  enddo
                enddo
              endif
            enddo
          enddo
        endif
      enddo

      call timer( 'dHinit', 2 )

      end subroutine
